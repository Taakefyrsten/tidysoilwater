#' Profile-likelihood confidence intervals for fitted VG parameters
#'
#' Computes confidence intervals for the Van Genuchten parameters returned by
#' [fit_swrc()] using the profile likelihood ratio (PLR) approach. For each
#' free parameter, one dimension of the parameter space is scanned while the
#' remaining free parameters are re-optimised at each step. The interval
#' endpoints are the values at which the profile F-statistic exceeds the
#' critical value:
#'
#' \deqn{F(p) = \frac{RSS(p) - RSS_0}{RSS_0 / (n - k)} \leq F_{1,\, n-k}(\text{level})}
#'
#' where RSS₀ is the minimum (fitted) residual sum of squares, n is the
#' number of observations, and k is the number of free parameters.
#'
#' PLR intervals are generally more reliable than normal-approximation intervals
#' (i.e., estimate ± z · SE) for the non-linear VG model, especially when
#' sample sizes are small or parameter distributions are skewed.
#'
#' @param object An object of class `fit_swrc`, as returned by [fit_swrc()].
#' @param parm Character vector of parameter names for which to compute
#'   intervals. Defaults to all free parameters (those not fixed in the
#'   original call to [fit_swrc()]).
#' @param level Confidence level. Defaults to `0.95`.
#' @param ... Not used; present for S3 compatibility.
#'
#' @return A tibble with columns:
#'   - Group keys (if the original fit was grouped)
#'   - `param` — parameter name
#'   - `estimate` — point estimate (MLE)
#'   - `lower` — lower confidence bound
#'   - `upper` — upper confidence bound
#'   - `level` — confidence level
#'
#' @seealso [fit_swrc()], [fit_swrc_hcc()]
#'
#' @examples
#' library(tibble)
#'
#' observed <- tibble(
#'   h     = c(0, 10, 100, 500, 1000, 5000, 15000),
#'   theta = c(0.44, 0.40, 0.32, 0.25, 0.20, 0.12, 0.08)
#' )
#'
#' fit <- fit_swrc(observed, theta_col = theta, h_col = h)
#' confint(fit)
#'
#' # 90% intervals for alpha and n only
#' confint(fit, parm = c("alpha", "n"), level = 0.90)
#'
#' # Grouped: returns CIs for each group
#' library(dplyr)
#' tibble(
#'   horizon = rep(c("A", "B"), each = 7),
#'   h       = rep(c(0, 10, 100, 500, 1000, 5000, 15000), 2),
#'   theta   = c(0.44, 0.40, 0.32, 0.25, 0.20, 0.12, 0.08,
#'               0.38, 0.33, 0.27, 0.21, 0.16, 0.10, 0.07)
#' ) |>
#'   group_by(horizon) |>
#'   fit_swrc(theta_col = theta, h_col = h) |>
#'   confint()
#'
#' @export
confint.fit_swrc <- function(object, parm = NULL, level = 0.95, ...) {
  data_list <- attr(object, ".data")
  nms       <- attr(object, ".nms")
  opts      <- attr(object, ".opts")

  if (is.null(data_list) || is.null(nms) || is.null(opts)) {
    cli::cli_abort(c(
      "{.arg object} does not contain fitting metadata.",
      "i" = "Re-fit using {.fn fit_swrc} to generate a valid {.cls fit_swrc} object."
    ))
  }

  theta_nm <- nms$theta
  h_nm     <- nms$h
  free_nms <- opts$free
  fixed    <- if (is.null(opts$fixed)) list() else as.list(setNames(unlist(opts$fixed), names(opts$fixed)))

  if (is.null(parm)) parm <- free_nms
  parm <- intersect(parm, free_nms)   # silently skip fixed params

  if (length(parm) == 0L) {
    cli::cli_warn("No free parameters to compute confidence intervals for.")
    return(tibble::tibble(param = character(), estimate = numeric(),
                          lower = numeric(), upper = numeric(), level = numeric()))
  }

  n_groups  <- length(data_list)
  # Group key columns = all columns except the 9 result columns
  n_result_cols <- 9L
  key_cols  <- if (ncol(object) > n_result_cols) {
    names(object)[seq_len(ncol(object) - n_result_cols)]
  } else character(0)

  ci_groups <- lapply(seq_len(n_groups), function(i) {
    df      <- data_list[[i]]
    theta   <- df[[theta_nm]]
    h_abs   <- abs(df[[h_nm]])
    n_obs   <- sum(!is.na(theta))
    n_free  <- length(free_nms)

    mle <- c(
      theta_r = object$theta_r[i], theta_s = object$theta_s[i],
      alpha   = object$alpha[i],   n       = object$n[i]
    )

    # RSS at MLE — handle weighted case if .weights stored (future: use it here)
    RSS_0 <- .vg_rss(mle, h_abs, theta, rep(1.0, length(theta)))

    if (!is.finite(RSS_0) || RSS_0 <= 0) {
      cli::cli_warn("Cannot compute profile CIs: RSS at MLE is non-positive or non-finite.")
      return(.ci_na_rows(parm, level))
    }

    F_crit <- stats::qf(level, df1 = 1L, df2 = n_obs - n_free)

    # Profile RSS: fix param p at value v, optimise remaining free params
    profile_rss <- function(p, v) {
      other_free  <- setdiff(free_nms, p)
      fixed_p     <- c(fixed, setNames(list(v), p))
      start_other <- mle[other_free]
      lower_other <- opts$lower[other_free]
      upper_other <- opts$upper[other_free]

      obj <- function(par) {
        params <- c(as.list(setNames(par, other_free)), fixed_p)
        .vg_rss(unlist(params)[c("theta_r","theta_s","alpha","n")],
                h_abs, theta, rep(1.0, length(theta)))
      }

      tryCatch(
        stats::optim(par = start_other, fn = obj, method = "L-BFGS-B",
                     lower = lower_other, upper = upper_other)$value,
        error = function(e) Inf
      )
    }

    # F-statistic as function of parameter value
    f_profile <- function(p, v) {
      (profile_rss(p, v) - RSS_0) * (n_obs - n_free) / RSS_0
    }

    # Find CI bounds for one parameter
    ci_one <- function(p) {
      v_mle <- mle[[p]]
      f0    <- f_profile(p, v_mle)    # should be ~0

      lo_bound <- opts$lower[[p]]
      hi_bound <- opts$upper[[p]]
      if (!is.finite(hi_bound)) hi_bound <- v_mle * 10 + 1

      # Lower: find where F first exceeds F_crit below MLE
      lo <- tryCatch({
        f_at_lo <- f_profile(p, lo_bound + .Machine$double.eps)
        if (f_at_lo < F_crit) {
          lo_bound  # flat profile at lower boundary → use boundary
        } else {
          stats::uniroot(
            function(v) f_profile(p, v) - F_crit,
            interval = c(lo_bound + .Machine$double.eps, v_mle),
            tol = .Machine$double.eps^0.5
          )$root
        }
      }, error = function(e) NA_real_)

      # Upper: find where F first exceeds F_crit above MLE
      hi <- tryCatch({
        f_at_hi <- f_profile(p, hi_bound)
        if (f_at_hi < F_crit) {
          hi_bound  # flat profile at upper boundary → use boundary
        } else {
          stats::uniroot(
            function(v) f_profile(p, v) - F_crit,
            interval = c(v_mle, hi_bound),
            tol = .Machine$double.eps^0.5
          )$root
        }
      }, error = function(e) NA_real_)

      tibble::tibble(param = p, estimate = v_mle, lower = lo, upper = hi,
                     level = level)
    }

    ci_df <- dplyr::bind_rows(lapply(parm, ci_one))

    # Prepend group keys
    if (length(key_cols) > 0L) {
      keys  <- object[i, key_cols, drop = FALSE]
      ci_df <- dplyr::bind_cols(
        keys[rep(1L, nrow(ci_df)), , drop = FALSE],
        ci_df
      )
    }
    ci_df
  })

  dplyr::bind_rows(ci_groups)
}

# Helper: NA CI rows for a parameter set
.ci_na_rows <- function(parm, level) {
  tibble::tibble(
    param    = parm,
    estimate = NA_real_,
    lower    = NA_real_,
    upper    = NA_real_,
    level    = level
  )
}
