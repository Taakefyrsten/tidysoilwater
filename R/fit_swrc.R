#' Fit Van Genuchten parameters to observed retention data
#'
#' Fits the Van Genuchten (1980) soil water retention model to observed (h, θ)
#' data. Returns a tidy tibble of estimated parameters. If the input data is
#' grouped (via [dplyr::group_by()]), the fit is performed independently for
#' each group.
#'
#' ## Fitting paths
#'
#' The function selects an optimisation back-end based on the arguments
#' supplied:
#'
#' | Arguments | Back-end |
#' |-----------|----------|
#' | Any combination without `fixed` | `nls(algorithm = "port")` — box-constrained NLS with physical default bounds; SE from NLS Jacobian |
#' | `fixed` specified (+ optional bounds + `weights`) | `stats::optim(method = "L-BFGS-B")` — reduced-dimension optimisation; SE from numerical Hessian (approximate) |
#'
#' ## Starting values
#'
#' Auto-initialised from the data:
#' - `theta_r`: 90% of minimum observed θ
#' - `theta_s`: 105% of maximum observed θ
#' - `alpha`: 0.1
#' - `n`: 2.0
#'
#' @param data A data frame or tibble, optionally grouped with
#'   [dplyr::group_by()].
#' @param theta_col Bare column name of observed volumetric water content
#'   (m³/m³).
#' @param h_col Bare column name of observed matric potential / pressure head.
#'   Absolute values are used internally.
#' @param workers Number of parallel workers. Defaults to `1` (sequential).
#'   Values > 1 use [parallel::mclapply()] on Unix. On Windows, silently
#'   reduced to `1`.
#' @param lower Named numeric vector of lower bounds for parameters. Names
#'   should be a subset of `c("theta_r", "theta_s", "alpha", "n")`. Unspecified
#'   parameters receive physical defaults (0, 1e-6, 1e-6, 1.001 respectively).
#' @param upper Named numeric vector of upper bounds. Same naming convention as
#'   `lower`. Unspecified parameters default to (1, 1, Inf, Inf).
#' @param fixed Named numeric vector of parameters to hold at fixed values
#'   rather than estimate. E.g. `fixed = c(theta_r = 0.05)`. Fixed parameters
#'   are excluded from the optimisation and reported with `NA` standard errors.
#' @param weights Bare column name or numeric scalar; per-observation weights
#'   for the least-squares objective. Larger weights increase influence of that
#'   observation. Defaults to uniform weights.
#'
#' @return An object of class `fit_swrc` (a tibble subclass) with one row per
#'   group containing:
#'   - Group keys (if input was grouped)
#'   - `theta_r`, `theta_s`, `alpha`, `n`: fitted parameter values (or the
#'     supplied value for fixed parameters)
#'   - `std_error_theta_r`, `std_error_theta_s`, `std_error_alpha`,
#'     `std_error_n`: standard errors (`NA` for fixed parameters)
#'   - `convergence`: `TRUE` if optimisation converged
#'
#'   The result stores the original data and fitting options as attributes,
#'   enabling [confint.fit_swrc()] to compute profile-likelihood confidence
#'   intervals.
#'
#' @seealso [confint.fit_swrc()], [fit_swrc_hcc()], [swrc_van_genuchten()]
#'
#' @references
#' Van Genuchten, M. Th. (1980). A closed-form equation for predicting the
#' hydraulic conductivity of unsaturated soils. *Soil Science Society of
#' America Journal*, 44(5), 892–898.
#' <https://doi.org/10.2136/sssaj1980.03615995004400050002x>
#'
#' @examples
#' library(tibble)
#'
#' observed <- tibble(
#'   h     = c(0, 10, 100, 500, 1000, 5000, 15000),
#'   theta = c(0.44, 0.40, 0.32, 0.25, 0.20, 0.12, 0.08)
#' )
#'
#' # Default unconstrained fit
#' fit_swrc(observed, theta_col = theta, h_col = h)
#'
#' # Box-constrained fit (n must be between 1.1 and 4)
#' fit_swrc(observed, theta_col = theta, h_col = h,
#'          lower = c(n = 1.1), upper = c(n = 4))
#'
#' # Fix theta_r at a known value; estimate the other three parameters
#' fit_swrc(observed, theta_col = theta, h_col = h,
#'          fixed = c(theta_r = 0.05))
#'
#' # Grouped fit
#' library(dplyr)
#' tibble(
#'   horizon = rep(c("A", "B"), each = 7),
#'   h       = rep(c(0, 10, 100, 500, 1000, 5000, 15000), 2),
#'   theta   = c(0.44, 0.40, 0.32, 0.25, 0.20, 0.12, 0.08,
#'               0.38, 0.33, 0.27, 0.21, 0.16, 0.10, 0.07)
#' ) |>
#'   group_by(horizon) |>
#'   fit_swrc(theta_col = theta, h_col = h)
#'
#' @export
fit_swrc <- function(data, theta_col, h_col, workers = 1L,
                     lower = NULL, upper = NULL, fixed = NULL,
                     weights = NULL) {
  theta_quo   <- rlang::enquo(theta_col)
  h_quo       <- rlang::enquo(h_col)
  weights_quo <- rlang::enquo(weights)

  theta_nm    <- rlang::as_name(theta_quo)
  h_nm        <- rlang::as_name(h_quo)
  has_weights <- !rlang::quo_is_null(weights_quo)

  for (nm in c(theta_nm, h_nm)) {
    if (!nm %in% names(data)) {
      cli::cli_abort(c(
        "Column {.val {nm}} not found in {.arg data}.",
        "i" = "Available columns: {.val {names(data)}}."
      ))
    }
  }

  workers <- as.integer(workers)
  if (workers < 1L) workers <- 1L
  if (workers > 1L && .Platform$OS.type == "windows") {
    cli::cli_warn(c("Parallel fitting not supported on Windows.",
                    "i" = "Falling back to {.code workers = 1}."))
    workers <- 1L
  }

  # Parameter names and bound resolution
  all_nms   <- c("theta_r", "theta_s", "alpha", "n")
  fixed_nms <- names(fixed)
  free_nms  <- setdiff(all_nms, fixed_nms)

  lower_full <- .vg_default_lower
  upper_full <- .vg_default_upper
  if (!is.null(lower)) lower_full[names(lower)] <- lower
  if (!is.null(upper)) upper_full[names(upper)] <- upper

  has_fixed  <- length(fixed_nms) > 0
  has_bounds <- !is.null(lower) || !is.null(upper)

  # ── Core fitter for one group ──────────────────────────────────────────────

  fit_one_group <- function(df, ...) {
    theta_vec <- df[[theta_nm]]
    h_vec     <- abs(df[[h_nm]])

    w_vec <- if (has_weights) {
      ww <- resolve_arg(weights_quo, df, "weights")
      if (length(ww) == 1L) rep(ww, length(theta_vec)) else ww
    } else rep(1.0, length(theta_vec))

    start_all <- c(
      theta_r = min(theta_vec, na.rm = TRUE) * 0.9,
      theta_s = max(theta_vec, na.rm = TRUE) * 1.05,
      alpha   = 0.1,
      n       = 2.0
    )
    if (has_fixed) start_all[fixed_nms] <- unlist(fixed)

    if (!has_fixed) {
      # ── Path 1: nls(algorithm = "port") — physical bounds prevent divergence
      # Uses PORT for all cases (with or without user-supplied bounds/weights).
      # The plain Gauss-Newton algorithm can diverge to unphysical theta_r < 0;
      # PORT with the default lower bound of 0 prevents this robustly.
      fit_df       <- data.frame(df)
      fit_df$h_abs <- h_vec
      vg_formula   <- stats::as.formula(
        paste0(theta_nm, " ~ theta_r + (theta_s - theta_r)",
               " / (1 + (alpha * h_abs)^n)^(1 - 1/n)")
      )
      tryCatch({
        fit <- stats::nls(
          formula   = vg_formula,
          data      = fit_df,
          start     = as.list(start_all),
          weights   = if (has_weights) w_vec else NULL,
          algorithm = "port",
          lower     = lower_full[all_nms],
          upper     = upper_full[all_nms],
          # warnOnly: return fit even for marginal convergence (code 8 = machine
          # precision; parameters are valid)
          control   = stats::nls.control(maxiter = 200, warnOnly = TRUE)
        )
        tidy_fit <- broom::tidy(fit)
        params   <- stats::setNames(tidy_fit$estimate, tidy_fit$term)
        se       <- stats::setNames(tidy_fit$std.error,  tidy_fit$term)
        tibble::tibble(
          theta_r = params[["theta_r"]], theta_s = params[["theta_s"]],
          alpha   = params[["alpha"]],   n       = params[["n"]],
          std_error_theta_r = se[["theta_r"]], std_error_theta_s = se[["theta_s"]],
          std_error_alpha   = se[["alpha"]],   std_error_n       = se[["n"]],
          convergence = TRUE
        )
      }, error = function(e) {
        cli::cli_warn(c("VG fit did not converge.",
                        "i" = "Original error: {conditionMessage(e)}"))
        .swrc_na_row()
      })

    } else {
      # ── Path 3: optim() — fixed params + optional bounds + weights ────────
      start_free  <- start_all[free_nms]
      lower_free  <- lower_full[free_nms]
      upper_free  <- upper_full[free_nms]
      fixed_list  <- as.list(setNames(unlist(fixed), fixed_nms))

      obj <- function(par) {
        params <- c(as.list(setNames(par, free_nms)), fixed_list)
        pred   <- .vg_predict(h_vec, params$theta_r, params$theta_s,
                              params$alpha, params$n)
        sum(w_vec * (theta_vec - pred)^2, na.rm = TRUE)
      }

      tryCatch({
        # nlminb (PORT) is more robust than L-BFGS-B for this problem:
        # handles infinite upper bounds and near-flat gradients without
        # ABNORMAL_TERMINATION errors.
        opt    <- stats::nlminb(start = start_free, objective = obj,
                                lower = lower_free, upper = upper_free,
                                control = list(iter.max = 500, eval.max = 1000))
        n_obs  <- sum(!is.na(theta_vec))
        n_free <- length(free_nms)
        sigma2   <- opt$objective / max(n_obs - n_free, 1L)
        hessian  <- tryCatch(stats::optimHess(opt$par, obj),
                             error = function(e) matrix(NA_real_, n_free, n_free))
        hess_inv <- tryCatch(solve(hessian),
                             error = function(e) matrix(NA_real_, n_free, n_free))
        se_free  <- sqrt(abs(diag(hess_inv)) * sigma2)

        params_all <- c(as.list(setNames(opt$par, free_nms)), fixed_list)
        se_all     <- stats::setNames(rep(NA_real_, 4L), all_nms)
        se_all[free_nms] <- se_free

        tibble::tibble(
          theta_r = params_all$theta_r, theta_s = params_all$theta_s,
          alpha   = params_all$alpha,   n       = params_all$n,
          std_error_theta_r = se_all[["theta_r"]],
          std_error_theta_s = se_all[["theta_s"]],
          std_error_alpha   = se_all[["alpha"]],
          std_error_n       = se_all[["n"]],
          convergence = (opt$convergence == 0L)
        )
      }, error = function(e) {
        cli::cli_warn(c("VG fit (fixed params) did not converge.",
                        "i" = "Original error: {conditionMessage(e)}"))
        .swrc_na_row()
      })
    }
  }

  # ── Dispatch ───────────────────────────────────────────────────────────────
  is_grouped <- dplyr::is_grouped_df(data)
  group_dfs  <- if (is_grouped) dplyr::group_split(data) else list(data)

  if (!is_grouped) {
    result <- fit_one_group(data)
  } else if (workers == 1L) {
    fits   <- lapply(group_dfs, fit_one_group)
    result <- dplyr::bind_cols(dplyr::group_keys(data), dplyr::bind_rows(fits))
  } else {
    fits   <- parallel::mclapply(group_dfs, fit_one_group, mc.cores = workers)
    result <- dplyr::bind_cols(dplyr::group_keys(data), dplyr::bind_rows(fits))
  }

  # ── Attach metadata for confint() ─────────────────────────────────────────
  class(result)           <- c("fit_swrc", class(result))
  attr(result, ".data")   <- group_dfs
  attr(result, ".nms")    <- list(theta = theta_nm, h = h_nm)
  attr(result, ".opts")   <- list(lower = lower_full, upper = upper_full,
                                  fixed = fixed, free  = free_nms)
  result
}
