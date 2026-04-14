#' Jointly fit Van Genuchten SWRC and hydraulic conductivity parameters
#'
#' Fits the Van Genuchten (1980) soil water retention curve and
#' Mualem-Van Genuchten hydraulic conductivity function simultaneously to
#' observed (h, θ) and (h, K) data. This joint fit shares the common
#' structural parameters (θ_r, θ_s, α, n) between both models, typically
#' improving identifiability compared to fitting them separately.
#'
#' The objective function minimises the combined weighted residual:
#' \deqn{RSS = w_\theta \sum ({\theta_{obs}} - \theta_{pred})^2
#'             + w_K \sum (\log_{10}K_{obs} - \log_{10}K_{pred})^2}
#'
#' K residuals are on the log₁₀ scale because K spans many orders of magnitude;
#' this prevents high-K observations from dominating the fit. The `wrc_weight`
#' and `hcc_weight` scalars allow the user to balance the two components.
#'
#' Missing values (`NA`) in `theta_col` or `K_col` are silently dropped from
#' the respective component, so WRC and HCC measurements need not be at the
#' same h values and can be in the same data frame with NAs where data is
#' absent.
#'
#' ## Parameters estimated
#'
#' | Parameter | Meaning |
#' |-----------|---------|
#' | `theta_r` | Residual water content (m³/m³) |
#' | `theta_s` | Saturated water content (m³/m³) |
#' | `alpha` | Van Genuchten shape (1/cm or 1/m) |
#' | `n` | Van Genuchten shape (dimensionless; > 1) |
#' | `Ks` | Saturated hydraulic conductivity (units of K_col) |
#' | `tau` | Mualem tortuosity (dimensionless; fixed at 0.5 by default) |
#'
#' @param data A data frame or tibble, optionally grouped with
#'   [dplyr::group_by()].
#' @param theta_col Bare column name of observed volumetric water content
#'   (m³/m³). Rows where this is `NA` are excluded from the WRC component.
#' @param K_col Bare column name of observed hydraulic conductivity. Rows
#'   where this is `NA` are excluded from the HCC component. Must be positive;
#'   the log₁₀ transformation is applied internally.
#' @param h_col Bare column name of matric potential / pressure head.
#'   Absolute values are used internally.
#' @param lower Named numeric vector of lower parameter bounds. Parameter
#'   names: `"theta_r"`, `"theta_s"`, `"alpha"`, `"n"`, `"Ks"`, `"tau"`.
#' @param upper Named numeric vector of upper parameter bounds.
#' @param fixed Named numeric vector of parameters to hold fixed. E.g.
#'   `fixed = c(tau = 0.5)` (the default behaviour if `tau` is not supplied).
#' @param wrc_weight Positive scalar weighting the WRC residual component in
#'   the joint objective. Defaults to `1`.
#' @param hcc_weight Positive scalar weighting the HCC (log K) residual
#'   component. Defaults to `1`.
#' @param workers Number of parallel workers. Defaults to `1`.
#'
#' @return An object of class `fit_swrc_hcc` (a tibble subclass) with one row
#'   per group containing:
#'   - Group keys (if grouped)
#'   - `theta_r`, `theta_s`, `alpha`, `n`, `Ks`, `tau`: fitted values
#'   - `std_error_*` for each free parameter (approximate, from numerical Hessian)
#'   - `convergence`: `TRUE` if `optim()` returned code 0
#'
#' @seealso [fit_swrc()], [hydraulic_conductivity()], [swrc_van_genuchten()]
#'
#' @references
#' Van Genuchten, M. Th. (1980). A closed-form equation for predicting the
#' hydraulic conductivity of unsaturated soils. *Soil Science Society of
#' America Journal*, 44(5), 892–898.
#' <https://doi.org/10.2136/sssaj1980.03615995004400050002x>
#'
#' Mualem, Y. (1976). A new model for predicting the hydraulic conductivity of
#' unsaturated porous media. *Water Resources Research*, 12(3), 513–522.
#' <https://doi.org/10.1029/WR012i003p00513>
#'
#' @examples
#' library(tibble)
#'
#' # Synthetic joint dataset — WRC and K at the same h values
#' joint <- tibble(
#'   h     = c(0, 10, 100, 500, 1000, 5000, 15000),
#'   theta = c(0.44, 0.40, 0.32, 0.25, 0.20, 0.12, 0.08),
#'   K     = c(25, 8.5, 0.9, 0.08, 0.012, 4e-4, 3e-5)
#' )
#'
#' fit_swrc_hcc(joint, theta_col = theta, K_col = K, h_col = h)
#'
#' # Fix tau at Mualem value; estimate the remaining 5 parameters
#' fit_swrc_hcc(joint, theta_col = theta, K_col = K, h_col = h,
#'              fixed = c(tau = 0.5))
#'
#' # Mixed h values: WRC and HCC measured at different points
#' mixed <- tibble(
#'   h     = c(0, 10, 100, 500, 1000, 5000,  10, 100, 1000),
#'   theta = c(0.44, 0.40, 0.32, 0.25, 0.20, 0.12, NA, NA, NA),
#'   K     = c(NA, NA, NA, NA, NA, NA, 8.5, 0.9, 0.012)
#' )
#'
#' fit_swrc_hcc(mixed, theta_col = theta, K_col = K, h_col = h)
#'
#' @export
fit_swrc_hcc <- function(data, theta_col, K_col, h_col,
                         lower = NULL, upper = NULL, fixed = NULL,
                         wrc_weight = 1, hcc_weight = 1,
                         workers = 1L) {
  theta_quo <- rlang::enquo(theta_col)
  K_quo     <- rlang::enquo(K_col)
  h_quo     <- rlang::enquo(h_col)

  theta_nm  <- rlang::as_name(theta_quo)
  K_nm      <- rlang::as_name(K_quo)
  h_nm      <- rlang::as_name(h_quo)

  for (nm in c(theta_nm, K_nm, h_nm)) {
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

  # Parameter names
  all_nms   <- c("theta_r", "theta_s", "alpha", "n", "Ks", "tau")
  fixed_nms <- names(fixed)
  free_nms  <- setdiff(all_nms, fixed_nms)

  # Default physical bounds
  default_lower <- c(theta_r = 0,        theta_s = 1e-6, alpha = 1e-6,
                     n = 1.001,           Ks      = 1e-12, tau  = -2 + 1e-6)
  default_upper <- c(theta_r = 1 - 1e-6, theta_s = 1.0,  alpha = Inf,
                     n = Inf,             Ks      = Inf,   tau  = 20)

  lower_full <- default_lower
  upper_full <- default_upper
  if (!is.null(lower)) lower_full[names(lower)] <- lower
  if (!is.null(upper)) upper_full[names(upper)] <- upper

  fixed_list <- if (is.null(fixed)) list() else as.list(setNames(unlist(fixed), fixed_nms))

  # ── Core fitter ────────────────────────────────────────────────────────────

  fit_one_group <- function(df, ...) {
    h_abs  <- abs(df[[h_nm]])
    theta  <- df[[theta_nm]]
    K_obs  <- df[[K_nm]]

    # Split into WRC and HCC subsets
    wrc_ok <- !is.na(theta) & !is.na(h_abs)
    hcc_ok <- !is.na(K_obs) & K_obs > 0 & !is.na(h_abs)

    if (sum(wrc_ok) < 2L && sum(hcc_ok) < 2L) {
      cli::cli_warn("Fewer than 2 valid observations for both WRC and HCC; skipping.")
      return(.swrc_hcc_na_row())
    }

    theta_wrc <- theta[wrc_ok];  h_wrc <- h_abs[wrc_ok]
    K_hcc     <- K_obs[hcc_ok];  h_hcc <- h_abs[hcc_ok]
    logK_hcc  <- log10(K_hcc)

    # Starting values
    start_all <- c(
      theta_r = min(theta_wrc, na.rm = TRUE) * 0.9,
      theta_s = max(theta_wrc, na.rm = TRUE) * 1.05,
      alpha   = 0.1,
      n       = 2.0,
      Ks      = max(K_hcc) * 1.1,
      tau     = if ("tau" %in% fixed_nms) fixed_list$tau else 0.5
    )
    start_free <- start_all[free_nms]

    # Joint objective function
    obj <- function(par) {
      p <- c(as.list(setNames(par, free_nms)), fixed_list)
      m <- 1 - 1 / p$n

      # WRC component
      wrc_rss <- 0
      if (sum(wrc_ok) > 0L) {
        pred_theta <- .vg_predict(h_wrc, p$theta_r, p$theta_s, p$alpha, p$n)
        wrc_rss    <- wrc_weight * sum((theta_wrc - pred_theta)^2)
      }

      # HCC component (log10 scale)
      hcc_rss <- 0
      if (sum(hcc_ok) > 0L) {
        Se        <- 1 / (1 + (p$alpha * h_hcc)^p$n)^m
        pred_K    <- p$Ks * Se^p$tau * (1 - (1 - Se^(1/m))^m)^2
        pred_K    <- pmax(pred_K, 1e-30)   # guard against log(0)
        hcc_rss   <- hcc_weight * sum((logK_hcc - log10(pred_K))^2)
      }

      wrc_rss + hcc_rss
    }

    tryCatch({
      opt <- stats::nlminb(
        start     = start_free,
        objective = obj,
        lower     = lower_full[free_nms],
        upper     = upper_full[free_nms],
        control   = list(iter.max = 500, eval.max = 1000)
      )

      n_obs  <- sum(wrc_ok) + sum(hcc_ok)
      n_free <- length(free_nms)
      sigma2   <- opt$objective / max(n_obs - n_free, 1L)
      hessian  <- tryCatch(stats::optimHess(opt$par, obj),
                           error = function(e) matrix(NA_real_, n_free, n_free))
      hess_inv <- tryCatch(solve(hessian),
                           error = function(e) matrix(NA_real_, n_free, n_free))
      se_free  <- sqrt(abs(diag(hess_inv)) * sigma2)

      p_all <- c(as.list(setNames(opt$par, free_nms)), fixed_list)
      se_all <- stats::setNames(rep(NA_real_, 6L), all_nms)
      se_all[free_nms] <- se_free

      tibble::tibble(
        theta_r = p_all$theta_r, theta_s = p_all$theta_s,
        alpha   = p_all$alpha,   n       = p_all$n,
        Ks      = p_all$Ks,      tau     = p_all$tau,
        std_error_theta_r = se_all[["theta_r"]],
        std_error_theta_s = se_all[["theta_s"]],
        std_error_alpha   = se_all[["alpha"]],
        std_error_n       = se_all[["n"]],
        std_error_Ks      = se_all[["Ks"]],
        std_error_tau     = se_all[["tau"]],
        convergence = (opt$convergence == 0L)   # nlminb: 0 = success
      )
    }, error = function(e) {
      cli::cli_warn(c("Joint WRC+HCC fit did not converge.",
                      "i" = "Original error: {conditionMessage(e)}"))
      .swrc_hcc_na_row()
    })
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

  class(result) <- c("fit_swrc_hcc", class(result))
  result
}

# NA row for failed joint fits
.swrc_hcc_na_row <- function() {
  tibble::tibble(
    theta_r = NA_real_, theta_s = NA_real_, alpha = NA_real_, n   = NA_real_,
    Ks      = NA_real_, tau     = NA_real_,
    std_error_theta_r = NA_real_, std_error_theta_s = NA_real_,
    std_error_alpha   = NA_real_, std_error_n       = NA_real_,
    std_error_Ks      = NA_real_, std_error_tau     = NA_real_,
    convergence = FALSE
  )
}
