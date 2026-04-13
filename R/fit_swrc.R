#' Fit Van Genuchten parameters to observed retention data
#'
#' Fits the Van Genuchten (1980) soil water retention model to observed (h, θ)
#' data using nonlinear least squares (`nls()`). Returns a tidy tibble of
#' estimated parameters. If the input data is grouped (via [dplyr::group_by()]),
#' the fit is performed independently for each group.
#'
#' Starting values for `nls()` are automatically initialised from the data:
#' - `theta_r`: 90% of the minimum observed θ
#' - `theta_s`: 105% of the maximum observed θ
#' - `alpha`: 0.1
#' - `n`: 2.0
#'
#' @param data A data frame or tibble, optionally grouped with
#'   [dplyr::group_by()].
#' @param theta_col Bare column name of observed volumetric water content
#'   (m³/m³).
#' @param h_col Bare column name of observed matric potential / pressure head.
#'   Absolute values are used internally (sign convention: suction positive).
#' @param workers Number of parallel workers for fitting multiple groups
#'   simultaneously. Defaults to `1` (sequential). Values greater than `1`
#'   use [parallel::mclapply()] on Unix-like systems. On Windows, any value
#'   greater than `1` is silently reduced to `1` because forking is
#'   unavailable; use [parallel::detectCores()] to choose an appropriate
#'   number on Unix. Parallelism is only beneficial when the data contains
#'   many groups (typically ≥ 50).
#'
#' @return A tibble with one row per group (or one row for ungrouped data)
#'   containing:
#'   - Group keys (if input was grouped)
#'   - `theta_r`, `theta_s`, `alpha`, `n`: fitted parameter values
#'   - `std_error_theta_r`, `std_error_theta_s`, `std_error_alpha`,
#'     `std_error_n`: standard errors of the fitted parameters
#'   - `convergence`: logical, `TRUE` if `nls()` converged successfully
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
#' fit_swrc(observed, theta_col = theta, h_col = h)
#'
#' # Grouped fit (one set of parameters per soil horizon)
#' library(dplyr)
#'
#' multi_horizon <- tibble(
#'   horizon = rep(c("A", "B"), each = 7),
#'   h       = rep(c(0, 10, 100, 500, 1000, 5000, 15000), 2),
#'   theta   = c(0.44, 0.40, 0.32, 0.25, 0.20, 0.12, 0.08,
#'               0.38, 0.33, 0.27, 0.21, 0.16, 0.10, 0.07)
#' )
#'
#' multi_horizon |>
#'   group_by(horizon) |>
#'   fit_swrc(theta_col = theta, h_col = h)
#'
#' # Parallel fit across many groups using all physical cores.
#' # Each group is fitted independently via parallel::mclapply(), giving
#' # near-linear speedup on Unix-like systems.  On Windows, workers is
#' # silently clamped to 1.
#' \dontrun{
#' n_cores <- parallel::detectCores(logical = FALSE)
#'
#' multi_horizon |>
#'   group_by(horizon) |>
#'   fit_swrc(theta_col = theta, h_col = h, workers = n_cores)
#' }
#'
#' @export
fit_swrc <- function(data, theta_col, h_col, workers = 1L) {
  theta_quo <- rlang::enquo(theta_col)
  h_quo     <- rlang::enquo(h_col)

  theta_nm <- rlang::as_name(theta_quo)
  h_nm     <- rlang::as_name(h_quo)

  if (!theta_nm %in% names(data)) {
    cli::cli_abort(c(
      "Column {.val {theta_nm}} not found in {.arg data}.",
      "i" = "Available columns: {.val {names(data)}}."
    ))
  }
  if (!h_nm %in% names(data)) {
    cli::cli_abort(c(
      "Column {.val {h_nm}} not found in {.arg data}.",
      "i" = "Available columns: {.val {names(data)}}."
    ))
  }

  workers <- as.integer(workers)
  if (workers < 1L) workers <- 1L

  # Windows cannot fork; mclapply silently falls back to sequential there
  # but we make it explicit to avoid a confusing mclapply warning.
  if (workers > 1L && .Platform$OS.type == "windows") {
    cli::cli_warn(c(
      "Parallel fitting is not supported on Windows.",
      "i" = "Falling back to {.code workers = 1} (sequential)."
    ))
    workers <- 1L
  }

  # Core fit for a single group data frame ---------------------------------
  fit_one_group <- function(df, ...) {
    theta_vec <- df[[theta_nm]]
    h_vec     <- abs(df[[h_nm]])

    fit_df        <- df
    fit_df$h_abs  <- h_vec

    vg_formula <- stats::as.formula(
      paste0(theta_nm, " ~ theta_r + (theta_s - theta_r) /",
             " (1 + (alpha * h_abs)^n)^(1 - 1/n)")
    )

    tryCatch(
      {
        fit      <- stats::nls(
          formula = vg_formula,
          data    = fit_df,
          start   = list(
            theta_r = min(theta_vec, na.rm = TRUE) * 0.9,
            theta_s = max(theta_vec, na.rm = TRUE) * 1.05,
            alpha   = 0.1,
            n       = 2.0
          ),
          control = stats::nls.control(maxiter = 200, tol = 1e-6)
        )
        tidy_fit <- broom::tidy(fit)
        params   <- stats::setNames(tidy_fit$estimate, tidy_fit$term)
        se       <- stats::setNames(tidy_fit$std.error, tidy_fit$term)

        tibble::tibble(
          theta_r           = params[["theta_r"]],
          theta_s           = params[["theta_s"]],
          alpha             = params[["alpha"]],
          n                 = params[["n"]],
          std_error_theta_r = se[["theta_r"]],
          std_error_theta_s = se[["theta_s"]],
          std_error_alpha   = se[["alpha"]],
          std_error_n       = se[["n"]],
          convergence       = TRUE
        )
      },
      error = function(e) {
        cli::cli_warn(c(
          "Van Genuchten fitting did not converge.",
          "i" = "Original error: {conditionMessage(e)}"
        ))
        tibble::tibble(
          theta_r           = NA_real_,
          theta_s           = NA_real_,
          alpha             = NA_real_,
          n                 = NA_real_,
          std_error_theta_r = NA_real_,
          std_error_theta_s = NA_real_,
          std_error_alpha   = NA_real_,
          std_error_n       = NA_real_,
          convergence       = FALSE
        )
      }
    )
  }

  # Dispatch ---------------------------------------------------------------
  is_grouped <- dplyr::is_grouped_df(data)

  if (!is_grouped) {
    # Single group — no parallelism to exploit
    return(fit_one_group(data))
  }

  if (workers == 1L) {
    # Sequential tidy path
    return(dplyr::group_modify(data, fit_one_group) |> dplyr::ungroup())
  }

  # Parallel path: split → mclapply → reassemble with group keys
  group_dfs  <- dplyr::group_split(data)
  group_keys <- dplyr::group_keys(data)

  fits <- parallel::mclapply(group_dfs, fit_one_group, mc.cores = workers)

  dplyr::bind_cols(group_keys, dplyr::bind_rows(fits))
}
