#' Compute the Van Genuchten soil water retention curve
#'
#' Calculates volumetric soil water content at given matric potential values
#' using the Van Genuchten (1980) model with the Mualem (1976) constraint
#' (`m = 1 - 1/n`).
#'
#' The model is:
#' \deqn{\theta(h) = \theta_r + \frac{\theta_s - \theta_r}{(1 + (\alpha |h|)^n)^m}}
#'
#' where `m = 1 - 1/n`.
#'
#' @param data A data frame or tibble. Passed as the first argument so the
#'   function is compatible with `|>` and `%>%` pipes.
#' @param theta_r Residual water content (m³/m³). Either a bare column name
#'   from `data` or a scalar numeric value.
#' @param theta_s Saturated water content (m³/m³). Either a bare column name
#'   from `data` or a scalar numeric value.
#' @param alpha Van Genuchten shape parameter (1/cm or 1/m, matching units of
#'   `h`). Either a bare column name or a scalar numeric value. Must be
#'   strictly positive.
#' @param n Van Genuchten shape parameter (dimensionless). Must be > 1 for the
#'   Mualem model. Either a bare column name or a scalar numeric value.
#' @param h Matric potential / pressure head (cm or m, must match units of
#'   `alpha`). Either a bare column name or a scalar numeric value. Sign
#'   convention: positive values are treated as suction (absolute value is
#'   used internally).
#'
#' @return The input `data` as a tibble with an additional column `.theta`
#'   containing the computed volumetric water content (m³/m³).
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
#' soil <- tibble(
#'   h     = c(0, 10, 100, 1000, 15000),
#'   alpha = 0.02,
#'   n     = 1.5
#' )
#'
#' # Using bare column names for alpha and n, scalars for theta_r / theta_s
#' swrc_van_genuchten(soil, theta_r = 0.05, theta_s = 0.45, alpha = alpha,
#'                   n = n, h = h)
#'
#' # All scalar values
#' swrc_van_genuchten(soil, theta_r = 0.05, theta_s = 0.45, alpha = 0.02,
#'                   n = 1.5, h = h)
#'
#' @export
swrc_van_genuchten <- function(data, theta_r, theta_s, alpha, n, h) {
  # Capture arguments as quosures for tidy evaluation
  theta_r_quo <- rlang::enquo(theta_r)
  theta_s_quo <- rlang::enquo(theta_s)
  alpha_quo   <- rlang::enquo(alpha)
  n_quo       <- rlang::enquo(n)
  h_quo       <- rlang::enquo(h)

  # Resolve each argument to a numeric vector or scalar
  .theta_r <- resolve_arg(theta_r_quo, data, "theta_r")
  .theta_s <- resolve_arg(theta_s_quo, data, "theta_s")
  .alpha   <- resolve_arg(alpha_quo,   data, "alpha")
  .n       <- resolve_arg(n_quo,       data, "n")
  .h       <- resolve_arg(h_quo,       data, "h")

  # Validate inputs
  check_theta_bounds(.theta_r, .theta_s)
  check_positive(.alpha, "alpha")
  check_n_gt_one(.n)

  # Mualem constraint: m = 1 - 1/n (scalars recycle in arithmetic below)
  .m <- 1 - 1 / .n

  # Van Genuchten equation — direct tibble assignment; R recycles scalars
  out        <- tibble::as_tibble(data)
  out$.theta <- .theta_r + (.theta_s - .theta_r) /
                (1 + (.alpha * abs(.h))^.n)^.m
  out
}
