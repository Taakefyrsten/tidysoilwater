#' Compute effective saturation (Van Genuchten)
#'
#' Calculates the effective (relative) saturation S_e(h) at given matric
#' potential values using the Van Genuchten (1980) model with the Mualem
#' constraint (m = 1 − 1/n):
#'
#' \deqn{S_e(h) = \frac{1}{(1 + (\alpha |h|)^n)^m}}
#'
#' S_e ranges from 0 (completely dry) to 1 (saturated). It is related to
#' volumetric water content by:
#'
#' \deqn{\theta(h) = \theta_r + (\theta_s - \theta_r) \cdot S_e(h)}
#'
#' Use this function when you need the dimensionless saturation state rather
#' than the absolute water content (e.g. as an intermediate for K(h) or when
#' comparing soils with different θ_r/θ_s).
#'
#' @param data A data frame or tibble. Passed as the first argument so the
#'   function is compatible with `|>` and `%>%` pipes.
#' @param alpha Van Genuchten shape parameter (1/cm or 1/m, matching units of
#'   `h`). Either a bare column name from `data` or a scalar numeric value.
#'   Must be strictly positive.
#' @param n Van Genuchten shape parameter (dimensionless). Must be > 1 for the
#'   Mualem model. Either a bare column name or a scalar numeric value.
#' @param h Matric potential / pressure head (cm or m, must match units of
#'   `alpha`). Either a bare column name or a scalar numeric value. Sign
#'   convention: absolute value is used internally.
#'
#' @return The input `data` as a tibble with an additional column `.Se`
#'   containing effective saturation (dimensionless, 0–1).
#'
#' @seealso [swrc_van_genuchten()] for absolute water content θ(h),
#'   [hydraulic_conductivity()] for K(h).
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
#' soil <- tibble(h = c(0, 10, 100, 1000, 15000))
#'
#' saturation_index(soil, alpha = 0.02, n = 1.5, h = h)
#'
#' # Confirm relationship to swrc_van_genuchten()
#' soil |>
#'   saturation_index(alpha = 0.02, n = 1.5, h = h) |>
#'   dplyr::mutate(.theta_check = 0.05 + (0.45 - 0.05) * .Se)
#'
#' @export
saturation_index <- function(data, alpha, n, h) {
  alpha_quo <- rlang::enquo(alpha)
  n_quo     <- rlang::enquo(n)
  h_quo     <- rlang::enquo(h)

  .alpha <- resolve_arg(alpha_quo, data, "alpha")
  .n     <- resolve_arg(n_quo,     data, "n")
  .h     <- resolve_arg(h_quo,     data, "h")

  check_positive(.alpha, "alpha")
  check_n_gt_one(.n)

  .m <- 1 - 1 / .n

  out     <- tibble::as_tibble(data)
  out$.Se <- 1 / (1 + (.alpha * abs(.h))^.n)^.m
  out
}
