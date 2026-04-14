#' Compute soil water capacity (Van Genuchten)
#'
#' Calculates the soil water capacity C(h) = −dθ/dh at given matric potential
#' values using the analytical derivative of the Van Genuchten (1980) model
#' with the Mualem constraint (m = 1 − 1/n).
#'
#' The soil water capacity is the slope of the soil water retention curve:
#' \deqn{C(h) = -\frac{d\theta}{dh} =
#'   \frac{(\theta_s - \theta_r) \cdot m \cdot n \cdot \alpha^n \cdot |h|^{n-1}}
#'        {(1 + (\alpha |h|)^n)^{m+1}}}
#'
#' C(h) is zero at saturation (h = 0) and at very high suction (h → ∞), and
#' reaches a maximum at the inflection point of the retention curve.
#'
#' C(h) is needed for computing soil water diffusivity D(h) = K(h)/C(h) and
#' appears in the Richards equation for unsaturated flow.
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
#'   convention: absolute value is used internally.
#'
#' @return The input `data` as a tibble with an additional column `.C`
#'   containing soil water capacity in units of 1/cm (or 1/m, matching
#'   the units of `h` and `alpha`). At h = 0, `.C` is 0.
#'
#' @seealso [swrc_van_genuchten()] for θ(h), [soil_water_diffusivity()] for
#'   D(h) = K(h)/C(h), [hydraulic_conductivity()] for K(h).
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
#' soil_water_capacity(soil, theta_r = 0.05, theta_s = 0.45,
#'                     alpha = 0.02, n = 1.5, h = h)
#'
#' # Pipeline with swrc_van_genuchten()
#' soil |>
#'   swrc_van_genuchten(theta_r = 0.05, theta_s = 0.45,
#'                      alpha = 0.02, n = 1.5, h = h) |>
#'   soil_water_capacity(theta_r = 0.05, theta_s = 0.45,
#'                       alpha = 0.02, n = 1.5, h = h)
#'
#' @export
soil_water_capacity <- function(data, theta_r, theta_s, alpha, n, h) {
  theta_r_quo <- rlang::enquo(theta_r)
  theta_s_quo <- rlang::enquo(theta_s)
  alpha_quo   <- rlang::enquo(alpha)
  n_quo       <- rlang::enquo(n)
  h_quo       <- rlang::enquo(h)

  .theta_r <- resolve_arg(theta_r_quo, data, "theta_r")
  .theta_s <- resolve_arg(theta_s_quo, data, "theta_s")
  .alpha   <- resolve_arg(alpha_quo,   data, "alpha")
  .n       <- resolve_arg(n_quo,       data, "n")
  .h       <- resolve_arg(h_quo,       data, "h")

  check_theta_bounds(.theta_r, .theta_s)
  check_positive(.alpha, "alpha")
  check_n_gt_one(.n)

  .m     <- 1 - 1 / .n
  .h_abs <- abs(.h)

  # Analytical derivative of Van Genuchten SWRC: C(h) = -dθ/dh
  # At h=0: h_abs^(n-1) = 0^(n-1) = 0 (for n > 1), so C(0) = 0 correctly
  x <- (.alpha * .h_abs)^.n

  out    <- tibble::as_tibble(data)
  out$.C <- (.theta_s - .theta_r) * .m * .n * .alpha^.n *
            .h_abs^(.n - 1) / (1 + x)^(.m + 1)
  out
}
