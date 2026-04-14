#' Compute soil water diffusivity (Van Genuchten / Mualem)
#'
#' Calculates soil water diffusivity D(h) at given matric potential values:
#' \deqn{D(h) = \frac{K(h)}{C(h)}}
#'
#' where K(h) is the Mualem-Van Genuchten hydraulic conductivity and C(h) is
#' the Van Genuchten soil water capacity (−dθ/dh). D(h) has units of length²
#' per time (e.g. cm²/day when `ks` is in cm/day and `h` in cm).
#'
#' D(h) appears in the diffusive form of the Richards equation:
#' \deqn{\frac{\partial \theta}{\partial t} =
#'   \nabla \cdot \left(D(\theta) \nabla \theta\right) - \frac{\partial K}{\partial z}}
#'
#' and is used in analytical solutions for infiltration and redistribution.
#'
#' @note D(h) is undefined at saturation (h = 0) because C(0) = 0.
#'   The function returns `Inf` at h = 0 unless K(0) is also zero, in which
#'   case R returns `NaN` via `0/0`. Handle h = 0 rows in advance if needed.
#'
#' @param data A data frame or tibble. Passed as the first argument so the
#'   function is compatible with `|>` and `%>%` pipes.
#' @param ks Saturated hydraulic conductivity (e.g. cm/day). Either a bare
#'   column name from `data` or a scalar numeric value. Must be strictly
#'   positive.
#' @param theta_r Residual water content (m³/m³). Either a bare column name
#'   or a scalar numeric value.
#' @param theta_s Saturated water content (m³/m³). Either a bare column name
#'   or a scalar numeric value.
#' @param alpha Van Genuchten shape parameter (1/cm or 1/m, matching units of
#'   `h`). Either a bare column name or a scalar numeric value. Must be
#'   strictly positive.
#' @param n Van Genuchten shape parameter (dimensionless). Must be > 1 for the
#'   Mualem model. Either a bare column name or a scalar numeric value.
#' @param m Van Genuchten shape parameter (dimensionless). Typically `1 - 1/n`.
#'   Must satisfy `0 < m < 1`. Either a bare column name or a scalar numeric
#'   value. Used only in the K(h) formula; C(h) derives m from n internally.
#' @param h Matric potential / pressure head (cm or m, must match units of
#'   `alpha`). Either a bare column name or a scalar numeric value. Sign
#'   convention: absolute value is used internally.
#' @param tau Tortuosity/connectivity parameter (dimensionless). Defaults to
#'   `0.5` (Mualem 1976). Must satisfy `tau > -2`. Either a bare column name
#'   or a scalar numeric value.
#'
#' @return The input `data` as a tibble with an additional column `.D`
#'   containing soil water diffusivity in units of length²/time (e.g.
#'   cm²/day when `ks` is in cm/day and `h` in cm).
#'
#' @seealso [hydraulic_conductivity()] for K(h),
#'   [soil_water_capacity()] for C(h), [swrc_van_genuchten()] for θ(h).
#'
#' @references
#' Mualem, Y. (1976). A new model for predicting the hydraulic conductivity of
#' unsaturated porous media. *Water Resources Research*, 12(3), 513–522.
#' <https://doi.org/10.1029/WR012i003p00513>
#'
#' Van Genuchten, M. Th. (1980). A closed-form equation for predicting the
#' hydraulic conductivity of unsaturated soils. *Soil Science Society of
#' America Journal*, 44(5), 892–898.
#' <https://doi.org/10.2136/sssaj1980.03615995004400050002x>
#'
#' @examples
#' library(tibble)
#'
#' # Exclude h = 0 (D is undefined at saturation)
#' soil <- tibble(h = c(1, 10, 100, 1000, 15000))
#'
#' soil_water_diffusivity(
#'   soil,
#'   ks      = 10,
#'   theta_r = 0.05,
#'   theta_s = 0.45,
#'   alpha   = 0.02,
#'   n       = 1.5,
#'   m       = 1 - 1/1.5,
#'   h       = h
#' )
#'
#' @export
soil_water_diffusivity <- function(data, ks, theta_r, theta_s, alpha, n, m, h,
                                   tau = 0.5) {
  ks_quo      <- rlang::enquo(ks)
  theta_r_quo <- rlang::enquo(theta_r)
  theta_s_quo <- rlang::enquo(theta_s)
  alpha_quo   <- rlang::enquo(alpha)
  n_quo       <- rlang::enquo(n)
  m_quo       <- rlang::enquo(m)
  h_quo       <- rlang::enquo(h)
  tau_quo     <- rlang::enquo(tau)

  .ks      <- resolve_arg(ks_quo,      data, "ks")
  .theta_r <- resolve_arg(theta_r_quo, data, "theta_r")
  .theta_s <- resolve_arg(theta_s_quo, data, "theta_s")
  .alpha   <- resolve_arg(alpha_quo,   data, "alpha")
  .n       <- resolve_arg(n_quo,       data, "n")
  .m       <- resolve_arg(m_quo,       data, "m")
  .h       <- resolve_arg(h_quo,       data, "h")
  .tau     <- resolve_arg(tau_quo,     data, "tau")

  check_positive(.ks,    "ks")
  check_theta_bounds(.theta_r, .theta_s)
  check_positive(.alpha, "alpha")
  check_n_gt_one(.n)

  if (any(!is.finite(.m)) || any(.m <= 0) || any(.m >= 1)) {
    cli::cli_abort(c(
      "{.arg m} must satisfy 0 < m < 1.",
      "x" = "Found values outside this range."
    ))
  }
  if (any(!is.finite(.tau)) || any(.tau <= -2)) {
    cli::cli_abort(c(
      "{.arg tau} must be greater than -2.",
      "x" = "Found values <= -2 or non-finite."
    ))
  }

  .h_abs <- abs(.h)
  .m_vg  <- 1 - 1 / .n    # Mualem m for C(h); distinct from explicit .m for K(h)

  # K(h) — Mualem-Van Genuchten with explicit m and tau
  Se  <- 1 / (1 + (.alpha * .h_abs)^.n)^.m
  K_h <- .ks * Se^.tau * (1 - (1 - Se^(1 / .m))^.m)^2

  # C(h) — analytical VG derivative using Mualem m = 1 - 1/n
  x   <- (.alpha * .h_abs)^.n
  C_h <- (.theta_s - .theta_r) * .m_vg * .n * .alpha^.n *
         .h_abs^(.n - 1) / (1 + x)^(.m_vg + 1)

  out    <- tibble::as_tibble(data)
  out$.D <- K_h / C_h
  out
}
