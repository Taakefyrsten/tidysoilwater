#' Compute unsaturated hydraulic conductivity (Mualem-Van Genuchten)
#'
#' Calculates unsaturated hydraulic conductivity at given matric potential
#' values using the Mualem (1976) - Van Genuchten (1980) model.
#'
#' The effective saturation is first computed as:
#' \deqn{S_e(h) = \frac{1}{(1 + (\alpha |h|)^n)^m}}
#'
#' Then hydraulic conductivity is given by:
#' \deqn{K(h) = K_s \cdot S_e^{\tau} \cdot \left(1 - \left(1 - S_e^{1/m}\right)^m\right)^2}
#'
#' The tortuosity/connectivity parameter `tau` (τ) defaults to `0.5`, which is
#' the value proposed by Mualem (1976). Soilhypfit and other tools treat it as
#' a free parameter; values between −2 and 2 are physically plausible.
#'
#' @param data A data frame or tibble. Passed as the first argument so the
#'   function is compatible with `|>` and `%>%` pipes.
#' @param ks Saturated hydraulic conductivity (any consistent unit, e.g.
#'   cm/day). Either a bare column name from `data` or a scalar numeric value.
#'   Must be strictly positive.
#' @param alpha Van Genuchten shape parameter (1/cm or 1/m). Either a bare
#'   column name or a scalar numeric value. Must be strictly positive.
#' @param n Van Genuchten shape parameter (dimensionless). Either a bare
#'   column name or a scalar numeric value. Must be strictly positive.
#' @param m Van Genuchten shape parameter (dimensionless). Typically set to
#'   `1 - 1/n`, but can be specified independently. Must satisfy `0 < m < 1`.
#'   Either a bare column name or a scalar numeric value.
#' @param h Matric potential / pressure head. Either a bare column name or a
#'   scalar numeric value. Sign convention: absolute value is used internally.
#' @param tau Tortuosity/connectivity parameter (dimensionless). Defaults to
#'   `0.5` (Mualem 1976). Must satisfy `tau > -2`. Either a bare column name
#'   or a scalar numeric value.
#'
#' @return The input `data` as a tibble with an additional column `.K`
#'   containing the computed hydraulic conductivity in the same units as `ks`.
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
#' soil <- tibble(h = c(0, 10, 100, 1000, 15000))
#'
#' # Default tau = 0.5 (Mualem)
#' hydraulic_conductivity(
#'   soil,
#'   ks    = 10,
#'   alpha = 0.02,
#'   n     = 1.5,
#'   m     = 1 - 1/1.5,
#'   h     = h
#' )
#'
#' # Custom tortuosity parameter
#' hydraulic_conductivity(
#'   soil,
#'   ks    = 10,
#'   alpha = 0.02,
#'   n     = 1.5,
#'   m     = 1 - 1/1.5,
#'   h     = h,
#'   tau   = 1.0
#' )
#'
#' @export
hydraulic_conductivity <- function(data, ks, alpha, n, m, h, tau = 0.5) {
  ks_quo    <- rlang::enquo(ks)
  alpha_quo <- rlang::enquo(alpha)
  n_quo     <- rlang::enquo(n)
  m_quo     <- rlang::enquo(m)
  h_quo     <- rlang::enquo(h)
  tau_quo   <- rlang::enquo(tau)

  .ks    <- resolve_arg(ks_quo,    data, "ks")
  .alpha <- resolve_arg(alpha_quo, data, "alpha")
  .n     <- resolve_arg(n_quo,     data, "n")
  .m     <- resolve_arg(m_quo,     data, "m")
  .h     <- resolve_arg(h_quo,     data, "h")
  .tau   <- resolve_arg(tau_quo,   data, "tau")

  check_positive(.ks,    "ks")
  check_positive(.alpha, "alpha")
  check_positive(.n,     "n")

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

  Se <- 1 / (1 + (.alpha * abs(.h))^.n)^.m

  out    <- tibble::as_tibble(data)
  out$.K <- .ks * Se^.tau * (1 - (1 - Se^(1 / .m))^.m)^2
  out
}
