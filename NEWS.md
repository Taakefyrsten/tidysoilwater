# tidysoilwater 1.2.0

## New functions

* `fit_swrc_hcc()` — joint simultaneous fitting of the Van Genuchten SWRC and
  Mualem-Van Genuchten K(h) function. Minimises a combined weighted residual
  (θ on linear scale, K on log₁₀ scale). Estimates θ_r, θ_s, α, n, Ks, and τ.
  Supports `lower`/`upper` bounds, `fixed` parameters, `wrc_weight`/`hcc_weight`
  scalars, grouping, and parallel fitting.
* `confint.fit_swrc()` — profile-likelihood confidence intervals for parameters
  fitted by `fit_swrc()`. Uses the F-statistic profile approach (more reliable
  than normal-approximation intervals for nonlinear models). Returns a tidy
  tibble of (param, estimate, lower, upper, level) per group.

## Changes to existing functions

* `fit_swrc()` — three new optional arguments:
  - `lower` / `upper` — named vectors of box constraints; uses
    `nls(algorithm = "port")` internally (SEs from Jacobian, unchanged).
  - `fixed` — named vector of parameters to hold at specified values;
    switches to `stats::optim(method = "L-BFGS-B")` with approximate SEs
    from the numerical Hessian.
  - `weights` — bare column name or scalar for per-observation weights.
  - Return type gains class `fit_swrc` (a tibble subclass) with metadata
    attributes enabling `confint()`. Fully backward-compatible.

---

# tidysoilwater 1.1.0

## New functions

* `soil_water_capacity()` — analytical soil water capacity C(h) = −dθ/dh from
  the Van Genuchten (1980) model. Covers `soilwater::cap()`.
* `saturation_index()` — effective saturation S_e(h) = [1 + (α|h|)^n]^(−m).
  Dimensionless (0–1); equivalent to `soilwater::swc(..., saturation_index = TRUE)`.
* `soil_water_diffusivity()` — soil water diffusivity D(h) = K(h)/C(h).
  Covers `soilwater::diffusivity()`.

## Changes

* `hydraulic_conductivity()` — new `tau` argument (tortuosity/connectivity
  parameter, default `0.5`). Previously hardcoded at 0.5 (Mualem 1976);
  now exposed as a free parameter matching `soilhypfit::hc_model()`.
  Default `tau = 0.5` is fully backward-compatible.

---

# tidysoilwater 1.0.0

Initial release.

## New functions

* `swrc_van_genuchten()` — Van Genuchten (1980) soil water retention curve.
  Returns volumetric water content θ(h) for any combination of scalar
  parameters and tidy data columns. Handles 1 M rows in < 20 ms.
* `hydraulic_conductivity()` — Mualem-Van Genuchten unsaturated hydraulic
  conductivity K(h). Handles 1 M rows in < 40 ms.
* `fit_swrc()` — fits Van Genuchten parameters to observed (h, θ) data via
  nonlinear least squares. Supports grouped tibbles (one fit per group) and
  optional parallel fitting via `parallel::mclapply()` on Unix systems.

## Notes

* All functions accept bare column names (tidy evaluation) or scalar numeric
  values interchangeably.
* `fit_swrc()` with `workers > 1` delivers near-linear parallel speedup on
  Unix; silently falls back to sequential on Windows.
* `broom::tidy()` output is used internally so fitted parameters, standard
  errors, and convergence status are all returned in a single tidy tibble.
