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
