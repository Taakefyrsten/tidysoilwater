# Internal utilities for tidysoilwater

# Van Genuchten prediction (vectorised; h_abs must already be abs(h))
.vg_predict <- function(h_abs, theta_r, theta_s, alpha, n) {
  m <- 1 - 1 / n
  theta_r + (theta_s - theta_r) / (1 + (alpha * h_abs)^n)^m
}

# Weighted RSS for VG model (used in optim / profile likelihood paths)
.vg_rss <- function(par, h_abs, theta_vec, w_vec) {
  pred <- .vg_predict(h_abs, par[["theta_r"]], par[["theta_s"]],
                      par[["alpha"]], par[["n"]])
  sum(w_vec * (theta_vec - pred)^2, na.rm = TRUE)
}

# Default NA row returned when any fitting path fails to converge
.swrc_na_row <- function() {
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

# Default parameter bounds (physical constraints for VG model)
.vg_default_lower <- c(theta_r = 0,      theta_s = 1e-6, alpha = 1e-6, n = 1.001)
.vg_default_upper <- c(theta_r = 1 - 1e-6, theta_s = 1.0,  alpha = Inf,  n = Inf)

# Resolve a quosure to either a column vector (if it names a column in `data`)
# or a scalar numeric value (if it evaluates to a number).
resolve_arg <- function(quo, data, arg_name) {
  expr <- rlang::get_expr(quo)
  env <- rlang::get_env(quo)

  # If the expression is a symbol, try to look it up as a column name first
  if (rlang::is_symbol(expr)) {
    nm <- rlang::as_string(expr)
    if (nm %in% names(data)) {
      return(data[[nm]])
    }
  }

  # Otherwise evaluate in caller environment
  val <- rlang::eval_tidy(quo, data = data)

  if (!is.numeric(val)) {
    cli::cli_abort(
      c(
        "{.arg {arg_name}} must be a numeric column name or a numeric scalar.",
        "x" = "Got an object of class {.cls {class(val)}}."
      )
    )
  }
  val
}

# Assert that all elements of x are strictly positive.
check_positive <- function(x, arg_name) {
  if (any(!is.finite(x)) || any(x <= 0)) {
    cli::cli_abort(
      c(
        "{.arg {arg_name}} must be strictly positive.",
        "x" = "Found non-positive or non-finite values."
      )
    )
  }
  invisible(x)
}

# Assert that theta_s > theta_r (element-wise or scalar).
check_theta_bounds <- function(theta_r, theta_s) {
  if (any(theta_s <= theta_r)) {
    cli::cli_abort(
      c(
        "{.arg theta_s} must be greater than {.arg theta_r} for all observations.",
        "x" = "Found {sum(theta_s <= theta_r)} row(s) where this is violated."
      )
    )
  }
  invisible(NULL)
}

# Assert that n > 1 (Van Genuchten constraint for Mualem model).
check_n_gt_one <- function(n, arg_name = "n") {
  if (any(!is.finite(n)) || any(n <= 1)) {
    cli::cli_abort(
      c(
        "{.arg {arg_name}} must be greater than 1 for the Mualem-Van Genuchten model.",
        "x" = "Found values <= 1."
      )
    )
  }
  invisible(n)
}
