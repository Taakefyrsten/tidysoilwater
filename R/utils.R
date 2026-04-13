# Internal utilities for tidysoilwater

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
