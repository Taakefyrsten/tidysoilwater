library(tibble)
library(dplyr)

make_obs <- function() {
  tibble(
    h     = c(0, 10, 100, 500, 1000, 5000, 15000),
    theta = c(0.44, 0.40, 0.32, 0.25, 0.20, 0.12, 0.08)
  )
}

# ── Default path still works ──────────────────────────────────────────────────

test_that("default fit_swrc() still returns fit_swrc class", {
  result <- fit_swrc(make_obs(), theta_col = theta, h_col = h)
  expect_s3_class(result, "fit_swrc")
  expect_s3_class(result, "tbl_df")
})

test_that("default fit_swrc() still produces finite parameter estimates", {
  result <- fit_swrc(make_obs(), theta_col = theta, h_col = h)
  expect_true(result$convergence)
  expect_true(all(is.finite(c(result$theta_r, result$theta_s, result$alpha, result$n))))
})

# ── lower/upper (port path) ───────────────────────────────────────────────────

test_that("bounds constrain parameters to specified range", {
  result <- fit_swrc(make_obs(), theta_col = theta, h_col = h,
                     lower = c(n = 1.2), upper = c(n = 2.5))
  expect_true(result$convergence)
  expect_gte(result$n, 1.2)
  expect_lte(result$n, 2.5)
})

test_that("bounds path result is consistent with default when bounds are loose", {
  tight  <- fit_swrc(make_obs(), theta_col = theta, h_col = h)
  loose  <- fit_swrc(make_obs(), theta_col = theta, h_col = h,
                     lower = c(n = 1.001), upper = c(n = 20))
  expect_equal(tight$theta_r, loose$theta_r, tolerance = 1e-4)
  expect_equal(tight$n,       loose$n,       tolerance = 1e-4)
})

# ── fixed params (optim path) ─────────────────────────────────────────────────

test_that("fixed theta_r is exactly recovered in result", {
  result <- fit_swrc(make_obs(), theta_col = theta, h_col = h,
                     fixed = c(theta_r = 0.05))
  expect_equal(result$theta_r, 0.05)
  expect_true(is.na(result$std_error_theta_r))
})

test_that("fixed path still estimates remaining parameters", {
  result <- fit_swrc(make_obs(), theta_col = theta, h_col = h,
                     fixed = c(theta_r = 0.05))
  expect_true(result$convergence)
  expect_true(all(is.finite(c(result$theta_s, result$alpha, result$n))))
})

test_that("two fixed params leaves two free params with finite SEs", {
  result <- fit_swrc(make_obs(), theta_col = theta, h_col = h,
                     fixed = c(theta_r = 0.05, theta_s = 0.45))
  expect_equal(result$theta_r, 0.05)
  expect_equal(result$theta_s, 0.45)
  expect_true(is.na(result$std_error_theta_r))
  expect_true(is.na(result$std_error_theta_s))
  expect_true(is.finite(result$std_error_alpha))
  expect_true(is.finite(result$std_error_n))
})

# ── weights ───────────────────────────────────────────────────────────────────

test_that("weights argument runs without error (port path)", {
  df <- make_obs() |> mutate(w = c(1, 1, 1, 1, 1, 2, 2))
  result <- fit_swrc(df, theta_col = theta, h_col = h,
                     lower = c(n = 1.001), weights = w)
  expect_true(result$convergence)
})

# ── metadata attributes ───────────────────────────────────────────────────────

test_that("fit_swrc result has required metadata attributes", {
  result <- fit_swrc(make_obs(), theta_col = theta, h_col = h)
  expect_false(is.null(attr(result, ".data")))
  expect_false(is.null(attr(result, ".nms")))
  expect_false(is.null(attr(result, ".opts")))
})

test_that("grouped fit returns one row per group and has correct class", {
  df <- tibble(
    grp   = rep(c("A","B"), each = 7),
    h     = rep(c(0, 10, 100, 500, 1000, 5000, 15000), 2),
    theta = c(0.44, 0.40, 0.32, 0.25, 0.20, 0.12, 0.08,
              0.38, 0.33, 0.27, 0.21, 0.16, 0.10, 0.07)
  )
  result <- df |> group_by(grp) |> fit_swrc(theta_col = theta, h_col = h)
  expect_equal(nrow(result), 2L)
  expect_s3_class(result, "fit_swrc")
  expect_true(all(result$convergence))
})
