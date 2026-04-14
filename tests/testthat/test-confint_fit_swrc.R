library(tibble)
library(dplyr)

make_obs <- function() {
  tibble(
    h     = c(0, 10, 100, 500, 1000, 5000, 15000),
    theta = c(0.44, 0.40, 0.32, 0.25, 0.20, 0.12, 0.08)
  )
}

test_that("confint() returns a tibble with expected columns", {
  fit    <- fit_swrc(make_obs(), theta_col = theta, h_col = h)
  result <- confint(fit)
  expect_s3_class(result, "tbl_df")
  expect_true(all(c("param", "estimate", "lower", "upper", "level") %in% names(result)))
})

test_that("confint() returns one row per free parameter by default", {
  fit    <- fit_swrc(make_obs(), theta_col = theta, h_col = h)
  result <- confint(fit)
  # 4 free params by default
  expect_equal(nrow(result), 4L)
  expect_setequal(result$param, c("theta_r", "theta_s", "alpha", "n"))
})

test_that("lower <= estimate <= upper for all parameters", {
  # lower == estimate is valid when a parameter converges to its lower bound
  fit    <- fit_swrc(make_obs(), theta_col = theta, h_col = h)
  result <- confint(fit)
  expect_true(all(result$lower <= result$estimate, na.rm = TRUE))
  expect_true(all(result$estimate <= result$upper, na.rm = TRUE))
})

test_that("parm argument subsets parameters correctly", {
  fit    <- fit_swrc(make_obs(), theta_col = theta, h_col = h)
  result <- confint(fit, parm = c("alpha", "n"))
  expect_equal(nrow(result), 2L)
  expect_setequal(result$param, c("alpha", "n"))
})

test_that("level argument is reflected in output", {
  fit    <- fit_swrc(make_obs(), theta_col = theta, h_col = h)
  result <- confint(fit, parm = "n", level = 0.90)
  expect_equal(unique(result$level), 0.90)
})

test_that("90% CI is narrower than 95% CI", {
  fit     <- fit_swrc(make_obs(), theta_col = theta, h_col = h)
  ci_90   <- confint(fit, parm = "n", level = 0.90)
  ci_95   <- confint(fit, parm = "n", level = 0.95)
  width90 <- ci_90$upper - ci_90$lower
  width95 <- ci_95$upper - ci_95$lower
  expect_lt(width90, width95)
})

test_that("fixed parameter is excluded from confint()", {
  fit    <- fit_swrc(make_obs(), theta_col = theta, h_col = h,
                     fixed = c(theta_r = 0.05))
  result <- confint(fit)
  expect_false("theta_r" %in% result$param)
  expect_equal(nrow(result), 3L)  # only 3 free params
})

test_that("grouped confint returns CIs for each group", {
  # Use data with theta_r clearly above zero so PORT converges without boundary
  df <- tibble(
    grp   = rep(c("A","B"), each = 7),
    h     = rep(c(0, 10, 100, 500, 1000, 5000, 15000), 2),
    theta = c(0.45, 0.42, 0.35, 0.28, 0.23, 0.15, 0.11,
              0.40, 0.37, 0.30, 0.23, 0.18, 0.11, 0.08)
  )
  fit    <- suppressWarnings(
    df |> group_by(grp) |> fit_swrc(theta_col = theta, h_col = h)
  )
  result <- confint(fit, parm = c("alpha", "n"))
  # 2 groups × 2 params = 4 rows
  expect_equal(nrow(result), 4L)
  expect_true("grp" %in% names(result))
})

test_that("confint() errors informatively on non-fit_swrc object", {
  expect_error(confint(tibble(x = 1)), class = "error")
})
