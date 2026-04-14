library(tibble)
library(dplyr)

make_joint <- function() {
  tibble(
    h     = c(0, 10, 100, 500, 1000, 5000, 15000),
    theta = c(0.44, 0.40, 0.32, 0.25, 0.20, 0.12, 0.08),
    K     = c(25, 8.5, 0.9, 0.08, 0.012, 4e-4, 3e-5)
  )
}

test_that("fit_swrc_hcc returns fit_swrc_hcc class", {
  result <- fit_swrc_hcc(make_joint(), theta_col = theta, K_col = K, h_col = h)
  expect_s3_class(result, "fit_swrc_hcc")
  expect_s3_class(result, "tbl_df")
})

test_that("fit_swrc_hcc returns all 6 parameters", {
  result <- fit_swrc_hcc(make_joint(), theta_col = theta, K_col = K, h_col = h)
  expect_true(all(c("theta_r","theta_s","alpha","n","Ks","tau") %in% names(result)))
})

test_that("fit_swrc_hcc converges on clean synthetic data", {
  result <- fit_swrc_hcc(make_joint(), theta_col = theta, K_col = K, h_col = h)
  expect_true(result$convergence)
  expect_true(all(is.finite(c(result$theta_r, result$theta_s,
                               result$alpha, result$n, result$Ks))))
})

test_that("theta_r < theta_s in fitted result", {
  result <- fit_swrc_hcc(make_joint(), theta_col = theta, K_col = K, h_col = h)
  expect_lt(result$theta_r, result$theta_s)
})

test_that("Ks is positive in fitted result", {
  result <- fit_swrc_hcc(make_joint(), theta_col = theta, K_col = K, h_col = h)
  expect_gt(result$Ks, 0)
})

test_that("fixed tau is exactly recovered", {
  result <- fit_swrc_hcc(make_joint(), theta_col = theta, K_col = K, h_col = h,
                         fixed = c(tau = 0.5))
  expect_equal(result$tau, 0.5)
  expect_true(is.na(result$std_error_tau))
})

test_that("WRC-only rows (NA in K) are handled correctly", {
  mixed <- tibble(
    h     = c(0, 10, 100, 500, 1000, 5000, 15000, 10, 100, 1000),
    theta = c(0.44, 0.40, 0.32, 0.25, 0.20, 0.12, 0.08, NA, NA, NA),
    K     = c(NA, NA, NA, NA, NA, NA, NA, 8.5, 0.9, 0.012)
  )
  result <- fit_swrc_hcc(mixed, theta_col = theta, K_col = K, h_col = h)
  expect_true(result$convergence)
  expect_true(is.finite(result$Ks))
})

test_that("grouped fit returns one row per group", {
  df <- tibble(
    grp   = rep(c("A","B"), each = 7),
    h     = rep(c(0, 10, 100, 500, 1000, 5000, 15000), 2),
    theta = c(0.44, 0.40, 0.32, 0.25, 0.20, 0.12, 0.08,
              0.38, 0.33, 0.27, 0.21, 0.16, 0.10, 0.07),
    K     = c(25, 8.5, 0.9, 0.08, 0.012, 4e-4, 3e-5,
              15, 5.0, 0.5, 0.05, 0.008, 2e-4, 1e-5)
  )
  result <- df |> group_by(grp) |>
    fit_swrc_hcc(theta_col = theta, K_col = K, h_col = h,
                 fixed = c(tau = 0.5))
  expect_equal(nrow(result), 2L)
  expect_true(all(result$convergence))
})

test_that("wrc_weight and hcc_weight are accepted", {
  result <- fit_swrc_hcc(make_joint(), theta_col = theta, K_col = K, h_col = h,
                         wrc_weight = 2, hcc_weight = 0.5)
  expect_true(result$convergence)
})

test_that("missing K column name produces informative error", {
  expect_error(
    fit_swrc_hcc(make_joint(), theta_col = theta, K_col = bad_col, h_col = h),
    regexp = "not found"
  )
})
