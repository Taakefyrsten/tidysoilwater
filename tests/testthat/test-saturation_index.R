library(tibble)

test_that("Se = 1 at h = 0 (saturated)", {
  df <- tibble(h = 0)
  result <- saturation_index(df, alpha = 0.02, n = 1.5, h = h)
  expect_equal(result$.Se, 1)
})

test_that("Se approaches 0 at very high h", {
  # Se = [1 + (alpha*h)^n]^(-m); for alpha=0.02, n=1.5, h=1e8: Se ~ 7e-4
  df <- tibble(h = 1e8)
  result <- saturation_index(df, alpha = 0.02, n = 1.5, h = h)
  expect_true(result$.Se < 0.01)
})

test_that("Se is strictly decreasing with increasing h", {
  df <- tibble(h = c(0, 1, 10, 100, 1000, 15000))
  result <- saturation_index(df, alpha = 0.02, n = 1.5, h = h)
  expect_true(all(diff(result$.Se) <= 0))
})

test_that("Se is in [0, 1] for all h", {
  df <- tibble(h = c(0, 0.1, 1, 100, 10000))
  result <- saturation_index(df, alpha = 0.02, n = 1.5, h = h)
  expect_true(all(result$.Se >= 0 & result$.Se <= 1))
})

test_that("Se recovers theta via theta_r + (theta_s - theta_r) * Se", {
  df <- tibble(h = c(10, 100, 1000))
  theta_r <- 0.05; theta_s <- 0.45
  Se    <- saturation_index(df, alpha = 0.02, n = 1.5, h = h)$.Se
  theta_from_Se <- theta_r + (theta_s - theta_r) * Se
  theta_direct  <- swrc_van_genuchten(df, theta_r = theta_r, theta_s = theta_s,
                                       alpha = 0.02, n = 1.5, h = h)$.theta
  expect_equal(theta_from_Se, theta_direct, tolerance = 1e-12)
})

test_that("sign convention: positive and negative h give same Se", {
  df_pos <- tibble(h =  100)
  df_neg <- tibble(h = -100)
  Se_pos <- saturation_index(df_pos, alpha = 0.02, n = 1.5, h = h)$.Se
  Se_neg <- saturation_index(df_neg, alpha = 0.02, n = 1.5, h = h)$.Se
  expect_equal(Se_pos, Se_neg)
})

test_that("column name and scalar inputs give identical results", {
  df <- tibble(h = c(10, 100, 1000), alpha = 0.02, n = 1.5)
  from_cols    <- saturation_index(df, alpha = alpha, n = n, h = h)
  from_scalars <- saturation_index(df, alpha = 0.02, n = 1.5, h = h)
  expect_equal(from_cols$.Se, from_scalars$.Se)
})

test_that("output is a tibble with .Se appended", {
  df <- tibble(h = c(10, 100), id = 1:2)
  result <- saturation_index(df, alpha = 0.02, n = 1.5, h = h)
  expect_s3_class(result, "tbl_df")
  expect_true(".Se" %in% names(result))
  expect_true("id"  %in% names(result))
})
