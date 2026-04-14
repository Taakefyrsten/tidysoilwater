library(tibble)

test_that("C(0) = 0 at saturation", {
  df <- tibble(h = 0)
  result <- soil_water_capacity(df, theta_r = 0.05, theta_s = 0.45,
                                 alpha = 0.02, n = 1.5, h = h)
  expect_equal(result$.C, 0)
})

test_that("C(h) is non-negative for all positive h", {
  df <- tibble(h = c(0, 1, 10, 100, 1000, 15000))
  result <- soil_water_capacity(df, theta_r = 0.05, theta_s = 0.45,
                                 alpha = 0.02, n = 1.5, h = h)
  expect_true(all(result$.C >= 0))
})

test_that("C(h) has a single maximum (unimodal) over h range", {
  df <- tibble(h = seq(0.1, 5000, length.out = 500))
  result <- soil_water_capacity(df, theta_r = 0.05, theta_s = 0.45,
                                 alpha = 0.02, n = 1.5, h = h)
  # One peak: values increase then decrease
  i_max <- which.max(result$.C)
  expect_true(i_max > 1 && i_max < nrow(result))
})

test_that("C(h) matches finite-difference approximation of dtheta/dh", {
  dh  <- 1e-4
  h0  <- 100
  df  <- tibble(h = c(h0 - dh, h0 + dh))
  th  <- swrc_van_genuchten(df, theta_r = 0.05, theta_s = 0.45,
                             alpha = 0.02, n = 1.5, h = h)$.theta
  C_fd <- (th[1] - th[2]) / (2 * dh)   # -dθ/dh via central difference

  df2 <- tibble(h = h0)
  C_an <- soil_water_capacity(df2, theta_r = 0.05, theta_s = 0.45,
                               alpha = 0.02, n = 1.5, h = h)$.C

  expect_equal(C_an, C_fd, tolerance = 1e-5)
})

test_that("column name and scalar inputs give identical results", {
  df <- tibble(h = c(10, 100, 1000), alpha = 0.02, n = 1.5,
               theta_r = 0.05, theta_s = 0.45)
  from_cols    <- soil_water_capacity(df, theta_r = theta_r, theta_s = theta_s,
                                       alpha = alpha, n = n, h = h)
  from_scalars <- soil_water_capacity(df, theta_r = 0.05, theta_s = 0.45,
                                       alpha = 0.02, n = 1.5, h = h)
  expect_equal(from_cols$.C, from_scalars$.C)
})

test_that("output is a tibble with .C appended", {
  df <- tibble(h = c(10, 100), id = 1:2)
  result <- soil_water_capacity(df, theta_r = 0.05, theta_s = 0.45,
                                 alpha = 0.02, n = 1.5, h = h)
  expect_s3_class(result, "tbl_df")
  expect_true(".C" %in% names(result))
  expect_true("id" %in% names(result))
})

test_that("theta_s <= theta_r produces error", {
  df <- tibble(h = 100)
  expect_error(
    soil_water_capacity(df, theta_r = 0.5, theta_s = 0.3,
                        alpha = 0.02, n = 1.5, h = h),
    regexp = "theta_s"
  )
})

test_that("sign convention: positive and negative h give same result", {
  df_pos <- tibble(h =  100)
  df_neg <- tibble(h = -100)
  C_pos <- soil_water_capacity(df_pos, theta_r = 0.05, theta_s = 0.45,
                                alpha = 0.02, n = 1.5, h = h)$.C
  C_neg <- soil_water_capacity(df_neg, theta_r = 0.05, theta_s = 0.45,
                                alpha = 0.02, n = 1.5, h = h)$.C
  expect_equal(C_pos, C_neg)
})
