library(tibble)

make_soil <- function(h = c(1, 10, 100, 1000, 15000)) {
  tibble(h = h)
}

test_that("D(h) = K(h) / C(h) identity holds", {
  df <- make_soil()
  params <- list(ks = 10, theta_r = 0.05, theta_s = 0.45,
                 alpha = 0.02, n = 1.5, m = 1 - 1/1.5, h_vals = df$h)

  K <- hydraulic_conductivity(df, ks = 10, alpha = 0.02, n = 1.5,
                               m = 1 - 1/1.5, h = h)$.K
  C <- soil_water_capacity(df, theta_r = 0.05, theta_s = 0.45,
                            alpha = 0.02, n = 1.5, h = h)$.C
  D_manual <- K / C

  D_fn <- soil_water_diffusivity(df, ks = 10, theta_r = 0.05, theta_s = 0.45,
                                  alpha = 0.02, n = 1.5, m = 1 - 1/1.5, h = h)$.D

  expect_equal(D_fn, D_manual, tolerance = 1e-10)
})

test_that("D(h) is positive for h > 0", {
  df <- make_soil()
  result <- soil_water_diffusivity(df, ks = 10, theta_r = 0.05, theta_s = 0.45,
                                    alpha = 0.02, n = 1.5, m = 1 - 1/1.5, h = h)
  expect_true(all(result$.D > 0))
})

test_that("tau parameter changes D(h) correctly", {
  df <- make_soil()
  D_default <- soil_water_diffusivity(df, ks = 10, theta_r = 0.05, theta_s = 0.45,
                                       alpha = 0.02, n = 1.5, m = 1 - 1/1.5,
                                       h = h, tau = 0.5)$.D
  D_tau1    <- soil_water_diffusivity(df, ks = 10, theta_r = 0.05, theta_s = 0.45,
                                       alpha = 0.02, n = 1.5, m = 1 - 1/1.5,
                                       h = h, tau = 1.0)$.D
  expect_false(isTRUE(all.equal(D_default, D_tau1)))
})

test_that("output is a tibble with .D appended", {
  df <- tibble(h = c(10, 100), id = 1:2)
  result <- soil_water_diffusivity(df, ks = 10, theta_r = 0.05, theta_s = 0.45,
                                    alpha = 0.02, n = 1.5, m = 1 - 1/1.5, h = h)
  expect_s3_class(result, "tbl_df")
  expect_true(".D"  %in% names(result))
  expect_true("id"  %in% names(result))
})

test_that("column name and scalar inputs give identical results", {
  df <- tibble(h = c(10, 100, 1000), ks = 10, theta_r = 0.05, theta_s = 0.45,
               alpha = 0.02, n = 1.5, m = 1 - 1/1.5)
  from_cols    <- soil_water_diffusivity(df, ks = ks, theta_r = theta_r,
                                          theta_s = theta_s, alpha = alpha,
                                          n = n, m = m, h = h)
  from_scalars <- soil_water_diffusivity(df, ks = 10, theta_r = 0.05,
                                          theta_s = 0.45, alpha = 0.02,
                                          n = 1.5, m = 1 - 1/1.5, h = h)
  expect_equal(from_cols$.D, from_scalars$.D)
})

test_that("invalid tau (<= -2) produces error", {
  df <- make_soil()
  expect_error(
    soil_water_diffusivity(df, ks = 10, theta_r = 0.05, theta_s = 0.45,
                            alpha = 0.02, n = 1.5, m = 1 - 1/1.5,
                            h = h, tau = -3),
    regexp = "tau"
  )
})
