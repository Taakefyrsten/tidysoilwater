test_that("at h = 0, K equals Ks", {
  df <- tibble::tibble(h = 0)
  result <- hydraulic_conductivity(df, ks = 10, alpha = 0.02,
                                   n = 1.5, m = 1 - 1/1.5, h = h)
  expect_equal(result$.K, 10)
})

test_that("K decreases monotonically as h increases", {
  df <- tibble::tibble(h = c(0, 1, 10, 100, 1000, 10000))
  result <- hydraulic_conductivity(df, ks = 10, alpha = 0.02,
                                   n = 1.5, m = 1 - 1/1.5, h = h)
  expect_true(all(diff(result$.K) <= 0))
})

test_that("K is non-negative for all finite h", {
  df <- tibble::tibble(h = c(0, 0.1, 1, 10, 100, 1000, 15000))
  result <- hydraulic_conductivity(df, ks = 10, alpha = 0.02,
                                   n = 1.5, m = 1 - 1/1.5, h = h)
  expect_true(all(result$.K >= 0))
})

test_that("column name inputs and scalar inputs return identical results", {
  df <- tibble::tibble(h = c(0, 100, 1000), ks = 10, alpha = 0.02,
                       n = 1.5, m = 1 - 1/1.5)

  from_cols <- hydraulic_conductivity(df, ks = ks, alpha = alpha, n = n,
                                      m = m, h = h)
  from_scalars <- hydraulic_conductivity(df, ks = 10, alpha = 0.02,
                                         n = 1.5, m = 1 - 1/1.5, h = h)
  expect_equal(from_cols$.K, from_scalars$.K)
})

test_that("output is a tibble with .K column appended", {
  df <- tibble::tibble(h = c(0, 100), sample_id = c(1L, 2L))
  result <- hydraulic_conductivity(df, ks = 10, alpha = 0.02,
                                   n = 1.5, m = 1 - 1/1.5, h = h)
  expect_s3_class(result, "tbl_df")
  expect_true(".K" %in% names(result))
  expect_true("sample_id" %in% names(result))
})

test_that("ks <= 0 produces an error", {
  df <- tibble::tibble(h = 100)
  expect_error(
    hydraulic_conductivity(df, ks = -5, alpha = 0.02, n = 1.5,
                           m = 1 - 1/1.5, h = h),
    regexp = "strictly positive"
  )
})

test_that("m outside (0, 1) produces an error", {
  df <- tibble::tibble(h = 100)
  expect_error(
    hydraulic_conductivity(df, ks = 10, alpha = 0.02, n = 1.5,
                           m = 1.5, h = h),
    regexp = "0 < m < 1"
  )
})

test_that("default tau = 0.5 matches previous hardcoded behaviour", {
  df <- tibble::tibble(h = c(0, 10, 100, 1000))
  explicit <- hydraulic_conductivity(df, ks = 10, alpha = 0.02, n = 1.5,
                                     m = 1 - 1/1.5, h = h, tau = 0.5)
  default  <- hydraulic_conductivity(df, ks = 10, alpha = 0.02, n = 1.5,
                                     m = 1 - 1/1.5, h = h)
  expect_equal(explicit$.K, default$.K)
})

test_that("tau != 0.5 changes K(h) values", {
  df <- tibble::tibble(h = c(10, 100, 1000))
  K_default <- hydraulic_conductivity(df, ks = 10, alpha = 0.02, n = 1.5,
                                      m = 1 - 1/1.5, h = h)$.K
  K_tau1    <- hydraulic_conductivity(df, ks = 10, alpha = 0.02, n = 1.5,
                                      m = 1 - 1/1.5, h = h, tau = 1.0)$.K
  expect_false(isTRUE(all.equal(K_default, K_tau1)))
})

test_that("tau <= -2 produces an error", {
  df <- tibble::tibble(h = 100)
  expect_error(
    hydraulic_conductivity(df, ks = 10, alpha = 0.02, n = 1.5,
                           m = 1 - 1/1.5, h = h, tau = -2),
    regexp = "tau"
  )
})
