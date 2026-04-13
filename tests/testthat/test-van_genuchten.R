test_that("at h = 0, theta equals theta_s", {
  df <- tibble::tibble(h = 0)
  result <- swrc_van_genuchten(df, theta_r = 0.05, theta_s = 0.45,
                               alpha = 0.02, n = 1.5, h = h)
  expect_equal(result$.theta, 0.45)
})

test_that("at very large h, theta approaches theta_r", {
  # At h = 1e9 cm the curve is numerically indistinguishable from theta_r
  # for any soil; tolerance 1e-4 is generous but correct for double precision.
  df <- tibble::tibble(h = 1e9)
  result <- swrc_van_genuchten(df, theta_r = 0.05, theta_s = 0.45,
                               alpha = 0.02, n = 1.5, h = h)
  expect_lt(abs(result$.theta - 0.05), 1e-4)
})

test_that("theta decreases monotonically as h increases", {
  df <- tibble::tibble(h = c(0, 1, 10, 100, 1000, 10000))
  result <- swrc_van_genuchten(df, theta_r = 0.05, theta_s = 0.45,
                               alpha = 0.02, n = 1.5, h = h)
  expect_true(all(diff(result$.theta) <= 0))
})

test_that("scalar column name and scalar numeric value give same result", {
  df <- tibble::tibble(h = c(0, 100, 1000), alpha = 0.02, n = 1.5)

  from_cols <- swrc_van_genuchten(df, theta_r = 0.05, theta_s = 0.45,
                                  alpha = alpha, n = n, h = h)
  from_scalars <- swrc_van_genuchten(df, theta_r = 0.05, theta_s = 0.45,
                                     alpha = 0.02, n = 1.5, h = h)
  expect_equal(from_cols$.theta, from_scalars$.theta)
})

test_that("output is a tibble with .theta column appended", {
  df <- tibble::tibble(h = c(0, 100), site = c("A", "B"))
  result <- swrc_van_genuchten(df, theta_r = 0.05, theta_s = 0.45,
                               alpha = 0.02, n = 1.5, h = h)
  expect_s3_class(result, "tbl_df")
  expect_true(".theta" %in% names(result))
  expect_true("site" %in% names(result))
})

test_that("theta_s <= theta_r produces an error", {
  df <- tibble::tibble(h = 100)
  expect_snapshot(
    swrc_van_genuchten(df, theta_r = 0.45, theta_s = 0.05,
                       alpha = 0.02, n = 1.5, h = h),
    error = TRUE
  )
})

test_that("n <= 1 produces an error", {
  df <- tibble::tibble(h = 100)
  expect_snapshot(
    swrc_van_genuchten(df, theta_r = 0.05, theta_s = 0.45,
                       alpha = 0.02, n = 0.8, h = h),
    error = TRUE
  )
})

test_that("negative h is treated identically to positive h (absolute value)", {
  df_pos <- tibble::tibble(h = 100)
  df_neg <- tibble::tibble(h = -100)
  res_pos <- swrc_van_genuchten(df_pos, theta_r = 0.05, theta_s = 0.45,
                                alpha = 0.02, n = 1.5, h = h)
  res_neg <- swrc_van_genuchten(df_neg, theta_r = 0.05, theta_s = 0.45,
                                alpha = 0.02, n = 1.5, h = h)
  expect_equal(res_pos$.theta, res_neg$.theta)
})
