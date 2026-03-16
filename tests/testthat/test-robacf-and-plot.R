test_that("robacf univariate output structure is correct", {
  set.seed(10)
  x <- rnorm(120)

  out <- robacf(x, lag.max = 12, type = "correlation", plot = FALSE)

  expect_s3_class(out, "robacf")
  expect_equal(out$type, "correlation")
  expect_equal(out$n.used, 120L)
  expect_equal(dim(out$acf), c(12L, 1L, 1L))
  expect_equal(dim(out$lag), c(12L, 1L, 1L))
})

test_that("robacf multivariate output structure and argument checks", {
  set.seed(11)
  mat <- cbind(rnorm(90), rnorm(90), rnorm(90))

  out <- robacf(mat, lag.max = 8, type = "covariance", plot = FALSE)
  expect_s3_class(out, "robacf")
  expect_equal(out$type, "covariance")
  expect_equal(dim(out$acf), c(8L, 3L, 3L))
  expect_equal(dim(out$lag), c(8L, 3L, 3L))

  # partial matching is supported via match.arg
  out_partial <- robacf(mat, lag.max = 5, type = "corr", plot = FALSE)
  expect_equal(out_partial$type, "correlation")

  expect_error(robacf(mat, lag.max = -1, plot = FALSE), "must be at least 0")
  expect_error(robacf(matrix(letters[1:6], ncol = 2), plot = FALSE), "'x' must be numeric")
})

test_that("plot.robacf works for univariate and multivariate objects", {
  set.seed(12)
  x <- rnorm(80)
  mat <- cbind(rnorm(80), rnorm(80))

  uni <- robacf(x, lag.max = 10, plot = FALSE)
  multi <- robacf(mat, lag.max = 6, plot = FALSE)

  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  }, add = TRUE)

  expect_invisible(tsqn:::plot.robacf(uni, main = "univariate"))
  expect_invisible(tsqn:::plot.robacf(multi, main = "multivariate"))
})
