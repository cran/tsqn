test_that("PerQn returns expected length and handles window argument", {
  set.seed(20)
  x <- as.numeric(arima.sim(list(ar = 0.4), n = 100))

  out_default <- PerQn(x)
  out_no_window <- PerQn(x, window = NULL)

  expect_true(is.numeric(out_default))
  expect_true(is.numeric(out_no_window))
  expect_equal(length(out_default), length(x) - 2L)
  expect_equal(length(out_no_window), length(x) - 2L)
  expect_true(all(is.finite(out_default)))
  expect_true(all(is.finite(out_no_window)))
})

test_that("PerQn and PerioMrob handle short inputs", {
  expect_equal(PerQn(c(1, 2)), numeric(0))
  expect_equal(PerQn(1), numeric(0))

  expect_equal(PerioMrob(1), numeric(0))
  short_out <- PerioMrob(c(1, 2))
  expect_equal(length(short_out), 1L)
  expect_true(is.finite(short_out))
  expect_true(short_out >= 0)
})

test_that("PerioMrob returns non-negative spectral estimates", {
  set.seed(21)
  x <- as.numeric(arima.sim(list(ar = 0.2), n = 80))

  out <- PerioMrob(x)
  expect_true(is.numeric(out))
  expect_equal(length(out), length(x) - 1L)
  expect_true(all(is.finite(out)))
  expect_true(all(out >= 0))
})

test_that("GPH_estimate supports all method branches and returns expected structure", {
  set.seed(22)
  x <- as.numeric(arima.sim(list(ar = 0.3), n = 120))

  out_gph <- GPH_estimate(x, method = "GPH")
  out_m <- GPH_estimate(x, method = "GPH-M")
  out_qn <- GPH_estimate(x, method = "GPH-Qn")
  out_fallback <- GPH_estimate(x, method = "something-else")

  for (out in list(out_gph, out_m, out_qn, out_fallback)) {
    expect_true(is.list(out))
    expect_equal(names(out)[1:3], c("method", "d", "sd.reg"))
    expect_equal(length(out), 4L)
    expect_true(is.numeric(out[[4]]))
    expect_true(is.finite(out$d))
    expect_true(is.finite(out$sd.reg))
  }

  expect_equal(out_gph$method, "GPH")
  expect_equal(out_m$method, "GPH-M")
  expect_equal(out_qn$method, "Qn")
  expect_equal(out_fallback$method, "Qn")
})

test_that("pm10 dataset can be loaded and has expected station names", {
  data("pm10", package = "tsqn")

  expect_true(is.data.frame(pm10))
  expect_equal(dim(pm10), c(1826L, 8L))
  expect_equal(
    names(pm10),
    c("Laranjeiras", "Carapina", "Camburi", "Sua", "VixCentro", "Ibes", "VVCentro", "Cariacica")
  )
})
