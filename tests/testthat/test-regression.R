test_that("pm10 dataset is available with expected dimensions", {
  data("pm10", package = "tsqn")
  expect_true(is.data.frame(pm10))
  expect_equal(nrow(pm10), 1826L)
  expect_equal(ncol(pm10), 8L)
})

test_that("core robust estimators remain stable on fixed seed", {
  set.seed(42)
  x <- rnorm(60)
  y <- rnorm(60)
  mat <- cbind(rnorm(60), rnorm(60), rnorm(60))
  ser <- as.numeric(arima.sim(list(ar = 0.25), n = 80))

  expect_equal(corQn(x, y), -0.10095025809082, tolerance = 1e-8)
  expect_equal(covQn(x, y), -0.108934907985392, tolerance = 1e-8)

  expect_equal(
    corMatQn(mat),
    structure(
      c(1, -0.176434211795965, 0.0519932431877509, -0.176434211795965,
        1, -0.177278906055383, 0.0519932431877509, -0.177278906055383,
        1),
      dim = c(3L, 3L)
    ),
    tolerance = 1e-8
  )

  expect_equal(
    covMatQn(mat),
    structure(
      c(0.849117262606964, -0.152855750145329, 0.0549394707370513,
        -0.152855750145329, 0.832082059130359, -0.179624187747065,
        0.0549394707370513, -0.179624187747065, 1.20249022202153),
      dim = c(3L, 3L)
    ),
    tolerance = 1e-8
  )

  expect_equal(
    robacf(x, lag.max = 5, type = "correlation", plot = FALSE)$acf[, 1, 1],
    c(1, -0.128687831152498, 0.0213846915659024, 0.0509693952235455, 0.0708672304610279),
    tolerance = 1e-8
  )

  expect_equal(
    robacf(mat, lag.max = 4, type = "covariance", plot = FALSE)$acf[, 1, 2],
    c(-0.152855750145328, 0.1143121898572, 0.100324874581005, -0.0679390273153274),
    tolerance = 1e-8
  )

  expect_equal(
    PerQn(ser)[1:8],
    c(0.0169010819338127, 0.212153296429965, 0.346746619201623, 0.320021148826367,
      0.207554970829092, 0.145510974284028, 0.140598132101839, 0.124468542798495),
    tolerance = 1e-8
  )

  expect_equal(
    PerioMrob(ser)[1:6],
    c(0.058132096367218, 0.0719656455815886, 0.648076704731737, 0.259249586943789,
      0.196552911673424, 0.164786748972315),
    tolerance = 1e-8
  )

  expect_equal(GPH_estimate(ser, method = "GPH")$d, 0.0741271121562797, tolerance = 1e-8)
  expect_equal(GPH_estimate(ser, method = "GPH-M")$d, 0.00980072784579627, tolerance = 1e-8)
  expect_equal(GPH_estimate(ser, method = "GPH-Qn")$d, -0.155590638966208, tolerance = 1e-8)
})
