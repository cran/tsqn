test_that("corQn and covQn basic properties and errors", {
  set.seed(123)
  x <- rnorm(80)
  y <- rnorm(80)

  expect_equal(corQn(x, x), 1)
  expect_equal(covQn(x, x), robustbase::Qn(x)^2, tolerance = 1e-12)

  expect_equal(corQn(x, y), corQn(y, x), tolerance = 1e-12)
  expect_equal(covQn(x, y), covQn(y, x), tolerance = 1e-12)

  expect_error(corQn(x, y[-1]), "unequal sizes")
  expect_error(covQn(x, y[-1]), "unequal sizes")
})

test_that("corMatQn and covMatQn return symmetric square matrices", {
  set.seed(321)
  mat <- cbind(rnorm(100), rnorm(100), rnorm(100), rnorm(100))

  cmat <- corMatQn(mat)
  vmat <- covMatQn(mat)

  expect_true(is.matrix(cmat))
  expect_true(is.matrix(vmat))
  expect_equal(dim(cmat), c(4L, 4L))
  expect_equal(dim(vmat), c(4L, 4L))

  expect_equal(cmat, t(cmat), tolerance = 1e-12)
  expect_equal(vmat, t(vmat), tolerance = 1e-12)
  expect_equal(diag(cmat), rep(1, 4), tolerance = 1e-12)

  qn_scales <- apply(mat, 2, robustbase::Qn)
  expect_equal(diag(vmat), qn_scales^2, tolerance = 1e-10)
})

test_that("matrix estimators handle one-column and non-numeric input", {
  set.seed(777)
  one_col <- matrix(rnorm(40), ncol = 1)

  expect_equal(corMatQn(one_col), matrix(1, nrow = 1, ncol = 1))
  expect_equal(covMatQn(one_col), matrix(robustbase::Qn(one_col[, 1])^2, nrow = 1, ncol = 1))

  bad <- matrix(letters[1:6], ncol = 2)
  expect_error(corMatQn(bad), "'x' must be numeric")
  expect_error(covMatQn(bad), "'x' must be numeric")
})
