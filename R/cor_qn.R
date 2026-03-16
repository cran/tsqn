#' Robust correlation between the variables \code{x} and \code{y}
#'
#' Computes the robust correlation of \code{x} and \code{y} proposed by Ma and Genton (2001) using the robust scale Qn of Rousseeuw and Croux (1993).
#' @param x a numeric vector
#' @param y a numeric vector
#'
#' @return a numerical value with the robust correlation between \code{x} and \code{y}
#' @references Ma, Y. and Genton, M. G. (2001) Highly robust estimation of dispersion matrices. \emph{Journal of Multivariate Analysis}, \bold{78}, 11--36.
#' @references Rousseeuw, P. J. and Croux, C. (1993) Alternatives to the median absolute deviation. \emph{Journal of the American Statistical Association}, \bold{88}, 1273--1283.
#' @import robustbase
#' @export
#' @examples
#' corQn(rnorm(100),rnorm(100))
corQn <- function(x, y){
  if(length(x) != length(y)) {
    stop("x and y are unequal sizes")
  }
  if(identical(x, y)) {
    return(1)
  }

  alpha.qn <- robustbase::Qn(x)
  beta.qn <- robustbase::Qn(y)
  qn1 <- robustbase::Qn(x / alpha.qn + y / beta.qn)
  qn2 <- robustbase::Qn(x / alpha.qn - y / beta.qn)
  qn1.sq <- qn1 * qn1
  qn2.sq <- qn2 * qn2

  (qn1.sq - qn2.sq) / (qn1.sq + qn2.sq)
}
