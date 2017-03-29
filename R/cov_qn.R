#' Robust covariance between the variables \code{x} and \code{y}
#'
#' Computes the robust covariance of \code{x} and \code{y} proposed by Ma and Genton (2001) using the robust scale Qn of Rousseeuw and Croux (1993).
#' @param x a numeric vector
#' @param y a numeric vector
#'
#' @return a numerical value with the robust covariance between \code{x} and \code{y}
#' @references Ma, Y. and Genton, M. G. (2001) Highly robust estimation of dispersion matrices. \emph{Journal of Multivariate Analysis}, \bold{78}, 11--36.
#' @references Rousseeuw, P. J. and Croux, C. (1993) Alternatives to the median absolute deviation. \emph{Journal of the American Statistical Association}, \bold{88}, 1273--1283.
#' @import robustbase
#' @export
#' @examples
#' covQn(rnorm(100),rnorm(100))
covQn <-function(x,y){
  leng.x <- length(x);
  leng.y <- length(y);
  alpha.qn <- robustbase::Qn(x);
  beta.qn <- robustbase::Qn(y);
  if(leng.x == leng.y){
    if(identical(x,y)){
      cov.qn <- robustbase::Qn(x)^2;
      return(cov.qn);
    } else {
      cov.qn <- 0;
      x.alpha <- x/alpha.qn;
      y.beta <- y/beta.qn;
      Qn1 <- robustbase::Qn(x.alpha + y.beta);
      Qn2 <- robustbase::Qn(x.alpha - y.beta);
      cov.qn <- ((alpha.qn*beta.qn)/4)*( (Qn1*Qn1) - (Qn2*Qn2));
      return(cov.qn);
    }
  } else {
    stop("x and y are unequal sizes");}
}
