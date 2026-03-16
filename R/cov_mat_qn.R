#' Robust covariance matrix
#'
#' Computes the robust covariance matrix of the matrix \code{x} proposed by Ma and Genton (2001) using the robust scale Qn of Rousseeuw and Croux (1993).
#' @param x a numeric matrix
#' @return a numeric matrix
#' @references Ma, Y. and Genton, M. G. (2001) Highly robust estimation of dispersion matrices. \emph{Journal of Multivariate Analysis}, \bold{78}, 11--36.
#' @references Rousseeuw, P. J. and Croux, C. (1993) Alternatives to the median absolute deviation. \emph{Journal of the American Statistical Association}, \bold{88}, 1273--1283.
#' @export
#' @examples
#' dataset <- cbind(rnorm(100),rnorm(100))
#' covMatQn(dataset)
covMatQn <- function(x){
  x <- as.matrix(x)
  if(!is.numeric(x)) {
    stop("'x' must be numeric")
  }

  n <- ncol(x)
  qn.scale <- apply(x, 2L, robustbase::Qn)
  cov.Mat.Qn <- matrix(0, nrow = n, ncol = n)
  diag(cov.Mat.Qn) <- qn.scale * qn.scale

  if(n <= 1L) {
    return(cov.Mat.Qn)
  }

  for(i in seq_len(n - 1L)){
    alpha.qn <- qn.scale[i]
    xi.scaled <- x[, i] / alpha.qn
    for(j in (i + 1L):n){
      beta.qn <- qn.scale[j]
      yj.scaled <- x[, j] / beta.qn
      qn1 <- robustbase::Qn(xi.scaled + yj.scaled)
      qn2 <- robustbase::Qn(xi.scaled - yj.scaled)
      value <- ((alpha.qn * beta.qn) / 4) * ((qn1 * qn1) - (qn2 * qn2))
      cov.Mat.Qn[i, j] <- value
      cov.Mat.Qn[j, i] <- value
    }
  }

  cov.Mat.Qn
}
