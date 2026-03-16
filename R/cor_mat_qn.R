#' Robust correlation matrix
#'
#' Computes the robust correlation matrix of the matrix \code{x} proposed by Ma and Genton (2001) using the robust scale Qn of Rousseeuw and Croux (1993).
#' @param x a numeric matrix
#' @return a numeric matrix
#' @references Ma, Y. and Genton, M. G. (2001) Highly robust estimation of dispersion matrices. \emph{Journal of Multivariate Analysis}, \bold{78}, 11--36.
#' @references Rousseeuw, P. J. and Croux, C. (1993) Alternatives to the median absolute deviation. \emph{Journal of the American Statistical Association}, \bold{88}, 1273--1283.
#' @export
#' @examples
#' dataset <- cbind(rnorm(100),rnorm(100))
#' corMatQn(dataset)
corMatQn <- function(x){
  x <- as.matrix(x)
  if(!is.numeric(x)) {
    stop("'x' must be numeric")
  }

  n <- ncol(x)
  cor.Mat.Qn <- diag(1, nrow = n, ncol = n)
  if(n <= 1L) {
    return(cor.Mat.Qn)
  }

  qn.scale <- apply(x, 2L, robustbase::Qn)

  for(i in seq_len(n - 1L)){
    xi.scaled <- x[, i] / qn.scale[i]
    for(j in (i + 1L):n){
      yj.scaled <- x[, j] / qn.scale[j]
      qn1 <- robustbase::Qn(xi.scaled + yj.scaled)
      qn2 <- robustbase::Qn(xi.scaled - yj.scaled)
      qn1.sq <- qn1 * qn1
      qn2.sq <- qn2 * qn2
      value <- (qn1.sq - qn2.sq) / (qn1.sq + qn2.sq)
      cor.Mat.Qn[i, j] <- value
      cor.Mat.Qn[j, i] <- value
    }
  }

  cor.Mat.Qn
}
