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
  n <- length(x[1,])
  cor.Mat.Qn <- matrix(data=NA,nrow=n,ncol=n)
  for(i in 1:n){
    for(j in 1:n){
      cor.Mat.Qn[i,j] <- corQn(x[ ,i],x[ ,j])
    }
  }
  return(cor.Mat.Qn)
}
