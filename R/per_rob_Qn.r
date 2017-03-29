#' Robust periodogram based on the Robust ACF
#'
#' Computes the robust pseudo-periodogram of Molinares et al (2009) based on the robust ACF by Ma and Genton (2000).
#' @param x univariate time series
#' @param window character string giving the type of the window. Allowed values are "truncated" (the default) or "\code{NULL}".
#' @param bandw.rob is a numeric value giving the truncation point.
#' @return a numeric vector containing the values of the robust periodogram proposed by Molinares (2009).
#' @author Valderio Reisen and Higor Cotta
#' @references Molinares, F. F. and Reisen, V. A., and Cribari-Neto, F. (2009) Robust estimation in long-memory processes under additive outliers. \emph{Journal of Statistical Planning and Inference}, \bold{139}, 2511--2525.
#' @references Ma, Y. and Genton, M. G. (2000) Highly robust estimation of the autocovariance function. \emph{Journal of Time Series Analysis}, \bold{21}, 663--684.
#' @export
#' @examples
#' PerQn(ldeaths)
PerQn <- function(x,window="truncated",bandw.rob=0.7){
  n <- length(x)
  kk <- 1:(n-1)
  g <- (n-1)
  j <- 1:g
  w <- (2*pi*j)/n
  M <- trunc(n^bandw.rob)
  pw <- numeric(n-2)
  if (is.null(window)){
    pw <- rep(1,n-2)}
  else { #(window=="truncated")
    for(i in 1:(n-2)){
      pw[i]<-ifelse(i<M,1.0,0.0)
    }
  }
  cov.aux <- robacf(x,lag.max = n,type="covariance",plot=FALSE)$acf
  cov.x0 <- cov.aux[,,1][1]
  cov.x <- cov.aux[,,1][2:(n-1)]
  per<- matrix(0, g, 1)
  for(i in 1:g){
    per[i] <- (1/(2*pi))*(cov.x0 + 2*sum(cov.x*pw*cos(w[i]*(1:(n-2)))))
  }
  return(per[1:(g-1)])
}
