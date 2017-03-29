#' Robust M-periodogram
#'
#' This function computes the robust M-periodogram proposed by Reisen et al. (2017).
#' @param series univariate time series
#' @return a numeric vector containing the robust estimates of the spectral density
#' @author Valderio Reisen, Céline Lévy-Leduc and Higor Cotta.
#' @references Reisen, V. A. and Lévy-Leduc, C. and Taqqu, M. (2017) An M-estimator for the long-memory parameter. \emph{To appear in Journal of Statistical Planning and Inference}.
#' @references Geweke, J. and Porter-Hudak, S. (1983) The estimation and application of long memory time series models. \emph{Journal of Time Series Analysis}, \bold{4}, 221--238.
#' @export
#' @import MASS
#' @examples
#' PerioMrob(ldeaths)
PerioMrob=function(series){
  n <- length(series)
  periorob <- FFT <- NULL
  g <- n-1
  for(j in 1:g){
    X1 <- X2<-NULL
    w<-2*pi*j/n
    for(i in 1:n) {
      X1[i]<-cos(w*i); X2[i]<-sin(w*i)
    }
  MX <- cbind(X1,X2)
  fitrob <- MASS::rlm(series~MX-1, method = "M", psi = MASS::psi.huber)
  FFT[j] <- sqrt(n/(8*pi))*complex(real=fitrob$coef[1], imaginary=-fitrob$coef[2])
  periorob[j] <- Mod(FFT[j])^2
  }
  return(periorob)
}
