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
PerioMrob <- function(series){
  n <- length(series)
  g <- n - 1L
  periorob <- numeric(g)
  if(g <= 0L) {
    return(periorob)
  }

  idx <- seq_len(n)
  freq.scale <- 2 * pi / n
  fft.scale <- sqrt(n / (8 * pi))

  for(j in seq_len(g)){
    wj <- freq.scale * j
    mx <- cbind(cos(wj * idx), sin(wj * idx))
    fitrob <- MASS::rlm(series ~ mx - 1, method = "M", psi = MASS::psi.huber)
    fft.j <- fft.scale * complex(real = fitrob$coef[1], imaginary = -fitrob$coef[2])
    periorob[j] <- Mod(fft.j)^2
  }

  periorob
}
