#' Classical and Robust Geweke and Porter-Hudak (GPH) estimators for the long-memory parameter d of a long-range dependent stationary processes
#'
#' Estimate the fractional (or “memory”) parameter d of long-range dependent stationary processes by the method of Geweke and Porter-Hudak (GPH). (GPH-M) and (GPH-Qn) correspond to the estimators devised by Reisen et al. (2017) and Molinares (2009), respectively.
#'
#' @param series univariate time series
#' @param bandw.exp the bandwidth used in the regression equation
#' @param method character string giving the type of GPH to be computed. Allowed values are "\code{GPH}" (the default), "\code{GPH-M}" or "\code{GPH-Qn}".
#' @return \code{d} GPH estimate
#' @return \code{sd.as} asymptotic standard deviation
#' @return \code{sd.reg} standard error deviation
#' @author Valderio Reisen, Céline Lévy-Leduc and Higor Cotta.
#' @references Reisen, V. A. and Lévy-Leduc, C. and Taqqu, M. (2017) An M-estimator for the long-memory parameter. \emph{To appear in Journal of Statistical Planning and Inference}.
#' @references Molinares, F. F. and Reisen, V. A., and Cribari-Neto, F. (2009) Robust estimation in long-memory processes under additive outliers. \emph{Journal of Statistical Planning and Inference}, \bold{139}, 2511--2525.
#' #' @references Geweke, J. and Porter-Hudak, S. (1983) The estimation and application of long memory time series models. \emph{Journal of Time Series Analysis}, \bold{4}, 221--238.
#' @export
#' @import fracdiff
#' @import stats
#' @examples
#' library(fracdiff)
#' simseries <- fracdiff.sim(1500, d = 0.3)
#' GPH_estimate(simseries$series,method="GPH")$d
#' \dontrun{
#' GPH_estimate(simseries$series,method="GPH-Qn")$d
#' GPH_estimate(simseries$series,method="GPH-M")$d
#' }
GPH_estimate <- function (series,bandw.exp = 0.7,method="GPH"){
  # There are 3 methods to compute GPH fractional estimates:
  # Method = GPH ( classical Periodogram), GPH-M ( M-spectral) GPH-Qn ( ACF robust Qn)
  n <- length(series)
  gn <- trunc(n^bandw.exp)
  if(method=='GPH') {
    # USE the periodogram classical- GPH- estimates.
    n <- length(series)
    x <- unclass(series) - mean(series)
    FFT <- fft(x)/sqrt(n)
    periodogramFFT <-((Mod(FFT))^2)/(2*pi)
    spectralest<-periodogramFFT[2:n]
  }
  else if(method =='GPH-M'){
       spectralest<-PerioMrob(series)
  }
  else {
    method <- "Qn"
    spectralest<-PerQn(series)
  }
  # Next steps: Estimation of d by regression method.
  contgn <- 0.0
  y.reg <-NULL
  x.reg <- NULL
  for(j in 1:gn) {
    w <- 2 * pi * j/n
    if(spectralest[j] > 0) {
      y.reg[j] <- log(spectralest[j]/(2 * pi))
      x.reg[j] <- 2 * log(2 * sin(w/2))
      contgn <- contgn +1
    }}
  fit <- stats::lm(y.reg ~ x.reg)
  d.GPH <- coef(fit)[2]
  names(d.GPH) <- NULL
  x.r2 <- sum((x.reg - mean(x.reg))^2)
  var.reg <- sum(resid(fit)^2)/((contgn - 1) * x.r2)
  results <- list(method=method,d = -d.GPH, sd.reg = sqrt(var.reg),bandw.exp)
}
