#' Robust autocorrelation or autocovariance function estimation
#'
#' This function computer and plots(by default) the robust estimates of the autocovariance or the autocorrelation function
#' based on the Qn.
#'
#' @param x a numeric vector or matrix.
#' @param lag.max maximum lag at which to calculate the acf. Default is 10*log10(N/m) where
#' N is the number of observations and m the number of series. Will be automatically limited
#' to one less than the number of observations in the series.
#' @param type character string giving the type of acf to be computed. Allowed values are "correlation" (the default) or "covariance".
#' Accepts parcial names.
#' @param plot logical. If TRUE (the default) the acf is plotted.
#' @param na.action function to be called to handle missing values. na.pass can be used.
#' @param demean logical. Should the covariances be about the sample means?
#' @param ... further arguments to be passed to plot.acf.
#' @return An object of class "robacf", which is a list with the following elements:
#' @return \code{lag} A three dimensional array containing the lags at which the acf is estimated.
#' @return \code{acf} An array with the same dimensions as lag containing the estimated acf.
#' @return \code{type} The type of correlation (same as the type argument).
#' @return \code{n.used} The number of observations in the time series.
#' @return \code{series} The name of the series x.
#' @return \code{snames} The series names for a multivariate time series.
#' @return The result is returned invisibly if plot is TRUE.
#' @author Higor Cotta, Valderio Reisen and Pascal Bondon
#' @references Cotta, H. and Reisen, V. A. and Bondon, P. and Stummer, W. (2017) Robust Estimation of Covariance and Correlation Functions of a Stationary Multivariate Process. \emph{To appear in 2017 25th European Signal Processing Conference (EUSIPCO 2017).}
#' @references Ma, Y. and Genton, M. G. (2000) Highly robust estimation of the autocovariance function. \emph{Journal of Time Series Analysis}, \bold{21}, 663--684.
#' @references Ma, Y. and Genton, M. G. (2001) Highly robust estimation of dispersion matrices. \emph{Journal of Multivariate Analysis}, \bold{78}, 11--36.
#' @references Rousseeuw, P. J. and Croux, C. (1993) Alternatives to the median absolute deviation. \emph{Journal of the American Statistical Association}, \bold{88}, 1273--1283.
#' @export
#' @import robustbase
#' @import stats
#' @examples
#' data.set <- cbind(fdeaths,mdeaths)
#' robacf(data.set)
#' robacf(data.set,type="covariance",lag.max=10)
robacf <- function(x, lag.max = NULL, type=c("correlation", "covariance"), plot=TRUE, na.action = na.fail, demean = TRUE, ...){
  type <- match.arg(type)
  series <- deparse(substitute(x))
  x <- na.action(as.ts(x))
  x.freq <- frequency(x)
  x <- as.matrix(x)
  if(!is.numeric(x))
    stop("'x' must be numeric")
  sampleT <- as.integer(nrow(x))
  nser <- as.integer(ncol(x))
  if(is.na(sampleT) || is.na(nser))
    stop("'sampleT' and 'nser' must be integer")
  if (is.null(lag.max))
    lag.max <- floor(10 * (log10(sampleT) - log10(nser)))
  lag.max <- as.integer(min(lag.max, sampleT - 1L))
  if (is.na(lag.max) || lag.max < 0)
    stop("'lag.max' must be at least 0")
  if(demean)
    x <- sweep(x, 2, colMeans(x, na.rm = TRUE), check.margin=FALSE)
  lag <- matrix(1, nser, nser)
  lag[lower.tri(lag)] <- -1
  acf.Qn <- matrix(1,lag.max,1)
  acf.Qn <- array(1, c(lag.max,nser,nser))
  if(nser==1L){#Univariate
    if(type=="correlation"){
      for(h in 1:(lag.max)){
        U <- x[h:sampleT]
        V <- x[1:(sampleT-(h-1))]
        Q.plus <- robustbase::Qn(U+V)^2
        Q.min <- robustbase::Qn(U-V)^2
        acf.Qn[h]<-(Q.plus-Q.min)/(Q.plus+Q.min)
      }
    }
    else{
      for(h in 1:(lag.max)){
        U <- x[h:sampleT]
        V <- x[1:(sampleT-(h-1))]
        Q.plus<- robustbase::Qn(U+V)^2
        Q.min <- robustbase::Qn(U-V)^2
        acf.Qn[h] <- (1/4)*(Q.plus-Q.min)
      }
    }
  }
    else{#Multivariate
      if(type=="correlation"){
        for(i in 1:nser){
          for(j in 1:nser){
            for(h in 1:(lag.max)){
              if(i==j){
                alpha.qn <- 1
                beta.qn <- 1
              }
              else{
                alpha.qn <- robustbase::Qn(x[,i])
                beta.qn <- robustbase::Qn(x[,j])
              }
              U <- x[h:sampleT,i]/alpha.qn
              V <- x[1:(sampleT-(h-1)),j]/beta.qn
              Q.plus <- robustbase::Qn(U+V)^2
              Q.min <- robustbase::Qn(U-V)^2
              acf.Qn[h,j,i] <- (Q.plus-Q.min)/(Q.plus+Q.min)
            }
          }
        }
      }
      else{
        for(i in 1:nser){
          for(j in 1:nser){
            for(h in 1:(lag.max)){
              if(i==j){
                alpha.qn <- 1
                beta.qn <- 1
              }
              else{
                alpha.qn <- robustbase::Qn(x[,i])
                beta.qn <- robustbase::Qn(x[,j])
              }
              U <- x[h:sampleT,i]/alpha.qn
              V <- x[1:(sampleT-(h-1)),j]/beta.qn
              Q.plus <- robustbase::Qn(U+V)^2
              Q.min <- robustbase::Qn(U-V)^2
              acf.Qn[h,j,i] <- (alpha.qn*beta.qn/4)*(Q.plus-Q.min)
            }
          }
        }
      }

    }
    lag <- outer(0:(lag.max-1), lag/x.freq)
  acf.out <- structure(list(acf = acf.Qn, type = type, n.used = sampleT,lag = lag, series = series, snames = colnames(x)),
                       class = "robacf")
  if (plot) {
    plot.robacf(acf.out)
    invisible(acf.out)
  }else acf.out
}
