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
robacf <- function(x, lag.max = NULL, type = c("correlation", "covariance"), plot = TRUE, na.action = na.fail, demean = TRUE, ...){
  type <- match.arg(type)
  series <- deparse(substitute(x))
  x <- na.action(as.ts(x))
  x.freq <- frequency(x)
  x <- as.matrix(x)

  if(!is.numeric(x)) {
    stop("'x' must be numeric")
  }

  sampleT <- as.integer(nrow(x))
  nser <- as.integer(ncol(x))
  if(is.na(sampleT) || is.na(nser)) {
    stop("'sampleT' and 'nser' must be integer")
  }

  if(is.null(lag.max)) {
    lag.max <- floor(10 * (log10(sampleT) - log10(nser)))
  }
  lag.max <- as.integer(min(lag.max, sampleT - 1L))
  if(is.na(lag.max) || lag.max < 0L) {
    stop("'lag.max' must be at least 0")
  }

  if(demean) {
    x <- sweep(x, 2L, colMeans(x, na.rm = TRUE), check.margin = FALSE)
  }

  lag <- matrix(1, nser, nser)
  lag[lower.tri(lag)] <- -1
  acf.Qn <- array(1, c(lag.max, nser, nser))
  lag.seq <- seq_len(lag.max)
  is.correlation <- identical(type, "correlation")

  if(lag.max > 0L) {
    if(nser == 1L) { # Univariate
      x.vec <- x[, 1L]
      for(h in lag.seq){
        idx.end <- sampleT - h + 1L
        u <- x.vec[h:sampleT]
        v <- x.vec[seq_len(idx.end)]
        q.plus <- robustbase::Qn(u + v)
        q.min <- robustbase::Qn(u - v)
        q.plus.sq <- q.plus * q.plus
        q.min.sq <- q.min * q.min
        if(is.correlation) {
          acf.Qn[h, 1L, 1L] <- (q.plus.sq - q.min.sq) / (q.plus.sq + q.min.sq)
        } else {
          acf.Qn[h, 1L, 1L] <- 0.25 * (q.plus.sq - q.min.sq)
        }
      }
    } else { # Multivariate
      qn.scale <- apply(x, 2L, robustbase::Qn)
      x.scaled <- sweep(x, 2L, qn.scale, "/", check.margin = FALSE)

      for(i in seq_len(nser)){
        xi <- x[, i]
        xi.scaled <- x.scaled[, i]
        for(j in seq_len(nser)){
          yj <- x[, j]
          yj.scaled <- x.scaled[, j]
          use.scaled <- i != j
          if(!is.correlation) {
            pair.scale <- if(use.scaled) (qn.scale[i] * qn.scale[j]) / 4 else 0.25
          }

          for(h in lag.seq){
            idx.end <- sampleT - h + 1L
            if(use.scaled) {
              u <- xi.scaled[h:sampleT]
              v <- yj.scaled[seq_len(idx.end)]
            } else {
              u <- xi[h:sampleT]
              v <- yj[seq_len(idx.end)]
            }

            q.plus <- robustbase::Qn(u + v)
            q.min <- robustbase::Qn(u - v)
            q.plus.sq <- q.plus * q.plus
            q.min.sq <- q.min * q.min

            if(is.correlation) {
              acf.Qn[h, j, i] <- (q.plus.sq - q.min.sq) / (q.plus.sq + q.min.sq)
            } else {
              acf.Qn[h, j, i] <- pair.scale * (q.plus.sq - q.min.sq)
            }
          }
        }
      }
    }
  }

  lag.values <- if(lag.max > 0L) 0:(lag.max - 1L) else numeric(0)
  lag <- outer(lag.values, lag / x.freq)
  acf.out <- structure(list(acf = acf.Qn, type = type, n.used = sampleT, lag = lag, series = series, snames = colnames(x)),
                       class = "robacf")

  if(plot) {
    plot.robacf(acf.out)
    invisible(acf.out)
  } else {
    acf.out
  }
}
