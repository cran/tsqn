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
PerQn <- function(x, window = "truncated", bandw.rob = 0.7){
  n <- length(x)
  g <- n - 1L
  if(g <= 1L) {
    return(numeric(0))
  }

  w <- (2 * pi * seq_len(g)) / n
  m <- trunc(n^bandw.rob)
  n.weights <- n - 2L
  if(is.null(window)) {
    pw <- rep(1, n.weights)
  } else { # (window == "truncated")
    pw <- as.numeric(seq_len(n.weights) < m)
  }

  cov.aux <- robacf(x, lag.max = n, type = "covariance", plot = FALSE)$acf
  cov.vec <- cov.aux[, 1L, 1L]
  cov.x0 <- cov.vec[1L]
  cov.x <- cov.vec[2:(n - 1L)]
  weighted.cov <- cov.x * pw
  lag.idx <- seq_len(n.weights)

  per <- numeric(g)
  for(i in seq_len(g)){
    per[i] <- (1 / (2 * pi)) * (cov.x0 + 2 * sum(weighted.cov * cos(w[i] * lag.idx)))
  }

  per[seq_len(g - 1L)]
}
