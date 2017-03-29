\name{TimeSeriesQn}
\alias{tsqn}
\docType{package}
\title{Applications of the Qn estimator to time series (univariate and multivariate)}
\description{This package contains applications of the Qn estimator of Rousseeuw and Croux(1993) to univariate and multivariate Time Series in time and frequency domains. More especifically, the robust estimation of autocorrelation or autocovariance matrix functions from Ma and Genton (2000,2001) and Cotta et. al. (2017). The robust periodogram of
    Molinares et. al. (2009) and the M-Periodogram of Reisen et. al. (2017). The robust GPH estimator of d considering
    robust periodogram approach. }

\details{
\tabular{ll}{
Package: \tab TimeSeriesQn\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2017-03-08\cr
License: \tab GPL (>= 2)\cr
}
}

\author{
Higor Cotta, Valderio Reisen, Pascal Bondon and Céline Lévy-Leduc


Maintainer: Higor Cotta \email{cotta.higor@gmail.com}
}
\references{
  H. Cotta, V. Reisen, P. Bondon, W. Stummer (2017):
  Robust Estimation of Covariance and Correlation Functions of a Stationary Multivariate Process. \emph{To appear in 2017 25th European Signal Processing Conference (EUSIPCO 2017)}

  V. Reisen, C. Lévy-Leduc, M. Taqqu (2017):
  An M-estimator for the long-memory parameter. \emph{To appear in Journal of Statistical Planning and Inference}.

  Y. Ma, M. Genton (2000):
  Highly robust estimation of the autocovariance function. \emph{Journal of Time Series Analysis}, 21, 663--684.

  Y. Ma, M. Genton (2001):
  Highly robust estimation of dispersion matrices. \emph{Journal of Multivariate Analysis}, 78, 11--36.

  P. J. Rousseuw, C. Croux (1993):
  Alternatives to the median absolute deviation.,\emph{Journal of the American Statistical Association}, Taylor & Francis Group,  88,
  n. 424, p. 1273–-1283.}



\examples{
x <- rnorm(1000,10,1)
y <- rnorm(1000,100,10)
xy <- cbind(x,y)

covQn(x,y)
corQn(x,y)
covMatQn(xy)
corMatQn(xy)


robacf(x)

dataset <- cbind(fdeaths,mdeaths)
robacf(dataset)
robacf(dataset,type="covariance",lag.max=10)

PerQn(ldeaths)

library(fracdiff)
simseries <- fracdiff.sim(1500, d = 0.3)
GPH_estimate(simseries$series,method="GPH")$d
GPH_estimate(simseries$series,method="GPH-Qn")$d

}
