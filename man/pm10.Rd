\name{pm10}
\docType{data}
\alias{pm10}
\title{PM10 Concentrations at Monitoring Stations}
\description{
Hourly PM10 concentration measurements collected at monitoring stations.
}
\usage{
data(pm10)
}
\format{
A data frame with 1826 observations on the following 8 variables.
  \describe{
    \item{\code{Laranjeiras}}{PM10 concentration at Laranjeiras station.}
    \item{\code{Carapina}}{PM10 concentration at Carapina station.}
    \item{\code{Camburi}}{PM10 concentration at Camburi station.}
    \item{\code{Sua}}{PM10 concentration at Sua station.}
    \item{\code{VixCentro}}{PM10 concentration at VixCentro station.}
    \item{\code{Ibes}}{PM10 concentration at Ibes station.}
    \item{\code{VVCentro}}{PM10 concentration at VVCentro station.}
    \item{\code{Cariacica}}{PM10 concentration at Cariacica station.}
  }
}
\source{
Packaged from \code{inst/extdata/pm10data.csv}.
}
\examples{
data(pm10)
str(pm10)
summary(pm10$VixCentro)
}
\keyword{datasets}
