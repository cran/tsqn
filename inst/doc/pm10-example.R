## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(tsqn)
data("pm10")

dim(pm10)
head(pm10)

## -----------------------------------------------------------------------------
pm10_subset <- as.matrix(pm10[1:365, ])

qn_cor <- corMatQn(pm10_subset)
qn_cov <- covMatQn(pm10_subset)

round(qn_cor, 3)
round(qn_cov, 1)

## -----------------------------------------------------------------------------
vix <- pm10_subset[, "VixCentro"]

acf_qn <- robacf(vix, lag.max = 24, type = "correlation", plot = FALSE)
head(acf_qn$acf[, 1, 1], 10)

per_qn <- PerQn(vix)
length(per_qn)
head(per_qn, 10)

GPH_estimate(vix, method = "GPH-Qn")

## ----fig.width=7, fig.height=4------------------------------------------------
tsqn:::plot.robacf(acf_qn, main = "PM10 Robust ACF (VixCentro)")

## ----fig.width=7, fig.height=4------------------------------------------------
stats::acf(vix, lag.max = 24, main = "PM10 Standard ACF (VixCentro)")

