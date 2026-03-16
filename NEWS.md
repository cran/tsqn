# tsqn News

## tsqn 1.2.0 (2026-03-16)

### Data and Documentation

- Added built-in dataset `pm10` (1826 x 8) from `dataset/pm10data.csv`.
- Added dataset help page `?pm10`.
- Added a new vignette in `vignettes/pm10-example.Rmd` showing robust covariance/correlation, robust ACF, and robust GPH estimation with `pm10`.

## tsqn 1.1.0 (2026-03-16)

### Performance

- Optimized `corQn()` and `covQn()` to reduce temporary allocations and redundant computations.
- Reworked `corMatQn()` and `covMatQn()` to cache per-column `Qn` scales and compute symmetric matrix entries once.
- Optimized `robacf()` internals to avoid repeated scale recomputation in multivariate settings.
- Reduced loop overhead in `PerioMrob()`, `PerQn()`, and `GPH_estimate()` while preserving formulas and returned values.
