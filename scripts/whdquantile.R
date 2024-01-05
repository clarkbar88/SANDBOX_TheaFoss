# UG package v1.0
# Last revised: February 25, 2021, rv02

# Compute weighted quantile estimator using Harrell-Davis method;
# citation: https://aakinshin.net/posts/weighted-quantiles/

# Revisions: none

## source('wquantile_generic.R')

whdquantile <- function(x, probs, weights = NA) {
  
  cdf_gen <- function(n, p) return(function(cdf.probs) {
    pbeta(cdf.probs, (n + 1) * p, (n + 1) * (1 - p))
  })
  wquantile_generic(x, probs, cdf_gen, weights)
}
