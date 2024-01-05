# UG package v1.0
# Last revised: June 16, 2022, rv03

# Compute weighted quantile estimator using Harrell-Davis method;
# citation: https://aakinshin.net/posts/weighted-quantiles/

# Revision history: rv03 uses Kish's effective sample size;
#  rv02 first removes any NaN, Inf, or NA values



# Weighted generic quantile estimator



wquantile_generic <- function(x, probs, cdf_gen, weights = NA) {
  
  bad.x <- x %in% c(NaN,NA_real_,Inf,-Inf)
  x <- x[!bad.x]; weights <- weights[!bad.x]
  
  n <- length(x)
  if (any(is.na(weights)))
    weights <- rep(1 / n, n)
  nw <- sum(weights)^2 / sum(weights^2)  # Kish's effective sample size
  
  indexes <- order(x)
  x <- x[indexes]
  weights <- weights[indexes]
  
  weights <- weights / sum(weights)
  cdf.probs <- cumsum(c(0, weights))
  
  sapply(probs, function(p) {
    cdf <- cdf_gen(nw, p)
    q <- cdf(cdf.probs)
    w <- tail(q, -1) - head(q, -1)
    sum(w * x)
  })
}
