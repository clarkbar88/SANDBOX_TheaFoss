# Author: Kirk Cameron, Ph.D., MacStat Consulting, Ltd.
# Last revised: October 12, 2023, rv02

# Purpose: generate Monte Carlo vector of random variates for use with censored data

# Revisions: rv02 minor update

mc_vector <- function(y,nd,mc.model=c('beta','uniform','triangle'),shape1=1.5,shape2=1.5) {
  
  N <- length(y)
  mc.model <- match.arg(mc.model)
  if (mc.model == 'beta') {
    rv <- rbeta(N,shape1,shape2)
  } else if (mc.model %in% c('uniform','triangle')) {
    rv <- runif(N)
    if (mc.model == 'triangle') rv <- sqrt(rv/2)*(rv <= 0.5) + (1-sqrt((1-rv)/2))*(rv > 0.5)
  }
  if_else(!nd,y,y*rv)
}
