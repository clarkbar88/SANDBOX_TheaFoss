# Author: Kirk Cameron, Ph.D., MacStat Consulting, Ltd.
# Last revised: November 3, 2022, rv02

# Purpose: generate Monte Carlo vector of random variates for use with censored data

# Revisions: none

mc_vector <- function(conc,nd,mc.model=c('beta','uniform','triangle'),shape1=1.5,shape2=1.5) {

  N <- length(conc)
  mc.model <- match.arg(mc.model)
  if (mc.model == 'beta') {
    rv <- rbeta(N,shape1,shape2)
  } else if (mc.model %in% c('uniform','triangle')) {
    rv <- runif(N)
    if (mc.model == 'triangle') rv <- sqrt(rv/2)*(rv <= 0.5) + (1-sqrt((1-rv)/2))*(rv > 0.5)
  }
  rv <- ifelse(!nd,conc,conc*rv)
  rv
}