# Author: Kirk Cameron, PhD, MacStat Consulting, Ltd
# Last Update: November 11, 2022, rv06

# Purpose: Compute confidence interval, using parametric, nonparametric, or bootstrap approaches

# Revisions: rv06 updates code;
#  rv05 completes code; updates name to <ci_compute> from <compute_ci>;
#  rv04 adds options for NP or TBOOT intervals;
#  rv03 minor update for consistency;
#  rv02 alters side parameter


ci_compute <- function(f,y=conc,nd=nd,clev=0.95,tlab='Normal',side=c('upr','lwr','both'),meth=c('P','NP','BOOT','TBOOT'),impute=c('KM','MC'),...) {
  library(tidyverse)
  library(rlang)
  library(purrr)
  
##  source('ci_p.R')
##  source('ci_boot.R')
##  source('ci_np.R')
##  source('ci_tboot.R')
  
  y <- ensym(y); nd <- ensym(nd)

# check for type of CI and method
  side <- match.arg(side)
  meth <- match.arg(meth)
  impute <- match.arg(impute)
  
# remove any rows with missing conc values; return if insufficient data
  n <- nrow(f)
  f <- f %>% filter(!is.na(!!y))
  f <- f %>% mutate(conc.half=case_when(!!nd == 0~!!y,!!nd == 1~!!y/2))
  
  nd.pct <- with(f,inject(round(100*mean(!!nd),1))); ndet <- with(f,inject(sum(1 - !!nd)))
  
  
  if (n < 4) return(tibble(model=meth,n,nd.pct,clev,ave=mean(conc.half),lcl=NA_real_,ucl=NA_real_))


# force type to BOOT if too few detects or too many NDs or tlab is NULL (no model)
  if (ndet < 4 || nd.pct >= 80 || is.null(tlab)) meth <- 'BOOT'

# compute confidence interval by type
  if (meth == 'P') {
    out <- inject(ci_p(f,!!y,!!nd,clev,tlab,side,impute,...))
  } else if (meth == 'NP') {
    out <- inject(ci_np(f,!!y,!!nd,clev,side,...))
  } else if (meth == 'BOOT') {
    out <- inject(ci_boot(f,!!y,!!nd,clev,side,impute,...))
  } else if (meth == 'TBOOT') {
    out <- inject(ci_tboot(f,!!y,!!nd,clev,tlab,side,...))
  }

  out
}

