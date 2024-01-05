# Author: Kirk Cameron, PhD, MacStat Consulting, Ltd
# Last revised: November 08, 2022, rv02

# Purpose: Use KM to estimate plotting positions (aka percentage points) for
# left-censored data; KM = Kaplan-Meier;
# defaults to Weibull plotting positions if no censoring is present

# Note: x is vector of observations and/or detection/reporting limits;
# nd is vector of censoring statuses: 1 = non-detect, 0 = detected/reported

# Revisions: rv02 changes name from <km_pp> to <pp_km>


pp_km <- function(x,nd,id,wt=NULL,eps=.Machine$double.eps) {
  library(tidyverse)
  
# first remove any missing data; normalize weights
  if (is.null(wt)) {
    f <- tibble(x,nd,id) %>% mutate(wt=1)
  } else {
    f <- tibble(x,nd,id,wt)
  }
  f <- f %>% filter(!is.na(x)) %>% mutate(n=n(),wt=n*wt/sum(wt))
  
# default cases of no censoring or complete censoring
  n <- nrow(f)
  if ((sum(nd)==0) || (sum(1-nd)==0)) {
    f.pp <- f %>% arrange(x) %>% mutate(pp=cumsum(wt)/(n+1))
    return(f.pp)
  }
  

# table by nd status; determine unique x values and number at risk; build KM CDF
  km.tab <- f %>% mutate(xx=x-nd*eps) %>% group_by(xx,nd) %>% summarise(km.lev=x[1],nn=n(),sw=sum(wt),.groups = 'drop') %>% mutate(km.rsk=cumsum(sw),tab1=case_when(nd==0~sw,nd==1~0),km.cdf=rev(c(1,cumprod(1 - rev(tab1[-1])/rev(km.rsk[-1])))),km.surv=1 - km.cdf)
  km.tab <- km.tab %>% select(-xx) %>% rename(sum.wt=sw) %>% relocate(km.lev,.before = nd)


# plotting positions for detects
  tmp <- km.tab %>% filter(nd == 0) %>% mutate(km.cdf=km.cdf*(n/(n+1))) %>% select(km.lev,nd,km.cdf) %>% rename(x=km.lev,pp=km.cdf)
  f.d <- f %>% filter(nd == 0) %>% left_join(.,tmp,by = c('x','nd'))
  
# plotting positions for NDs
  tmp <- km.tab %>% filter(nd == 1) %>% mutate(km.cdf=km.cdf*(n/(n+1)),km.inc=km.cdf/(nn+1)) %>% select(km.lev,nd,km.inc) %>% rename(x=km.lev)
  f.nd <- f %>% filter(nd == 1) %>% arrange(x) %>% left_join(.,tmp,by = c('x','nd')) %>% group_by(x) %>% mutate(pp=cumsum(km.inc)) %>% ungroup() %>% select(-km.inc)


# return tibble with concatenated data and plotting positions
  f.pp <- rbind(f.d,f.nd) %>% arrange(pp,x,desc(nd))
  f.pp
}
