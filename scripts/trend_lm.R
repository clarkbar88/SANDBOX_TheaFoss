# Author: Kirk Cameron, PhD, MacStat Consulting, Ltd
# Last Update: December 13, 2022, rv02

# Purpose: Compute linear model trend, with or without NDs

# Revisions: rv02 renames from <lm_trend>; updates code


trend_lm <- function(f,y=conc,nd=nd,ns=100,...) {
  library(tidyverse)
  library(rlang)
  library(purrr)
  library(broom)

  ys <- ensym(y); nds <- ensym(nd)
  
# remove any rows with missing values in date or conc; set seed
  n <- nrow(f); nd.pct <- with(f,round(100*inject(mean(!!nds)),1))
  any.nd <- with(f,inject(any(!!nds == 1)))
  f <- f %>% filter(!(is.na(date) | is.na(!!ys)))

  set.seed(60021151)
  
# case of insufficient data
  if (n < 3) return(tibble(n,nd.pct,intercept=NA_real_,slope=NA_real_,pval=NA_real_))
  
# convert dates to numeric if necessary
  f <- f %>% mutate(x=as.integer(date))
  

# check for and remove any high influence pts or extreme standardized outliers
  if (!any.nd) {
    lm.tst <- inject(lm(!!ys~x,data = f))
  } else if (any.nd) {
    f <- f %>% mutate(y.half=case_when(!!nds == 1~!!ys/2,TRUE~!!ys))
    lm.tst <- lm(y.half~x,data = f)
  }
  pval.tst <- tidy(lm.tst)$p.value[2]
  # lm.tst.sum <- summary(lm.tst)
  # pval.tst <- lm.tst.sum$coeff[2,4]
  
# check residuals for high leverage and extreme standardized outliers;
# but first check for no variation/perfect fit and hence unusable regression summary
  if (is.nan(pval.tst)) {
    f <- f %>% mutate(lev.vec=0,resid=0,esi=0,hi.lev=F,esi.out=F)
  } else if (!is.nan(pval.tst)) {
    p <- 2; se <- summary(lm.tst)$sigma
    f <- f %>% mutate(lev.vec=influence(lm.tst)$hat,resid=lm.tst$residuals,esi=resid/(se*sqrt(1-lev.vec)),hi.lev=lev.vec > 3*p/n,esi.out=abs(esi) > 3)
  }
  
# remove high leverage pts and/or extreme outliers
  f <- f %>% filter(!(hi.lev | esi.out))

  
  
  
  
# if any censored values, set up loop to impute NDs, compute linear model, then average model results
  if (any.nd) {
    
    fit0.tmp <- map_dfr(1:ns,.f=function(.x,xf,...) {
      xf <- xf %>% mutate(y.r=mc_vector(!!ys,!!nds,...))
      lm0 <- lm(y.r~x,data=xf)
      tmp <- tidy(lm0)
      tibble(intercept=tmp$estimate[1],slope=tmp$estimate[2],pval=tmp$p.value[2])
    },xf=f,...)
    
    fit0 <- fit0.tmp %>% summarise(across(c(intercept:pval),~ mean(.x,na.rm=T)))


# case with no NDs
  } else if (!any.nd) {

    lm0 <- inject(lm(!!ys~x,data=f))
    tmp <- tidy(lm0)
    fit0 <- tibble(intercept=tmp$estimate[1],slope=tmp$estimate[2],pval=tmp$p.value[2])
  }

# output estimated trend line
  tibble(n,nd.pct,fit0)
}

