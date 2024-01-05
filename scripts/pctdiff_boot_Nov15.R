# Author: Kirk Cameron, PhD, MacStat Consulting, Ltd
# Last Update: December 16, 2022, rv01

# Purpose: Compute percent change bootstrap CI, with or without NDs

# Revisions: none


pctdiff_boot <- function(f,y=conc,nd=nd,clev=0.95,side=c('upr','lwr','both'),impute=c('KM','MC'),ns=50,nb=500,...) {
  library(tidyverse)
  library(magrittr)
  library(boot)
  library(purrr)
  library(rlang)

  y <- ensym(y); nd <- ensym(nd)
  
# remove any rows with missing values in y; check for type of confidence interval
  f <- f %>% filter(!is.na(!!y)) %>% mutate(conc.half=case_when(!!nd == 0~!!y,!!nd == 1~!!y/2))
  side <- match.arg(side)
  impute <- match.arg(impute)

  n <- nrow(f)
  nd.pct <- with(f,inject(round(100*mean(!!nd,na.rm = T),1)))
  any.nd <- with(f,inject(any(!!nd == 1)))

# if <wt> does not exist, set to 1
  wtd <- with(f,exists('wt'))
  if (!wtd) f <- f %>% mutate(wt=1)
  
# case of insufficient data
  if (n < 4) return(tibble(model='BOOT',side,n,nd.pct,clev,impute='NONE',pctdiff=NA_real_,lcl=NA_real_,ucl=NA_real_))

  ix <- 1:n
  
# convert dates to numeric if necessary
  f <- f %>% mutate(x=as.integer(date))
  xrng <- f %>% pull(x) %>% range() %>% as_tibble() %>% rename(x=value)
  
  boot_pctdiff <- function(f,ix,y) {
    y <- ensym(y)
    
    td <- f[ix,] %>% arrange(x)
    lm0 <- inject(lm(!!y ~ x,data=td))
    fit0 <- predict(lm0,newdata=xrng)
    pctdiff <- 100*((fit0[2] - fit0[1])/fit0[1])
    pctdiff
  }
  
  
  assemble_bootci <- function(boot.out) {
    
    if (side == 'upr') {
      ci0 <- whdquantile(boot.out$t,p=c(0.5,clev))
      ci0 <- tibble(pctdiff=ci0[1],lcl=NA_real_,ucl=ci0[2])
    } else if (side == 'lwr') {
      ci0 <- whdquantile(boot.out$t,p=c(0.5,1-clev))
      ci0 <- tibble(pctdiff=ci0[1],lcl=ci0[2],ucl=NA_real_)
    } else if (side == 'both') {
      tlev <- (1 + clev)/2
      ci0 <- whdquantile(boot.out$t,p=c(0.5,1-tlev,tlev))
      ci0 <- tibble(pctdiff=ci0[1],lcl=ci0[2],ucl=ci0[3])
    }
    ci0
  }
  

  
# compute pct difference bootstrap CI estimate
# first case: all detects
  if (!any.nd) {
    
    impute <- 'NONE'
    boot.out <- inject(boot::boot(f,boot_pctdiff,R=nb,stype='i',y=!!y,parallel = 'multicore',ncpus = 10))
    tmp <- assemble_bootci(boot.out)
    
# second case: any NDs imputed using Kaplan-Meier and assuming Normal model
  } else if (any.nd) {
##    yv <- f %>% pull({{y}}); ndv <- f %>% pull({{nd}})
    
    if (impute == 'KM') {
      yt <- with(f,inject(impute_km(!!y,!!nd,id,tlab='Normal',wt))$f.pp$xhat)
      f <- f %>% mutate(y.r=yt)
      boot.out <- boot::boot(f,boot_pctdiff,R=nb,stype='i',y=y.r,parallel = 'multicore',ncpus = 10)
      tmp <- assemble_bootci(boot.out)
      
# third case: MC imputation of NDs
    } else if (impute == 'MC') {
    
      tmp.out <- map_dfr(1:ns,.f=function(.x,xf,...) {
        xf <- xf %>% mutate(row=.x,y.r=mc_vector(!!y,!!nd,...))
        boot.out <- boot::boot(xf,boot_pctdiff,R=nb,stype='i',y=y.r,parallel = 'multicore',ncpus = 10)
        assemble_bootci(boot.out)
      },xf=f,...)


      if (side == 'upr') {
        tmp <- tmp.out %>% summarise(across(c(pctdiff,ucl),~ mean(.x,na.rm=T)),lcl=NA_real_) %>% relocate(lcl,.before = ucl)
      } else if (side == 'lwr') {
        tmp <- tmp.out %>% summarise(across(pctdiff:lcl,~ mean(.x,na.rm=T)),ucl=NA_real_)
      } else if (side == 'both') {
        tmp <- tmp.out %>% summarise(across(pctdiff:ucl,~ mean(.x,na.rm=T)))
      }
    }
  }


# output estimated CI
  out <- tibble(model='BOOT',side,n,nd.pct,clev,impute,tmp)
  out
}

