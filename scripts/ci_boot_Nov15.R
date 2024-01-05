# Author -- Kirk Cameron, PhD, MacStat Consulting, Ltd
# Last Update -- November 14, 2022, rv06

# Purpose: Compute percentile bootstrap CI, with or without NDs

# Revisions: rv06 installs flexible inputs;
#  rv05 improves ND imputation options; if KM option used, assumes Normal model;
#  rv04 switches to multicore, parallel bootstrap;
#  rv03 minor update for consistency;
#  rv02 adds median estimate of bootstrap distribution to output


ci_boot <- function(f,y=conc,nd=nd,clev=0.95,side=c('upr','lwr','both'),impute=c('KM','MC'),ns=50,nb=500,...) {
  library(tidyverse)
  library(magrittr)
  library(boot)
  library(purrr)
  library(rlang)
  
## source('whdquantile.R')

  
  side <- match.arg(side)
  impute <- match.arg(impute)

  y <- ensym(y); nd <- ensym(nd)
##  set.seed(rseed,kind = "L'Ecuyer-CMRG")
  
# remove any rows with missing values in conc; check for type of confidence band
  f <- f %>% filter(!is.na(!!y)) %>% mutate(conc.half=case_when(!!nd == 0~!!y,!!nd==1~!!y/2))
  

  n <- nrow(f)
  nd.pct <- with(f,inject(round(100*mean(!!nd,na.rm = T),1)))
  any.nd <- with(f,inject(any(!!nd == 1)))

# if <wt> does not exist, set to 1
  wtd <- with(f,exists('wt'))
  if (!wtd) f <- f %>% mutate(wt=1)
  
# case of insufficient data
  if (n < 4) return(tibble(model='BOOT',side,n,nd.pct,clev,attain=clev,stat='MEAN',impute='NONE',ave=mean(conc.half),lcl=NA_real_,ucl=NA_real_))

  ix <- 1:nrow(f)
  
  boot_ave <- function(f,ix,y) {
    y <- ensym(y)
    td <- f[ix,]
    tv <- td %>% pull(!!y)
    mean(tv,na.rm = T)
  }
  
  assemble_bootci <- function(boot.out) {
    
    if (side == 'upr') {
      ci0 <- whdquantile(boot.out$t,p=c(0.5,clev))
      ci0 <- tibble(ave=ci0[1],lcl=NA_real_,ucl=ci0[2])
    } else if (side == 'lwr') {
      ci0 <- whdquantile(boot.out$t,p=c(0.5,1-clev))
      ci0 <- tibble(ave=ci0[1],lcl=ci0[2],ucl=NA_real_)
    } else if (side == 'both') {
      tlev <- (1 + clev)/2
      ci0 <- whdquantile(boot.out$t,p=c(0.5,1-tlev,tlev))
      ci0 <- tibble(ave=ci0[1],lcl=ci0[2],ucl=ci0[3])
    }
    ci0
  }
  
# compute percentile bootstrap CI estimate
# first case: all detects
  if (!any.nd) {
    
    impute <- 'NONE'
    boot.out <- boot::boot(f,boot_ave,R=nb,stype='i',y=!!y,parallel = 'multicore',ncpus = 10)
    tmp <- assemble_bootci(boot.out)
    
# second case: any NDs imputed using Kaplan-Meier and assuming Normal model
  } else if (any.nd) {
##    yv <- f %>% pull(!!ys); ndv <- f %>% pull(!!nds)
    
    if (impute == 'KM') {
      yt <- with(f,inject(impute_km(!!y,!!nd,id,tlab='Normal',wt)$f.pp$xhat))
      f <- f %>% mutate(conc.r=yt)
      boot.out <- boot::boot(f,boot_ave,R=nb,stype='i',y=conc.r,parallel = 'multicore',ncpus = 10)
      tmp <- assemble_bootci(boot.out)
      
# third case: MC imputation of NDs
    } else if (impute == 'MC') {
    
      tmp.out <- map_dfr(1:ns,.f=function(.x,xf,...) {
        xf <- xf %>% mutate(row=.x,conc.r=mc_vector(!!y,!!nd,...))
        boot.out <- boot::boot(xf,boot_ave,R=nb,stype='i',y=conc.r,parallel = 'multicore',ncpus = 10)
        assemble_bootci(boot.out)
      },xf=f,...)


      if (side == 'upr') {
        tmp <- tmp.out %>% summarise(across(c(ave,ucl),~ mean(.x,na.rm=T)),lcl=NA_real_) %>% relocate(lcl,.before = ucl)
      } else if (side == 'lwr') {
        tmp <- tmp.out %>% summarise(across(ave:lcl,~ mean(.x,na.rm=T)),ucl=NA_real_)
      } else if (side == 'both') {
        tmp <- tmp.out %>% summarise(across(ave:ucl,~ mean(.x,na.rm=T)))
      }
    }
  }


# output estimated CI
  out <- tibble(model='BOOT',side,n,nd.pct,clev,attain=clev,stat='MEAN',impute,tmp)
  out
}

