# Author: Kirk Cameron, PhD, MacStat Consulting, Ltd
# Last Update: November 15, 2022, rv10

# Purpose: compute linear model simultaneous confidence band around linear regression line,
# with or without NDs; default to confidence interval when no significant trend exists

# Revisions: rv10 removes confidence band loop in favor of <map_dfr>;
#  rv09 updates code; changes some defaults; updates name to <ci_band_lm> from <lm_band>;
#  rv08 adjusts for boundary case of perfect fit;
#  rv07 adds residuals checks for high leverage and standardized outliers;
#  rv06 alters method to construct simultaneous bands; also minor updates;
#  rv05 defaults to bootstrapped confidence interval on imputed NDs when trend
#  is non-significant; rv04 adds p-value and Rsq stats to output;
#  rv03 combines <lm_band> and <lm_impute_band> functions; adds significance test


ci_band_lm <- function(f,y=conc,nd=nd,clev=0.95,side=c('upr','lwr','both'),ns=50,...) {
  library(tidyverse)
  library(rlang)
  library(purrr)
  
## source('ci_compute.R')
## source('mc_vector.R')

# remove any rows with missing values in date or conc; set seed; check for type of confidence band
  side <- match.arg(side)
  
  n <- nrow(f)
  
  y <- ensym(y); nd <- ensym(nd)
  
  nd.pct <- with(f,inject(round(100*mean(!!nd,na.rm=T),1)))
  any.nd <- with(f,inject(any(!!nd == 1)))
  
  f <- f %>% filter(!(is.na(date) | is.na(!!y)))
  f <- f %>% mutate(conc.half=case_when(!!nd==0~!!y,!!nd==1~!!y/2))
  
##  set.seed(rseed,kind = "L'Ecuyer-CMRG")
  
# case of insufficient data
  if (n < 4) return(tibble(n,nd.pct,intercept=NA_real_,slope=NA_real_,pval=NA_real_,R2=NA_real_,xs=NA_real_,dd=NA,yhat=mean(f$conc.half),lwr=NA_real_,upr=NA_real_))
  
# convert dates to numeric if necessary
  f <- f %>% mutate(x=as.integer(date))

# compute cuts pts along range of x and predicted y values;
# also make sure confidence band is computed at each sample date so as to allow tests at those dates
  xs <- f$x %>% range %>% (function(rg) seq(rg[1],rg[2],length.out=101))
  xs <- unique(sort(c(xs,f$x)))

  
  
  assemble_fit0 <- function(xf,y,pval) {

    y <- ensym(y)
    lm.tmp <- inject(lm(!!y~x,data=xf))
    lm.tmp.sum <- summary(lm.tmp)

# if perfect linear fit, compute <fit0> specially
    if (is.nan(pval)) {
      fit0 <- tibble(intercept=lm.tmp.sum$coeff[1,1],slope=lm.tmp.sum$coeff[2,1],pval=NA_real_,R2=NA_real_)
    } else if (!is.nan(pval)) {
      fit0 <- tibble(intercept=lm.tmp.sum$coeff[1,1],slope=lm.tmp.sum$coeff[2,1],pval=lm.tmp.sum$coeff[2,4],R2=lm.tmp.sum$r.squared)
    }
    list(fit0=fit0,lmm=lm.tmp)
  }
  

  
  assemble_ci_fit <- function(xf,y,nd,xs) {
    
    y <- ensym(y); nd <- ensym(nd)
      
    if (side == 'upr') {
      tlev <- 2*clev - 1
      fit.ci <- inject(ci_compute(xf,!!y,!!nd,clev=tlev,side='both',meth='BOOT',impute = 'KM'))
      fit <- tibble(xs,dd=as.Date(xs,origin='1970-01-01'),yhat=fit.ci$ave,lwr=fit.ci$lcl,upr=fit.ci$ucl)
    } else if (side == 'lwr') {
      tlev <- 2*clev - 1
      fit.ci <- inject(ci_compute(xf,!!y,!!nd,clev=tlev,side='both',meth='BOOT',impute = 'KM'))
      fit <- tibble(xs,dd=as.Date(xs,origin='1970-01-01'),yhat=fit.ci$ave,lwr=fit.ci$lcl,upr=fit.ci$ucl)
    } else if (side == 'both') {
      fit.ci <- inject(ci_compute(xf,!!y,!!nd,clev=clev,side='both',meth='BOOT',impute = 'KM'))
      fit <- tibble(xs,dd=as.Date(xs,origin='1970-01-01'),yhat=fit.ci$ave,lwr=fit.ci$lcl,upr=fit.ci$ucl)
    }
    fit
  }

  
  assemble_cb_fit <- function(lmm,xs) {
    if (side %in% c('lwr','upr')) {
      Fval <- qf(clev,2,n-2)
      W <- sqrt(2*Fval)
      CI <- predict(lmm,tibble(x=xs),se.fit = T,interval = 'confidence',level = 2*clev-1)
      half.width <- W*CI$se.fit
      fit <- tibble(xs,dd=as.Date(xs,origin='1970-01-01'),yhat=CI$fit[,1],lwr=yhat - half.width,upr=yhat + half.width)

    } else if (side == 'both') {
      tlev <- (1 + clev)/2
      Fval <- qf(tlev,2,n-2)
      W <- sqrt(2*Fval)
      CI <- predict(lmm,tibble(x=xs),se.fit=T,interval = 'confidence',level=clev)
      half.width <- W*CI$se.fit
      fit <- tibble(xs,dd=as.Date(xs,origin='1970-01-01'),yhat=CI$fit[,1],lwr=yhat - half.width,upr=yhat + half.width)
    }
    fit
  }


  if (!any.nd) {
    lm.tst <- inject(lm(!!y~x,data = f))
  } else if (any.nd) {
    lm.tst <- lm(conc.half~x,data = f)
  }
  lm.tst.sum <- summary(lm.tst)
  pval.tst <- lm.tst.sum$coeff[2,4]
    
# check residuals for high leverage and extreme standardized outliers;
# but first check for no variation/perfect fit and hence unusable regression summary
  if (is.nan(pval.tst)) {
    hi.lev <- F; esi.out <- F
    f <- f %>% mutate(resid=0,esi=0)
  } else if (!is.nan(pval.tst)) {
    p <- 2
    lev.vec <- influence(lm.tst)$hat
    hi.lev <- lev.vec > 3*p/n
    se <- summary(lm.tst)$sigma
    f <- f %>% mutate(resid=lm.tst$residuals,esi=resid/(se*sqrt(1-lev.vec)))
    esi.out <- abs(f$esi) > 3
  }

# remove high leverage pts and/or extreme outliers and recompute linear fit
  if (any(hi.lev) || any(esi.out)) f <- f %>% filter(!(hi.lev | esi.out))
    
  if (!any.nd) {
    lm0 <- inject(lm(!!y~x,data=f))
  } else if (any.nd) {
    lm0 <- lm(conc.half~x,data=f)
  }
  lm0.sum <- summary(lm0)
  pval <- lm0.sum$coeff[2,4]

# construct confidence band components: <fit0> and <fit>
# case of all detects
  if (!any.nd) {
    out0 <- inject(assemble_fit0(f,!!y,pval))
    fit0 <- out0$fit0; lmm <- out0$lmm
    
# case of no trend; fit extruded confidence interval
    if (is.nan(pval) || pval >= 0.1) {
      fit <- inject(assemble_ci_fit(f,!!y,!!nd,xs))
# case of mild or better trend; fit simultaneous confidence band
    } else if (pval < 0.1) {
      fit <- assemble_cb_fit(lmm,xs)
    }
      
# case of any NDs
  } else if (any.nd) {
    tmp.out <- map(1:ns,.f=function(.x,xf,...) {
      xf <- xf %>% mutate(ix=.x,nd.r=0,conc.r=mc_vector(!!y,!!nd,...))
      out.r <- assemble_fit0(xf,conc.r,pval)
      fit0 <- out.r$fit0; lm.r <- out.r$lmm
        
# case of no trend; fit extruded confidence interval
      if (is.nan(pval) || pval >= 0.1) {
        fit <- assemble_ci_fit(xf,conc.r,nd.r,xs)
# case of mild or better trend; fit simultaneous confidence band
      } else if (pval < 0.1) {
        fit <- assemble_cb_fit(lm.r,xs)
      }
      list(fit0=fit0,fit=fit)
    },xf=f,...)

      
    fit0.tmp <- map_dfr(1:ns,.f=function(.x,xlst) {
      xlst[[.x]]$fit0
    },xlst=tmp.out)
      
    fit.tmp <- map_dfr(1:ns,.f=function(.x,xlst) {
      xlst[[.x]]$fit
    },xlst=tmp.out)
      
      
    fit0 <- fit0.tmp %>% summarise(across(c(intercept:R2),~ mean(.x,na.rm=T)))
      
    fit <- fit.tmp %>% group_by(xs,dd) %>% summarise(across(c(yhat:upr),~ mean(.x,na.rm=T)),.groups = 'drop')
  }
    


# output estimated trend line and confidence band fit
  tibble(side,n,nd.pct,clev,fit0,fit)
}

