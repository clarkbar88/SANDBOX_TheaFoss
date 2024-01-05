# Author: Kirk Cameron, PhD, MacStat Consulting, Ltd
# Last update: April 22, 2023, r09


# Purpose: impute NDs using Kaplan-Meier (KM) on left-censored data according to robust model;
# note: this version only considers normal-derived ladder of powers 
# models (need to add Weibull and gamma)

# Revisions: rv09 updates outputs; adds separate branch for all detects;
#  rv08 updates code; expands edge cases; for edge cases, passes through 
#  Monte Carlo imputation from <pp_km>;
#  rv07 adds more error checks; returns NULL for robust model if no convergent fit;
#. rv06 changes output to list; adds robust model;
#  rv05 increases tolerance in <lmrob> in cases of convergence error;
#  rv04 fixes case of too few detects;
#  rv03 changes name from <impute_rob_km> to <impute_km>;
#  rv02 changes name from <km_rob_impute> to <impute_km>


gof_impute_km <- function(f,y=conc,nd=nd,wt,tlab='Normal',...) {
  
  library(tidyverse)
  library(robustbase)
  library(rlang)
  
##  source('pp_km.R')
##  source('trans.R')
##  source('backtrans.R')
  
  y <- ensym(y); nd <- ensym(nd)
  
# ensure weights are set
  wt.exists <- hasName(f,'wt')
  
  if (missing(wt)) {
    if (!wt.exists) f <- f |> mutate(wt=1)
  } else {
    wt <- ensym(wt)
    f <- f |> mutate(wt=!!wt)
  }

# compute KM plotting positions
  f.pp <- inject(pp_km(f,!!y,!!nd,wt=wt,...))
  

# define set of possible transformations and their labels
  tf.f <- tibble::tribble(
    ~tf, ~tf.lab,
    'log', 'Log',
    '0.125', 'Eighth Root',
    '0.1428571', 'Seventh Root',
    '0.1666666', 'Sixth Root',
    '0.2', 'Fifth Root',
    '0.25', 'Fourth Root',
    '0.3333333', 'Cube Root',
    '0.5', 'Square Root',
    '1', 'Normal',
    '2', 'Square',
    '3', 'Cube',
    '4', 'Fourth Power',
    '5', 'Fifth Power',
    '6', 'Sixth Power',
    '7', 'Seventh Power',
    '8', 'Eighth Power'
  )
  
  
# pull selected transformation
  tf <- tf.f |> filter(tf.lab == tlab) |> pull(tf)
  pow <- ifelse(tf=='log',NA_real_,as.numeric(tf))
  

# define edge cases
  f.pp <- f.pp |> mutate(yv=!!y)
  
  any.nd <- with(f.pp,inject(any(!!nd == 1)))
  all.nd <- with(f.pp,inject(all(!!nd == 1)))
  no.nd <- with(f.pp,inject(all(!!nd == 0)))
  nd.pct <- with(f.pp,inject(round(100*mean(!!nd),1)))
  ndet <- with(f.pp,inject(sum(1 - !!nd)))
  
  nd.hi <- F
  if (!(all.nd || no.nd)) nd.hi <- with(f.pp,inject(min(yv[!!nd == 1]) >= max(yv[!!nd == 0])))
  
  edge.nd <- (all.nd || nd.hi || (ndet < 4 && any.nd) || nd.pct >= 70)
  
# handle case of all detects
  if (no.nd) {
    
    f.pp <- f.pp |> mutate(tf.lab=tlab,yhat=!!y,yt=trans(!!y,tlab),yt.hat=trans(yhat,tlab),qd=qnorm(pp)) |> dplyr::select(-yv)
    
    rmod <- try(lmrob(qd~yt.hat,data = f.pp,weights = f.pp$pp.wt,method = 'SMDM',setting = 'KS2014'),silent = T)
    rfit <- try(predict(rmod,f.pp,interval = 'confidence',level = 0.99,weights = f.pp$pp.wt),silent = T)
    if ('try-error' %in% class(rmod) || 'try-error' %in% class(rfit) || any(is.na(rfit[,2])) || any(is.na(rfit[,3]))) {
      rmod <- try(lmrob(qd~yt.hat,data = f.pp,weights = f.pp$pp.wt,method = 'SM',refine.tol=1e-5,k.max=300),silent = T)
      rfit <- try(predict(rmod,f.pp,interval = 'confidence',level = 0.99,weights=f.pp$pp.wt),silent = T)
    }
    if ('try-error' %in% class(rmod) || 'try-error' %in% class(rfit) || any(is.na(rfit[,2])) || any(is.na(rfit[,3]))) {
      rmod <- lm(qd~yt.hat,data = f.pp,weights = pp.wt)
      rfit <- predict(rmod,f.pp,interval = 'confidence',level = 0.99,weights=f.pp$pp.wt)
    }
    
    f.pp <- f.pp |> mutate(fit=rfit[,1],lwr=rfit[,2],upr=rfit[,3])


# handle edge cases; pass through Monte Carlo imputations from <pp_km>
  } else if (edge.nd) {
    
    f.pp <- f.pp |> mutate(tf.lab=tlab,yhat=y.r,yt=trans(!!y,tlab),yt.hat=trans(yhat,tlab),qd=qnorm(pp)) |> dplyr::select(-y.r)

    rmod <- try(lmrob(qd~yt.hat,data = f.pp,weights = f.pp$pp.wt,method = 'SMDM',setting = 'KS2014'),silent = T)
    rfit <- try(predict(rmod,f.pp,interval = 'confidence',level = 0.99,weights = f.pp$pp.wt),silent = T)
    if ('try-error' %in% class(rmod) || 'try-error' %in% class(rfit) || any(is.na(rfit[,2])) || any(is.na(rfit[,3]))) {
      rmod <- try(lmrob(qd~yt.hat,data = f.pp,weights = f.pp$pp.wt,method = 'SM',refine.tol=1e-5,k.max=300),silent = T)
      rfit <- try(predict(rmod,f.pp,interval = 'confidence',level = 0.99,weights=f.pp$pp.wt),silent = T)
    }
    if ('try-error' %in% class(rmod) || 'try-error' %in% class(rfit) || any(is.na(rfit[,2])) || any(is.na(rfit[,3]))) {
      rmod <- lm(qd~yt.hat,data = f.pp,weights = pp.wt)
      rfit <- predict(rmod,f.pp,interval = 'confidence',level = 0.99,weights=f.pp$pp.wt)
    }

    f.pp <- f.pp |> mutate(fit=rfit[,1],lwr=rfit[,2],upr=rfit[,3])


# handle non-edge cases; fit censored normal for transformed detects depending on given transformation
  } else if (!edge.nd) {

    f.pp <- f.pp |> mutate(tf.lab=tlab,yhat=yv,yt=trans(yv,tlab),yt.hat=yt,qd=qnorm(pp)) |> dplyr::select(-yv)
    x.pp <- f.pp |> filter(!!nd == 0)

# regress detects on z-scores; then compute imputed values for NDs
	  ld <- try(lmrob(yt.hat~qd,data = x.pp,weights = x.pp$pp.wt,method = 'SMDM',setting = 'KS2014'),silent = T)
	  if ('try-error' %in% class(ld) || is.na(ld$coef[1] || is.na(ld$coef[2]))) {
	    ld <- try(lmrob(yt.hat~qd,data=x.pp,weights = x.pp$pp.wt,method = 'SM',refine.tol=1e-5,k.max=300),silent = T)
	  }
	  if ('try-error' %in% class(ld) || is.na(ld$coef[1] || is.na(ld$coef[2]))) {
	    ld <- lm(yt.hat~qd,data = x.pp,weights = pp.wt)
	  }
	  
	  ahat <- ld$coef[1]
	  bhat <- ld$coef[2]
	
# back-transform imputed non-detects
    f.tmp1 <- f.pp |> filter(!!nd == 1) |> mutate(yt.hat=ahat + bhat*qd,yhat=backtrans(yt.hat,tlab))

# construct dataframe of imputed data, unscaled and scaled, sorted by plotting position

    f.tmp2 <- f.pp |> filter(!!nd == 0)
    f.pp <- bind_rows(f.tmp1,f.tmp2) |> arrange(pp,yhat,desc(!!nd))

	
# fit robust regression model of z-scores on detects (flipped from above)
    x.pp <- f.pp |> filter(!!nd == 0)
    rmod <- try(lmrob(qd~yt.hat,data = x.pp,weights = x.pp$pp.wt,method = 'SMDM',setting = 'KS2014'),silent = T)
    rfit <- try(predict(rmod,f.pp,interval = 'confidence',level = 0.99,weights=f.pp$pp.wt),silent = T)
    if ('try-error' %in% class(rmod) || 'try-error' %in% class(rfit) || any(is.na(rfit[,2])) || any(is.na(rfit[,3]))) {
      rmod <- try(lmrob(qd~yt.hat,data = x.pp,weights = x.pp$pp.wt,method = 'SM',refine.tol=1e-5,k.max=300),silent = T)
      rfit <- try(predict(rmod,f.pp,interval = 'confidence',level = 0.99,weights=f.pp$pp.wt),silent = T)
    }
    if ('try-error' %in% class(rmod) || 'try-error' %in% class(rfit) || any(is.na(rfit[,2])) || any(is.na(rfit[,3]))) {
      rmod <- lm(qd~yt.hat,data = x.pp,weights = pp.wt)
      rfit <- predict(rmod,f.pp,interval = 'confidence',level = 0.99,weights=f.pp$pp.wt)
    }

    f.pp <- f.pp |> mutate(fit=rfit[,1],lwr=rfit[,2],upr=rfit[,3])
  }
	  
  list(f.pp=f.pp,rmod=rmod)
}
