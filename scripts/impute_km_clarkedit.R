# Author: Kirk Cameron, PhD, MacStat Consulting, Ltd
# Last revised: December 03, 2022, r06


# Purpose: impute NDs using Kaplan-Meier (KM) on left-censored data according to robust model;
# note: this version only considers normal-derived ladder of powers 
# models (need to add Weibull and gamma)

# Revisions: rv06 changes output to list; adds robust model;
#  rv05 increases tolerance in <lmrob> in cases of convergence error;
#  rv04 fixes case of too few detects;
#  rv03 changes name from <impute_rob_km> to <impute_km>;
#  rv02 changes name from <km_rob_impute> to <impute_km>


impute_km <- function(x,nd,id,tlab='Normal',wt=NULL) {
  library(tidyverse)
  library(robustbase)
  
 #source('scripts/pp_km.R')
 # source('scripts_original/trans.R')
 # source('scripts_original/backtrans.R')
 # source('scripts_original/pp_tricube.R')
  
  
  tf.df <- tibble(tf=c('log','0.125','0.1428571','0.1666666','0.2','0.25','0.3333333','0.5','1','2','3','4','5','6','7','8'),tf.lab=c('Log','Eighth Root','Seventh Root','Sixth Root','Fifth Root','Fourth Root','Cube Root','Square Root','Normal','Square','Cube','Fourth Power','Fifth Power','Sixth Power','Seventh Power','Eighth Power'))
  
# pull selected transformation
  tf <- tf.df %>% filter(tf.lab == tlab) %>% pull(tf)
  pow <- ifelse(tf=='log',NA_real_,as.numeric(tf))
  
# compute plotting positions using <km_pp>
  f.pp <- pp_km(x,nd,id,wt)
  f.pp <- f.pp %>% mutate(xhat=x,tf.lab=tlab)

  ndet <- sum(1-nd)
  

# if no NDs, all NDs, or < 4 detects, return data 'as is', along with data scaled according to selected transformation
  if (all(f.pp$nd == 0) || all(f.pp$nd == 1) || ndet < 4 ) {
    f.pp <- f.pp %>% mutate(xt=trans(x,tlab),xt.nom=xt)
    f.pp <- f.pp %>% arrange(pp,xhat)
    f.pp <- f.pp %>% mutate(qd=qnorm(pp))

    rmod <- try(with(f.pp,lmrob(qd~xt,method = 'SMDM',setting = 'KS2014')),silent=T)
    rfit <- try(predict(rmod,f.pp,interval = 'prediction',level = 0.99),silent = T)
    if ('try-error' %in% class(rfit) || any(is.na(rfit[,2])) || any(is.na(rfit[,3]))) {
      rmod <- try(with(f.pp,lmrob(qd~xt,method = 'SM',refine.tol=1e-5)),silent=T)
      rfit <- predict(rmod,f.pp,interval = 'prediction',level = 0.99)
    }
    if ('try-error' %in% class(rmod) || 'try-error' %in% class(rfit) || any(is.na(rfit[,2])) || any(is.na(rfit[,3]))) {
      rmod <- lm(qd~xt,data = f.pp)
      rfit <- predict(rmod,f.pp,interval = 'prediction',level = 0.99)
    }
    f.pp <- f.pp %>% mutate(fit=rfit[,1],lwr=rfit[,2],upr=rfit[,3],resid=qd-fit,adj.resid=case_when((qd >= lwr & qd <= upr)~0,TRUE~resid),rwt=pp_tricube(adj.resid,median((upr-lwr)/2)))
    
    list(f.pp=f.pp,rmod=rmod)
    
# case of mixture of detected, censored data
  } else {

# fit censored normal for transformed detects depending on given transformation
    f.pp <- f.pp %>% mutate(xt=trans(x,tlab),xt.nom=xt,qd=qnorm(pp))
    qd <- f.pp %>% filter(nd == 0) %>% pull(qd)
    xt <- f.pp %>% filter(nd == 0) %>% pull(xt)

# regress detects on z-scores; then compute imputed values for NDs
	  ld <- try(lmrob(xt~qd,method = 'SMDM',setting = 'KS2014'),silent = T)
	  if ('try-error' %in% class(ld) || is.na(ld$coef[1] || is.na(ld$coef[2]))) {
	    ld <- try(lmrob(xt~qd,method = 'SM',refine.tol=1e-5),silent=T)
	  }
	  if ('try-error' %in% class(ld) || is.na(ld$coef[1] || is.na(ld$coef[2]))) {
	    ld <- lm(xt~qd)
	  }
	  ahat <- ld$coef[1]
	  bhat <- ld$coef[2]
	
    f.tmp <- f.pp %>% filter(nd == 1) %>% mutate(qc=qnorm(pp),ndthat=ahat + bhat*qc)

# back-transform imputed non-detects
    f.tmp <- f.tmp %>% mutate(ndhat=backtrans(ndthat,tlab))
  
    ndhat <- f.tmp %>% pull(ndhat)
    ndthat <- f.tmp %>% pull(ndthat)

# construct dataframe of imputed data, unscaled and scaled, sorted by plotting position
	  f.pp$xhat[f.pp$nd == 1] <- ndhat
	  f.pp$xt[f.pp$nd == 1] <- ndthat
	  f.pp <- f.pp %>% arrange(pp,xhat)
	
# compute residuals from robust regression model; then weight residuals by
# degree of discrepancy from model using tricube weighting function
##	  f.pp <- f.pp %>% mutate(qd=qnorm(pp))
    f.tmp <- f.pp %>% filter(nd == 0)
    rmod <- try(with(f.tmp,lmrob(qd~xt,method = 'SMDM',setting = 'KS2014')),silent=T)
    rfit <- try(predict(rmod,f.pp,interval = 'prediction',level = 0.99),silent = T)
    if ('try-error' %in% class(rmod) ||'try-error' %in% class(rfit) || any(is.na(rfit[,2])) || any(is.na(rfit[,3]))) {
      rmod <- try(lmrob(qd~xt,data=f.tmp,method = 'SM',refine.tol=1e-5,k.max=300),silent=T)
      rfit <- try(predict(rmod,f.pp,interval = 'prediction',level = 0.99),silent = T)
    }
    if ('try-error' %in% class(rmod) || 'try-error' %in% class(rfit) || any(is.na(rfit[,2])) || any(is.na(rfit[,3]))) {
      rmod <- lm(qd~xt,data = f.tmp)
      rfit <- predict(rmod,f.pp,interval = 'prediction',level = 0.99)
    }
    f.pp <- f.pp %>% mutate(fit=rfit[,1],lwr=rfit[,2],upr=rfit[,3],resid=qd-fit,adj.resid=case_when((qd >= lwr & qd <= upr)~0,TRUE~resid),rwt=pp_tricube(adj.resid,median((upr-lwr)/2)))
	  
	  list(f.pp=f.pp,rmod=rmod)
  }
}
