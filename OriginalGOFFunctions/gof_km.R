# Author: Kirk Cameron, PhD, MacStat Consulting, Ltd
# Last update: March 20, 2023, rv05

# Purpose: use KM with ladder of powers and weighted correlation measure
# to pick best-fitting data transformation
# KM = Kaplan-Meier for left-censored data; 
# NOTE: this version only screens for normal-derived ladder of powers models 
# (excludes Weibull and gamma, though cube root is approximation for gamma)

# Revisions: rv05 updates code; name change from <gof_rob_km>;
#  rv04 updates code;
#  rv03 removes printing of probability plot; simply includes in list;
#  also checks for convergence of robust fit;
#  rv02 adds confidence band around probability plot; also plots all data on
#  probability plot, not just detects, using the KM imputations for NDs; also changes name
#  from <km_rob_gof> to <gof_rob_km>


gof_km <- function(f,y=conc,nd=nd,wt,hdr,pt.size=2) {
  
  library(tidyverse)
  library(robustbase)
#  library(DescTools)
  library(rlang)
  
  #  source('lop.R')
  #  source('tricube.R')
  #  source('impute_km.R')
  #  source('trans.R')
  
  y <- ensym(y); nd <- ensym(nd)
  
  if (missing(hdr)) hdr <- 'TEST'
  
  f <- f |> filter(!is.na(!!y))
  
  # ensure weights are set
  wt.exists <- hasName(f,'wt')
  
  if (missing(wt)) {
    if (!wt.exists) f <- f |> mutate(wt=1)
  } else {
    wt <- ensym(wt)
    f <- f |> mutate(wt=!!wt)
  }
  
  n <- nrow(f)
  
  # set conditions to define edge cases
  f <- f |> mutate(yv=!!y)
  
  any.nd <- with(f,inject(any(!!nd == 1)))
  all.nd <- with(f,inject(all(!!nd == 1)))
  no.nd <- with(f,inject(all(!!nd == 0)))
  nd.pct <- with(f,inject(round(100*mean(!!nd),1)))
  ndet <- with(f,inject(sum(1 - !!nd)))
  
  nd.hi <- F
  if (!(all.nd || no.nd)) nd.hi <- with(f,inject(min(yv[!!nd == 1]) >= max(yv[!!nd == 0])))
  
  edge.nd <- (all.nd || nd.hi || (ndet < 4 && any.nd) || nd.pct >= 70)
  
  
  # create censored probability plot for detects after checking ladder of powers,
  # but include nominal NDs on plot
  tf.list <- inject(lop(f,!!y,!!nd,wt=wt,hdr=hdr))
  
  tmax <- tf.list$corr.place$tf[1]
  corr.max <- tf.list$corr.place$tcorr[1]
  lab.max <- tf.list$corr.place$tf.lab[1]
  if(all(f$y==0)){corr.max=1;lab.max="Normal"} #quickfix to allow to not cause a log-transformation of zero when no goodness fits
  pow <- ifelse(tmax=='log',NA_real_,as.numeric(tmax))
  
  # transform detects and impute any non-detects according to best robust model
  # f.pp <- tf.list$f.pp
  impute.list <- inject(gof_impute_km(f,!!y,!!nd,wt=wt,tlab = lab.max))
  f.impute <- impute.list$f.pp
  rmod <- impute.list$rmod
  
  # create probability plot with robust fit; must fit model outside <ggplot2>
  ff <- f.impute
  ff <- ff |> mutate(nd2= case_when(nd == 0~'Detect',nd==1~'ND'))
  
  p <- ggplot(ff,aes(yt.hat,qd))
  p <- p + labs(x=paste0('Scaled Conc., Model = ',lab.max),y='Z-Score',title = hdr)
  ##  p <- p + geom_abline(slope=rmod$coef[2],intercept = rmod$coef[1],color='blue',linewidth=1.25) + theme_bw()
  p <- p + geom_line(aes(y=fit),color='blue',linewidth=1) + theme_bw()
  p <- p + geom_ribbon(aes(x=yt.hat,ymin=lwr,ymax=upr),color='lightblue',fill='lightblue',alpha=0.3)
  
  # check for all NDs or all detects; adjust symbol plotting accordingly
  if (all(ff$nd==1)) {
    p <- p + geom_point(data=ff,aes(color='ND'),size=pt.size,shape=0) + scale_color_manual('',breaks=c('ND'),values=c('red'))
  } else if (all(ff$nd==0)) {
    p <- p + geom_point(data=ff,aes(color='Detect'),size=pt.size,shape=16) + scale_color_manual('',breaks=c('Detect'),values=c('black'))
  } else {
    p <- p + geom_point(data=ff,aes(shape=nd2,color=nd2),size=pt.size) + scale_shape_manual('',values=c(Detect=16,ND=0)) + scale_color_manual('',values=c(Detect='black',ND='red'))
  }
  
  
  # output list of results
  rtab <- tibble(n=n,nd.pct,ndet=ndet,tf.max=lab.max,corr.max=round(corr.max,4),tf.2nd=tf.list$corr.place$tf.lab[2],corr.2nd=round(tf.list$corr.place$tcorr[2],4),tf.3rd=tf.list$corr.place$tf.lab[3],corr.3rd=round(tf.list$corr.place$tcorr[3],4))
  
  list(rtab=rtab,f.impute=f.impute,qplot=p,dplot=tf.list$dplot)
  
}
