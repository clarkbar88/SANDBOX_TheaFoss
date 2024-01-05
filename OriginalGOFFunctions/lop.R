# Author: Kirk Cameron, PhD, MacStat Consulting, Ltd
# Last revised: March 20, 2023, rv11

# Purpose: Compute best augmented ladder of powers transformation via a robust correlation
# measure, given a data set and set of plotting positions for use with censored data; NOTE:
# this version only considers normal-derived ladder of powers models (excludes gamma and Weibull)

# Revisions: rv11 computes correlation on each of series of Monte Carlo imputation instances;
#  rv10 updates code; name change from <lop_rob>; replaces <robcor> with weighted correlation;
#  rv09 minor update;
#  rv08 eliminates printing of dotplot; uses Monte
#  Carlo imputation for edge cases;
#  rv07 changes name from <lop_rob_cens> to <lop_rob>;
#  rv06 switches to KM imputed NDs;
#  rv05 adds option of pre-weighted data;
#  rv04 improves handling of negative values;
#  rv03 adds index variable to track specific records;
#  rv02 name change


lop <- function(f,y=conc,nd=nd,wt,hdr,ns=20,...) {
  
  library(tidyverse)
  library(robcor)
  library(rlang)
  
## source('pp_km.R')
## source('qqcor_dotchart.R')
## source('mc_vector.R')
  
  y <- ensym(y); nd <- ensym(nd)
  
  if (missing(hdr)) hdr <- 'TEST'

## define possible transformations and labels
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
  
## remove missing data
  f <- f |> filter(!is.na(!!y)) |> mutate(yv=!!y)
  
## check weights
  wt.exists <- hasName(f,'wt')
  
  if (missing(wt)) {
    if (!wt.exists) f <- f |> mutate(wt=1)
  } else {
    wt <- ensym(wt)
    f <- f |> mutate(wt=!!wt)
  }
  
## define edge cases
  n <- nrow(f)
  ndet <- with(f,inject(sum(1 - !!nd)))
  nd.pct <- with(f,inject(round(100*mean(!!nd),1)))
  any.nd <- with(f,inject(any(!!nd == 1)))
  all.nd <- with(f,inject(all(!!nd == 1)))
  no.nd <- with(f,inject(all(!!nd == 0)))

  nd.hi <- F
  if (!(all.nd || no.nd)) nd.hi <- with(f,inject(min(yv[!!nd == 1]) >= max(yv[!!nd == 0])))

  edge.nd <- (all.nd || nd.hi || (ndet < 4 && any.nd) || nd.pct >= 70)

  
## helper function to compute vector of robust correlations
  assemble_cor <- function(fpp,edge.nd=F,...) {
    tlab <- tf.f$tf.lab
    td <- fpp |> mutate(qd=qnorm(pp))
    
    if (edge.nd) {
      cor.v <- map_dbl(tlab,.f=function(.x,xf) {
        xf |> mutate(yt=trans(y.r,.x),bad=is.na(yt) | is.nan(yt) | is.infinite(yt)) |> summarise(corr=case_when(any(bad)~NA_real_,sd(y.r)==0~0,.default = cov.wt(tibble(yt[!bad],qd[!bad]),wt=pp.wt[!bad],cor=T)$cor[2,1])) |> pull(corr)
      },xf=td)
    } else if (!edge.nd) {
      cor.v <- map_dbl(tlab,.f=function(.x,xf) {
        xf |> mutate(yt=trans(!!y,.x),bad=is.na(yt) | is.nan(yt) | is.infinite(yt)) |> summarise(corr=case_when(any(bad)~NA_real_,sd(!!y)==0~0,.default = cov.wt(tibble(yt[!bad],qd[!bad]),wt=pp.wt[!bad],cor=T)$cor[2,1])) |> pull(corr)
      },xf=td)
    }
    cor.v
  }

  
## handle ND edge cases by using all of <f.pp> to compute correlations;
## determine maximum, 2nd, and 3rd largest correlations between detects and each model
  if (edge.nd) {
    tmp.lst <- tibble(ix=1:ns) |> rowwise() |> mutate(tab=list({
      f.pp <- inject(pp_km(f,!!y,!!nd,wt=wt))
      out <- tf.f |> mutate(tcorr=assemble_cor(f.pp,edge.nd=T),io=1:nrow(tf.f))
      list(f.pp=f.pp,out=out)
    })) |> ungroup() |> hoist(tab,'f.pp','out')
    
    cor.f <- tmp.lst |> select(out) |> unnest(cols = c(out)) |> group_by(io) |> summarise(tf.lab=tf.lab[1],tf=tf[1],tcorr=mean(tcorr,na.rm = T),.groups = 'drop') |> select(-io)
    f.pp <- tmp.lst |> select(f.pp) |> unnest(cols = c(f.pp)) |> group_by(ir) |> summarise(across(c(!!y,!!nd,y.r,pp.wt,pp),~mean(.x,na.rm=T)),.groups = 'drop') |> select(-ir) |> arrange(pp,desc(!!nd))
    
    corr.place <- cor.f |> arrange(desc(tcorr)) |> slice_head(n=3)
    
## non-edge cases
  } else if (!edge.nd) {
    f.pp <- inject(pp_km(f,!!y,!!nd,wt=wt))
    x.pp <- f.pp |> filter(nd == 0)
    cor.f <- tf.f |> mutate(tcorr=assemble_cor(x.pp,edge.nd=F))
    corr.place <- cor.f |> arrange(desc(tcorr)) |> slice_head(n=3)
  }
  
  
## construct dotplot of robust correlation results across transformations
  dplot <- with(cor.f,qqcor_dotchart(tcorr,tf.lab,hdr))

  
  return(list(corr.place=corr.place,f.pp=f.pp,dplot=dplot))
}
