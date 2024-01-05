# Author: Kirk Cameron, PhD, MacStat Consulting, Ltd
# Last revised: April 13, 2023, rv06

# Purpose: Use KM to estimate plotting positions (aka percentage points) for
# left-censored data; KM = Kaplan-Meier;
# defaults to Weibull plotting positions if no censoring is present


# Revisions: rv06 outputs both <wt> and <pp.wt> variables;
#  rv05 constructs single Monte Carlo instance if edge case;
#  rv04 updates code;
#  rv03 switches x to y for consistency;
#  rv02 changes name from <km_pp> to <pp_km>


pp_km <- function(f,y=conc,nd=nd,wt,eps=.Machine$double.eps,...) {
  
  library(tidyverse)
  library(rlang)
  
  y <- ensym(y); nd <- ensym(nd)
  

# remove any missing data
  f <- f |> filter(!is.na(!!y)) |> mutate(yv=!!y)
  
  
# normalize weights
  wt.exists <- hasName(f,'wt')
  
  if (missing(wt)) {
    if (wt.exists) {
      f <- f |> mutate(N=n(),pp.wt=N*wt/sum(wt))
    } else {
      f <- f |> mutate(N=n(),wt=1,pp.wt=1)
    }
  } else {
    wt <- ensym(wt)
    f <- f |> mutate(N=n(),pp.wt=N*!!wt/sum(!!wt))
  }
  
# set conditions for defining edge cases
  any.nd <- with(f,inject(any(!!nd == 1)))
  all.nd <- with(f,inject(all(!!nd == 1)))
  no.nd <- with(f,inject(all(!!nd == 0)))
  nd.pct <- with(f,inject(round(100*mean(!!nd),1)))
  ndet <- with(f,inject(sum(1 - !!nd)))
  
  nd.hi <- F
  if (!(all.nd || no.nd)) nd.hi <- with(f,inject(min(yv[!!nd == 1]) >= max(yv[!!nd == 0])))
  
  nr <- nrow(f)
  
# default/edge cases of no censoring, complete censoring, <= 3 detects, all NDs > all detects, or ND pct >= 70%;
# use Monte Carlo imputation if any NDs to compute plotting positions

  edge.nd <- (all.nd || nd.hi || (ndet < 4 && any.nd) || nd.pct >= 70)

# case of all detects (no censoring)
  if (no.nd) {
    f.pp <- f |> arrange(yv,desc(!!nd)) |> select(c(!!y,!!nd,wt,pp.wt)) |> mutate(pp=cumsum(pp.wt)/(nr+1))

# edge cases with NDs; construct single instance of Monte Carlo imputation
  } else if (edge.nd) {
    f.pp <- f |> mutate(y.r=mc_vector(yv,!!nd,...),ir=1:nr) |> arrange(y.r,desc(!!nd)) |> select(c(ir,!!y,!!nd,y.r,wt,pp.wt)) |> mutate(pp=cumsum(pp.wt)/(nr+1))

# non-default/non-edge cases
  } else {

# table by nd status; determine unique x values and number at risk; build KM CDF
    km.tab <- f |> mutate(yy=yv-!!nd*eps) |> group_by(yy,!!nd) |> summarise(km.lev=yv[1],nn=n(),sw=sum(pp.wt),.groups = 'drop') |> mutate(km.rsk=cumsum(sw),tab1=case_when(!!nd==0~sw,!!nd==1~0),km.cdf=rev(c(1,cumprod(1 - rev(tab1[-1])/rev(km.rsk[-1])))),km.surv=1 - km.cdf)
    km.tab <- km.tab |> select(-yy) |> rename(sum.wt=sw) |> relocate(km.lev,.before = !!nd)


# plotting positions for detects
    tmp <- km.tab |> filter(!!nd == 0) |> mutate(km.cdf=km.cdf*(nr/(nr+1))) |> select(c(km.lev,!!nd,km.cdf)) |> rename(yv=km.lev,pp=km.cdf)
    f.d <- f |> filter(!!nd == 0) |> left_join(tmp,join_by(yv,!!nd)) |> select(c(!!y,!!nd,wt,pp.wt,pp))
  
# plotting positions for NDs
    tmp <- km.tab |> filter(!!nd == 1) |> mutate(km.cdf=km.cdf*(nr/(nr+1)),km.inc=km.cdf/(nn+1)) |> select(c(km.lev,!!nd,km.inc)) |> rename(yv=km.lev)
    f.nd <- f |> filter(!!nd == 1) |> arrange(yv) |> left_join(tmp,join_by(yv,!!nd)) |> group_by(yv) |> mutate(pp=cumsum(km.inc)) |> ungroup() |> select(-km.inc) |> select(c(!!y,!!nd,wt,pp.wt,pp))


# return tibble with concatenated data and plotting positions
    f.pp <- bind_rows(f.d,f.nd) |> arrange(pp,!!y,desc(!!nd))
  }
  
  f.pp
}
