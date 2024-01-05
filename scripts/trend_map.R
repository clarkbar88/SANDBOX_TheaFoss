# Author: Kirk Cameron, PhD, MacStat Consulting, Ltd
# Last Update: December 14, 2022, rv13

# Purpose: Compute trendmap values at one location; 
# <period> governs time window used: hist(orical) for all data; recent for
# latest <span> yrs (>= 5 events minimum); new for latest <span> sampling events (minimum 5+ events)

# Revisions: rv13 updates code to match subroutines; alters name to <trend_map> from <trendmap>;
#  rv12 dumps ATS in favor of LM trend modeling; adds flexibility to period of interest;
#  rv11 updates code;
#  rv10 converts slope estimate to annual rate instead of daily rate of change;
#  rv09 adds external sources; allows variable span for new trends; incorporates
#  updated datafile structure;
#  rv08 computes slope for type FLAT since this only indicates tie.pct=100;
#  rv07 aligns GTS trend type classifications;
#  rv06 renames dependencies and function name;
#  rv05 renames function from <trendmap.loc>;
#  rv04 uses Akritas-Theil-Sen (ATS) to estimate linear slope of
#  censored data; other minor tweaks, including option for how long recent period should be;
#  rv03 changes options to ensure output dataframe keeps strings as character;
#  rv02 changes variable names to match current GTS implementation


trend_map <- function(f,y=conc,nd=nd,east=east,north=north,period='hist',ns=100,...) {
  library(tidyverse)
  library(rlang)
  library(purrr)
  library(broom)

##  source('trend_lm.R')
##  source('pus.R')
  
  ys <- ensym(y); nds <- ensym(nd)
  easts <- ensym(east); norths <- ensym(north)
  east <- f %>% slice(1L) %>% pull(!!easts)
  north <- f %>% slice(1L) %>% pull(!!norths)
  
  # remove any rows with missing values in date or conc; set seed
  f <- f %>% filter(!(is.na(date) | is.na(!!ys)))

  n <- nrow(f)
  d <- f %>% pus(date)
  ndate <- length(d)
  
  if (period != 'hist') {
    nper <- length(period)
    if (nper == 1) {
      start <- as.Date(period,origin='1970-01-01'); end <- d[ndate]
      ndate <- sum(d >= start)
      if (ndate == 0) end <- start
      f.tmp <- f %>% filter(date >= start)
    } else if (nper == 2) {
      start <- as.Date(period[1],origin='1970-01-01')
      end <- as.Date(period[2],origin='1970-01-01')
      ndate <- sum(d >= start & d <= end)
      f.tmp <- f %>% filter(date >= start,date <= end)
    }
  } else if (period == 'hist') {
    start <- d[1]; end <- d[ndate]
    f.tmp <- f %>% filter(date >= start)
  }
  
  n <- nrow(f.tmp)
  if (n == 0) {
    nd.pct <- 0; all.nd <- F
  } else if (n > 0) {
    nd.pct <- with(f.tmp,round(100*inject(mean(!!nds)),1))
    all.nd <- with(f.tmp,inject(all(!!nds == 1)))
  }
    
# compute approximate linear slope and statistical significance, with Monte Carlo imputation of any NDs

  if (ndate < 4 || all.nd) {
    slope <- 0; pval <- 1
    trend_f <- tibble(east,north,n,nd.pct,start,end,intercept=NA_real_,slope,units=paste0(f$units[1],'/yr'),pval,trend.tag=case_when(all.nd~'ND',TRUE~'N/A'))
  } else {
    tmp <- inject(trend_lm(f,!!ys,!!nds,...))
# convert slope from daily rate to annual rate of change
    tmp <- tmp %>% mutate(slope=slope*365,units=paste0(f$units[1],'/yr'),start=start,end=end,east=east,north=north)
  
# use slope significance to classify each trend
    tmp <- tmp %>% mutate(trend.tag=case_when(pval >= 0.1~'FLAT',pval < 0.05 && slope > 0~'INCR',pval >= 0.05 && pval < 0.1 && slope > 0~'PROB-INCR',pval < 0.05 && slope < 0~'DECR',pval >= 0.05 && pval < 0.1 && slope < 0~'PROB-DECR'))
    trend_f <- tmp %>% select(c(east,north,n,nd.pct,start,end,intercept,slope,units,pval,trend.tag))
  }

  trend_f <- trend_f %>% mutate(trend.tag=factor(trend_f$trend.tag,levels=c('DECR','PROB-DECR','ND','FLAT','PROB-INCR','INCR','N/A'),ordered=T))
  
  trend_f
}
