# Author: Kirk Cameron, PhD, MacStat Consulting, Ltd
# Last Update: December 14, 2022, rv02

# Purpose: compute trendmap for user selected time period, either 'hist' for all
# historical data (default), or period = c(p.start,p.end), where p.start, p.end are date strings,
# but p.end is optional

# Revisions: rv02 changes name from <compute_trendmap>; updates code


trend_map_compute <- function(f,y=conc,nd=nd,east=east,north=north,locid=locid,period='hist',ns=100,...) {
  library(tidyverse)
  library(rlang)
  library(purrr)
  
  
##  source('trend_map.R')
  
  ys <- ensym(y); nds <- ensym(nd)
  easts <- ensym(east); norths <- ensym(north)
  locids <- ensym(locid)

# compute trendmap for selected time period
  trendmap.f <- f %>% group_by(!!locids) %>% do(trend_map(.,!!ys,!!nds,!!easts,!!norths,period,ns,...)) %>% ungroup()

# compute slope groupings
  if (all(trendmap.f$slope == 0)) {
      slope.vec <- c(-0.1,0.01,0.02,0.03,0.04,0.05)
  } else {
      slope.vec <- c(-0.1,quantile(abs(trendmap.f$slope[trendmap.f$slope != 0]),p=c(0.2,0.4,0.6,0.8,1),na.rm = T))
  }

# add slope tags to trendmap_df based on slope groupings
  trendmap.f <- trendmap.f %>% mutate(slope.tag=findInterval(abs(slope),slope.vec,rightmost.closed=T))
  
  trendmap.f
}
