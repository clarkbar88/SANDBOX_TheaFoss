# Tacoma Stormwater Analysis, December 2022


# SETUP ----

library(tidyverse)
library(readxl)
library(lubridate)
library(R.utils)
library(coin)
library(rlang)
library(purrr)
library(broom)
library(stats)
library(lmPerm)
library(boot)
library(magrittr)
library(RColorBrewer)


source("scripts/backtrans.R")
source("scripts/ci_band_lm_Nov15.R")
source("scripts/ci_boot_Nov15.R")
source("scripts/ci_band_plot_Nov15.R")
source("scripts/ci_compute.R")
source("scripts/convex_hull.R")
source("scripts/eplot.R")
source("scripts/impute_km_clarkedit.R")
source("scripts/mc_vector.R")
source("scripts/Mode.R")
source("scripts/oneway_mc.R")
source("scripts/oneway_mc_plot.R")
source("scripts/pp_km.R")
source("scripts/pp_tricube.R")
source("scripts/pus.R")
source("scripts/trans.R")
source("scripts/trend_lm.R")
source("scripts/trend_map.R")
source("scripts/trend_map_compute.R")
source("scripts/trend_map_plot.R")
source("scripts/whdquantile.R")
source("scripts/wquantile_generic.R")


# source('add_gis.R')
# source('ci_compute.R')
# source('pctdiff_boot.R')


site <- 'Tacoma'; site.tag <- 'Tacoma'; site.filetag <- 'Tacoma_Stormwater'
user.cut <- '2022-01-01'
new.yr <- '2022'
yr.cut <- '2022-01-01'

fdate <- format.Date(Sys.Date(),format='%y%m%d')

# WRANGLE ----

a0 <- read_excel('data_raw/WY2021_lab_data_kc.xlsx',col_types = c(rep('text',3),rep('numeric',2),rep('text',3),rep('numeric',3),'date',rep('text',9),rep('numeric',9),rep('text',2)))
names(a0) <- tolower(names(a0))

a00 <- a0 %>% select(-c(sys_loc_code,loc_number,`_237a_comp`,task_code,sys_sample_code,cipp_lining,halfnd,loghalfnd,totphthsum,number,forstats,forstats_5yr,forstats_2yr,validator_qualifiers))

a00 <- a00 %>% mutate(id=1:n(),wt=1)
a00 <- a00 %>% mutate(nd=1-detect) %>% rename(conc=result_numeric,units=result_unit,qual=lab_qualifiers,cas=cas_rn) %>% select(-c(detect,result_text,fraction_label))

a00 <- a00 %>% mutate(date=as.Date(logdate)) %>% select(-c(logdate)) %>% rename(matrix=matrix_code)

a01 <- a00 %>% mutate(locid=case_when(locid == '237A New'~'237A',TRUE~locid))

a01 <- a01 %>% mutate(units=case_when(units == '#/100mL'~'#/100ml',units %in% c('CFU/100 ml','CFU/100 mL','CFU/100mL')~'cfu/100ml',units=='mg/l'~'mg/L',units=='ug/l'~'ug/L',units=='pH Units'~'SU',TRUE~units))

a01 <- a01 %>% group_by(coc) %>% mutate(units=Mode(units)) %>% ungroup()

tst.coc <- a01 %>% group_by(coc) %>% summarise(min.d=min(date),max.d=max(date),.groups = 'drop')
coc.elim <- tst.coc %>% filter(max.d < ymd('2021-01-01') | min.d > ymd('2019-01-01')) %>% pull(coc)

a02 <- a01 %>% filter(!coc %in% coc.elim)
a02 <- a02 %>% mutate(locid=factor(locid,ordered = T))

f <- a02

# EXPLORE ----

## Box Plots ----
fn <- paste0(site.filetag,'_boxplots_ByOutfall_',fdate,'.pdf')
pdf(file=fn,w=11,h=8.5)
box.outfall <- f %>% group_by(coc) %>% do(plot={
  tcoc <- .$coc[1]
  print(tcoc)
  hdr <- paste0('Box Plots for ',tcoc,' by Outfall')
  print(eplot(.,thead=hdr,'box',vx=locid,vfill=locid,rotate.lab = T,pt.size = 1.5,varwidth=T))
})
dev.off()

fn <- paste0(site.filetag,'_boxplots5y_ByOutfall_',fdate,'.pdf')
pdf(file=fn,w=11,h=8.5)
box.outfall.5y <- f %>% group_by(coc) %>% filter(date >= ymd('2017-01-01')) %>% do(plot={
  tcoc <- .$coc[1]
  print(tcoc)
  hdr <- paste0('Box Plots for ',tcoc,' by Outfall Since 2017')
  print(eplot(.,thead=hdr,'box',vx=locid,vfill=locid,rotate.lab = T,pt.size = 1.5,varwidth=T))
})
dev.off()


## Time Series Plots ----
fn <- paste0(site.filetag,'_tsplots_ByOutfall_',fdate,'.pdf')
pdf(file=fn,w=11,h=8.5)
tsplots.outfall <- f %>% group_by(coc) %>% do(plot={
  tcoc <- .$coc[1]
  print(tcoc)
  hdr <- paste0('Time Series Plots for ',tcoc,' by Outfall')
  print(eplot(.,thead=hdr,'ts',rotate.lab = T,pt.size = 1.5))
}) %>% ungroup()
dev.off()

fn <- paste0(site.filetag,'_tsplots5y_ByOutfall_',fdate,'.pdf')
pdf(file=fn,w=11,h=8.5)
tsplots.outfall.5y <- f %>% group_by(coc) %>% filter(date >= ymd('2017-01-01')) %>% do(plot={
  tcoc <- .$coc[1]
  print(tcoc)
  hdr <- paste0('Time Series Plots for ',tcoc,' by Outfall Since 2017')
  print(eplot(.,thead=hdr,'ts',rotate.lab = T,pt.size = 1.5))
}) %>% ungroup()
dev.off()


## Time Series Plots with LOESS ----
fn <- paste0(site.filetag,'_tsplots_loess_ByOutfall_',fdate,'.pdf')
pdf(file=fn,w=11,h=8.5)
tsplots.outfall.loess <- f %>% group_by(coc) %>% do(plot={
  tcoc <- .$coc[1]
  print(tcoc)
  hdr <- paste0('Time Series Plots for ',tcoc,' By Outfall With LOESS Confidence Band')
  print(eplot(.,thead=hdr,'ts',showlines = F,rotate.lab = T,loess = T,pt.size = 1.5))
}) %>% ungroup()
dev.off()


fn <- paste0(site.filetag,'_tsplots5y_loess_ByOutfall_',fdate,'.pdf')
pdf(file=fn,w=11,h=8.5)
tsplots.outfall.loess.5y <- f %>% group_by(coc) %>% filter(date >= ymd('2017-01-01')) %>% do(plot={
  tcoc <- .$coc[1]
  print(tcoc)
  hdr <- paste0('Time Series Plots for ',tcoc,' By Outfall Since 2017 With LOESS Confidence Band')
  print(eplot(.,thead=hdr,'ts',showlines = F,rotate.lab = T,loess = T,pt.size = 1.5))
}) %>% ungroup()
dev.off()


# TREND TESTING/MAPPING ----

## import coordinates
outfall.loc0 <- read_excel("data_raw/EIMLocation2021_kc.xlsx")
outfall.loc <- outfall.loc0 %>% select(c(loc_id,loc_name,lon,lat)) %>% mutate(outfall=c('OF237B','OF235','OF245','OF243','OF254','OF230','OF237A','RG15CUW','FD2_237A','FD6_235','FD3NEW_230'))

## approximate location of outfall OF248
outfall.loc <- bind_rows(outfall.loc,tibble(loc_id='UNKNOWN',loc_name='APPROX_OF248',lon=-122.43,lat=47.248,outfall='OF248'))
outfall.loc <- outfall.loc %>% mutate(locid=sub('OF','',outfall))

## construct trend maps by outfall and sediment trap
tmp <- outfall.loc %>% select(locid,lon,lat)
f <- f %>% left_join(.,tmp)
f <- f %>% mutate(locid=factor(locid,ordered = T))

## create outfall trend maps
outfall.trend.map.hist <- f %>% group_by(coc) %>% do({
  td <- .
  tcoc <- td$coc[1]
  print(tcoc)
  trend_map_compute(td,east = lon,north = lat)
}) %>% ungroup()

outfall.trend.map.5y <- f %>% group_by(coc) %>% do({
  td <- .
  tcoc <- td$coc[1]
  print(tcoc)
  trend_map_compute(td,east = lon,north = lat,period = '2017-01-01')
}) %>% ungroup()


fn <- paste0(site.filetag,'_TrendMapsHist_ByOutfall_',fdate,'.pdf')
pdf(file=fn,w=11,h=8.5)
outfall.trendmap.plots.hist <- outfall.trend.map.hist %>% group_by(coc) %>% do(plot={
  tmap <- .
  tcoc <- tmap$coc[1]
  print(tcoc)
  hdr <- paste0('Historical Trend Map by Outfall for ',tcoc)
  trend_map_plot(tmap,period='hist',hdr=hdr,gis=F,gis.lst=NULL,bex=0.01,interact=F,wlabs=T)
}) %>% ungroup()
dev.off()


fn <- paste0(site.filetag,'_TrendMaps5y_ByOutfall_',fdate,'.pdf')
pdf(file=fn,w=11,h=8.5)
outfall.trendmap.plots.5y <- outfall.trend.map.5y %>% group_by(coc) %>% do(plot={
  tmap <- .
  tcoc <- tmap$coc[1]
  print(tcoc)
  hdr <- paste0('Trend Map by Outfall for ',tcoc,' Since 2017')
  trend_map_plot(tmap,period='2017-01-01',hdr=hdr,gis=F,gis.lst=NULL,bex=0.01,interact=F,wlabs=T)
}) %>% ungroup()
dev.off()


# ONE-WAY COMPARISONS ----
outfall.oneway.tests <- f %>% group_by(coc) %>% do({
  td <- . 
  tcoc <- td$coc[1]
  print(tcoc)
  oneway_mc(td,x=locid)
}) %>% ungroup()


fn <- paste0(site.filetag,'_OneWayTests_ByOutfall_',fdate,'.pdf')
pdf(file=fn,w=11,h=8.5)
outfall.oneway.plots <- outfall.oneway.tests %>% group_by(coc) %>% do(plot={
  td <- .
  tcoc <- td$coc[1]
  hdr <- paste0('Pairwise Mean Difference Contrasts for ',tcoc)
  print(oneway_mc_plot(td,hdr))
}) %>% ungroup()
dev.off()


outfall.oneway.tests.5y <- f %>% group_by(coc) %>% filter(date >= ymd('2017-01-01')) %>% do({
  td <- . 
  tcoc <- td$coc[1]
  print(tcoc)
  oneway_mc(td,x=locid)
}) %>% ungroup()


fn <- paste0(site.filetag,'_OneWayTests5y_ByOutfall_',fdate,'.pdf')
pdf(file=fn,w=11,h=8.5)
outfall.oneway.plots.5y <- outfall.oneway.tests.5y %>% group_by(coc) %>% do(plot={
  td <- .
  tcoc <- td$coc[1]
  hdr <- paste0('Pairwise Mean Difference Contrasts for ',tcoc,' Since 2017')
  print(oneway_mc_plot(td,hdr))
}) %>% ungroup()
dev.off()

fn <- paste0(site.filetag,'_oneway_tests_',fdate,'.csv')
write_excel_csv(outfall.oneway.tests,file=fn)

fn <- paste0(site.filetag,'_oneway_tests_5y_',fdate,'.csv')
write_excel_csv(outfall.oneway.tests.5y,file=fn)



# CI BANDS, PCT REDUCTIONS ----
outfall.pair.cband <- f %>% group_by(coc,locid) %>% do({
  tcoc <- .$coc[1]; tloc <- .$locid[1]
  print(paste(tcoc,tloc))
  ci_band_lm(.,clev=0.99,side='both')
}) %>% ungroup()


## plot confidence bands
fn <- paste0(site.filetag,'_cband_',fdate,'.pdf')
pdf(file=fn,w=11,h=8.5)
outfall.cband.plots <- outfall.pair.cband %>% group_by(coc) %>% do(plot={
  tcoc= .$coc[1]
  print(tcoc)
  hdr <- paste0('Confidence Bands by Outfall for ',tcoc,' With 99% Confidence')
  td <- f %>% filter(coc==tcoc)
  ci_band_plot(td,.,hdr,show.lim=F,pt.size = 2,lab.size = 3)
})
dev.off()


## Last 5 year bands and plots
outfall.pair.cband.5y <- f %>% filter(date >= ymd('2017-01-01')) %>% group_by(coc,locid) %>% do({
  tcoc <- .$coc[1]; tloc <- .$locid[1]
  print(paste(tcoc,tloc))
  ci_band_lm(.,clev=0.99,side='both')
}) %>% ungroup()



fn <- paste0(site.filetag,'_cband5y_',fdate,'.pdf')
pdf(file=fn,w=11,h=8.5)
outfall.cband.plots.5y <- outfall.pair.cband.5y %>% group_by(coc) %>% do(plot={
  tcoc= .$coc[1]
  print(tcoc)
  hdr <- paste0('Confidence Bands by Outfall Since 2017 for ',tcoc,' With 99% Confidence')
  td <- f %>% filter(coc==tcoc,date >= ymd('2017-01-01'))
  ci_band_plot(td,.,hdr,show.lim=F,pt.size = 2,lab.size = 3)
})
dev.off()


## pct changes across range of each coc-outfall pair
outfall.pctdiff <- f %>% group_by(coc,locid) %>% do({
  tcoc <- .$coc[1]; tloc <- .$locid[1]
  print(paste(tcoc,tloc))
  pctdiff_boot(.,side = 'both',impute = 'MC',ns=20,nb=100)
})


outfall.pctdiff.5y <- f %>% filter(date >= ymd('2017-01-01')) %>% group_by(coc,locid) %>% do({
  tcoc <- .$coc[1]; tloc <- .$locid[1]
  print(paste(tcoc,tloc))
  pctdiff_boot(.,side = 'both',impute = 'MC',ns=20,nb=100)
})


fn <- paste0(site.filetag,'_outfall_pctdiff_',fdate,'.csv')
write_excel_csv(outfall.pctdiff,file=fn)

fn <- paste0(site.filetag,'_outfall_pctdiff_5y_',fdate,'.csv')
write_excel_csv(outfall.pctdiff.5y,file=fn)





# SAVE PROJECT ----
fdate <- format.Date(Sys.Date(),format='%y%m%d')
fn <- paste0(tolower(site.filetag),'_',fdate,'.rda')
save.image(fn)

