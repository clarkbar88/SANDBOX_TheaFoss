# Author: Kirk Cameron, PhD, MacStat Consulting, Ltd
# Last Update: December 15, 2022, rv01

# Purpose: Compute one-way permutation test and multiple comparisons with or without NDs

# Revisions: none

oneway_mc <- function(f,x=locid,y=conc,nd=nd,ns=50,...) {
  
  library(lmPerm)
  library(stats)
  library(rlang)
  library(purrr)
  
  ys <- ensym(y); nds <- ensym(nd)
  xs <- ensym(x)
  
  N <- nrow(f); nd.pct <- with(f,inject(round(100*mean(!!nds,na.rm=T),1)))
  
  out.tmp <- map_dfr(1:ns,.f=function(.x,xf,...) {
    xf <- xf %>% mutate(y.r = mc_vector(!!ys,!!nds,...))
    out1 <- inject(aovp(y.r ~ !!xs,data = xf))
    tmp <- summary(out1)[[1]]
    pval <- tmp[1,5]
    F.stat <- tmp[1,3]/tmp[2,3]
    df.num <- tmp[1,1]; df.denom <- tmp[2,1]
    out2 <- TukeyHSD(out1)[[1]] %>% as_tibble(rownames='contrast') %>% rename(p.adj=`p adj`)
    tibble(N,nd.pct,F.stat,df.num,df.denom,pval,out2)
  },xf=f,...)
  
  out <- out.tmp %>% group_by(contrast) %>% summarise(across(c(N:pval,diff:p.adj),~ mean(.x,na.rm=T)),.groups = 'drop') %>% relocate(contrast,.before = diff)
  out
  
}