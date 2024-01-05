# Author -- Kirk Cameron, PhD, MacStat Consulting, Ltd
# Last Update -- December 16, 2022, rv16

# Purpose: Plot confidence bands with possible limit overlay

# Revisions: rv16 adds option for no limit;
#  rv15 adds <check_overlap> option to avoid label overprinting; also adds lim.type = 'tl';
#  rv14 allows for limits that vary by location;
#  rv13 updates and streamlines code; changes name to <ci_band_plot> from <confband_plot>;
#  rv12 minor update;
#  rv11 adds logscale option;
#  rv10 adds scales option;
#  rv09 adds MCL option;
#  rv08 minor fix to variable names;
#  rv07 combines <confband_plot> and <confband_test_grid_plot> functionality;
#  rv06 changes function name; makes gridding optional; allows limit options;
#  rv05 adds options for label and point size;
#  rv04 updates name; fixes cases where confidence band cannot be computed;
#  also handles cases with two-sided GWPS;
#  rv03 adds option to rotate x-axis labels;
#  rv02 minor cleanup

ci_band_plot <- function(f,fit,hdr,vloc=locid,showlines=T,grid=T,show.lim=T,lim.type=c('gwps','btv','mcl','tl'),x.section=F,rotate.lab=T,lab.size=3,pt.size=3,logscale=F,grid.scale=c('fixed','free','free_x','free_y'),grid.row = F,aspect=NULL,axis.lab.size=11) {
  
  library(tidyverse)
  library(purrr)
  library(rlang)
  
  lim.type <- match.arg(lim.type)
  grid.scale <- match.arg(grid.scale)

  f$nd2 <- with(f,ifelse(nd==0,'Detect','ND'))
  
  if (show.lim) {
    if (lim.type == 'gwps') {
      clvec <- c(f$gwps.lo,f$gwps.hi,f$conc,fit$yhat,fit$lwr,fit$upr)
      lim.lo.na <- all(is.na(f$gwps.lo)); lim.hi.na <- all(is.na(f$gwps.hi))
    } else if (lim.type == 'btv') {
      clvec <- c(f$lbtv,f$ubtv,f$conc,fit$yhat,fit$lwr,fit$upr)
      lim.lo.na <- all(is.na(f$lbtv)); lim.hi.na <- all(is.na(f$ubtv))
    } else if (lim.type == 'mcl') {
      clvec <- c(f$mcl,f$conc,fit$yhat,fit$lwr,fit$upr)
      lim.lo.na <- TRUE; lim.hi.na <- all(is.na(f$mcl))
    } else if (lim.type == 'tl') {
      clvec <- c(f$ltl,f$utl,f$conc,fit$yhat,fit$lwr,fit$upr)
      lim.lo.na <- all(is.na(f$ltl)); lim.hi.na <- all(is.na(f$utl))
    }
  } else if (!show.lim) {
    clvec <- c(f$conc,fit$yhat,fit$lwr,fit$upr)
    lim.lo.na <- T; lim.hi.na <- T
  }
  
  any.nd <- sum(f$nd,na.rm=T) > 0
  all.same <- var(clvec,na.rm=T) == 0
  if (!all.same) {
    crng <- c(min(clvec,na.rm=T),min(clvec,na.rm=T) + 1.1*diff(range(clvec,na.rm=T)))
  } else if (all.same) {
    crng <- c(0,min(clvec,na.rm=T))
  }

# create labels for limits when they exist; add onto data
  if (show.lim) {
    if (lim.type == 'gwps') {
      f <- f %>% mutate(lwr.text_limit=case_when(is.na(gwps.lo)~NA_character_,TRUE~paste0('Lwr GWPS = ',as.character(gwps.lo),' ',units)),upr.text_limit=case_when(is.na(gwps.hi)~NA_character_,TRUE~paste0('Upr GWPS = ',as.character(gwps.hi),' ',units)),z.lo=gwps.lo,z.hi=gwps.hi)
    } else if (lim.type == 'btv') {
      f <- f %>% mutate(lwr.text_limit=case_when(is.na(lbtv)~NA_character_,TRUE~paste0('Lwr BTV = ',as.character(lbtv),' ',units)),upr.text_limit=case_when(is.na(ubtv)~NA_character_,TRUE~paste0('Upr BTV = ',as.character(ubtv),' ',units)),z.lo=lbtv,z.hi=ubtv)
    } else if (lim.type == 'mcl') {
      f <- f %>% mutate(lwr.text_limit=NA_character_,upr.text_limit=case_when(is.na(mcl)~NA_character_,TRUE~paste0('MCL = ',as.character(mcl),' ',units)),z.lo=NA_real_,z.hi=mcl)
    } else if (lim.type == 'tl') {
      f <- f %>% mutate(lwr.text_limit=case_when(is.na(ltl)~NA_character_,TRUE~paste0('LTL = ',as.character(ltl),' ',units)),upr.text_limit=case_when(is.na(utl)~NA_character_,TRUE~paste0('UTL = ',as.character(utl),' ',units)),z.lo=ltl,z.hi=utl)
    }
  } else if (!show.lim) {
    f <- f %>% mutate(z.lo=NA_character_,z.hi=NA_character_)
  }

  f <- f %>% mutate(x1=max(date))
  
  

# approximate key columns from <fit> and add onto <f>
  
  f.tmp <- f %>% group_by({{vloc}}) %>% do({
    td <- .
    tloc <- td %>% slice(1) %>% pull({{vloc}}) %>% as.character() #change from n=1 to 1
    tfit <- fit %>% filter({{vloc}} == tloc)
    if (nrow(tfit) <= 1) {
      td <- td %>% mutate(yhat=NA_real_,lwr=NA_real_,upr=NA_real_)
    } else if (nrow(tfit) > 1) {
      td <- td %>% mutate(yhat=approx(tfit$dd,tfit$yhat,xout=date)$y,lwr=approx(tfit$dd,tfit$lwr,xout=date)$y,upr=approx(tfit$dd,tfit$upr,xout=date)$y)
    }
    td
  }) %>% ungroup()


# if logscale = T, clip plotting region to be non-negative
  if (logscale) {
    clip <- 1e-6
    f.tmp <- f.tmp %>% mutate(yhat=pmax(yhat,clip),lwr=pmax(lwr,clip),upr=pmax(upr,clip))
    crng <- pmax(crng,clip)
  }

# specify plot range if using fixed scale
  if (grid.scale == 'fixed') {
    p <- ggplot(data=f.tmp,aes(x=date)) + ylim(crng) + theme_bw(base_size = axis.lab.size)
  } else {
    p <- ggplot(data=f.tmp,aes(x=date)) + theme_bw()
  }

  p <- p + labs(title=hdr,x='Sampling Date',y=paste('Concentration (',as.character(f.tmp$units[1]),')',sep=''))
  p <- p + labs(color='')
  
# overlay confidence band or confidence band cross-sections
  if (!x.section) {
    p <- p + geom_ribbon(aes(ymin=lwr,ymax=upr),fill='gray80',alpha=0.5) + geom_line(aes(y=yhat),color='brown',linewidth=1.25)
  } else if (x.section) {
    p <- p + geom_line(aes(y=yhat),color='gray80',linewidth=0.75) + geom_linerange(aes(ymin=lwr,ymax=upr),color='red')
  }


# add limits (if non-missing) to plot
  if (!lim.lo.na) p <- p + geom_hline(aes(yintercept=z.lo),linewidth=1,linetype='dashed',colour='blue')
  if (!lim.hi.na) p <- p + geom_hline(aes(yintercept=z.hi),linewidth=1,linetype='dashed',colour='blue')

  
# connect dots if requested
  if (showlines) p <- p + geom_line(aes(y=conc),linewidth=0.5,color='gray80')


# check for all NDs or all detects; adjust symbol plotting accordingly
  if (all(f.tmp$nd == 1)) {
    p <- p + geom_point(data=f.tmp,aes(y=conc,color='ND'),size=pt.size,shape=0) + scale_color_manual('',breaks=c('ND'),values=c('red'))
  } else if (all(f.tmp$nd == 0)) {
    p <- p + geom_point(data=f.tmp,aes(y=conc,color='Detect'),size=pt.size,shape=16) + scale_color_manual('',breaks=c('Detect'),values=c('black'))
  } else {
    p <- p + geom_point(aes(y=conc,shape=nd2,color=nd2),data=f.tmp,size=pt.size) + scale_shape_manual('',values=c(16,0)) + scale_color_manual('',values=c(Detect='black',ND='red'))
  }
  
# rotate axis date labels if desired
  if (rotate.lab) p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))

  
# print labels if non-missing
  if (!lim.lo.na) p <- p + geom_text(aes(x=x1,y=z.lo,label=lwr.text_limit),hjust=1,vjust=1,size=lab.size,check_overlap = T)
  if (!lim.hi.na) p <- p + geom_text(aes(x=x1,y=z.hi,label=upr.text_limit),hjust=1,vjust=-0.5,size=lab.size,check_overlap = T)
  
# display on log base 10 scale if desired
  if (logscale) p <- p + scale_y_log10()
  
# change aspect ratio if requested
  if (!is.null(aspect)) p <- p + theme(aspect.ratio = aspect)
  
# change axis label size if requested
##  if (!is.null(axis.lab.size)) p <- p + theme(axis.text=element_text(size=axis.lab.size))
  vlocs <- ensym(vloc)

  if (grid) {
    ftot <- f.tmp %>% pus({{vloc}}) %>% length()
##    ftot <- length(unique(df.tmp$locid))
    nc <- ceiling(sqrt(ftot))
    if (!grid.row) {
      p <- p + inject(facet_wrap(~!!vlocs,ncol= nc,scales = grid.scale))
    } else {
      p <- p + inject(facet_wrap(~!!vlocs,nrow= nc,scales = grid.scale))
    }
  }

  print(p)
}
