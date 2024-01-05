# Author: Kirk Cameron, PhD, MacStat Consulting, Ltd
# Last Update: December 6, 2022, rv02

# Purpose: Construct exploratory time series, box, or violin plot(s)

# Revisions: rv02 minor update


eplot <- function(f,thead,etype=c('ts','box','violin'),vx=NULL,vlim=NULL,vfill=NULL,vfacet=locid,showlines=T,grid=T,logscale=F,loess=F,rotate.lab=T,pt.size=3,lab.size=3,...) {
  library(tidyverse)
  
  etype <- match.arg(etype)
  
  if (etype != 'ts') {
# set seed for Monte Carlo draw; replace NDs with Monte Carlo imputations
    set.seed(41217002)
    f <- f %>% mutate(N=n(),conc.r=ifelse(nd==0,conc,conc*runif(N)))
  }


  if (!is.null(substitute(vlim))) {
    lim.lab <- toupper(quo_name(enquo(vlim)))
    lim <- f %>% pull({{vlim}})
    if (etype != 'ts') {
      clvec <- c(lim[1],f$conc,f$conc.r)
    } else if (etype == 'ts') {
      clvec <- c(lim[1],f$conc)
    }
  } else if (is.null(substitute(vlim))) {
    if (etype != 'ts') {
      clvec <- c(f$conc,f$conc.r)
    } else if (etype == 'ts') {
      clvec <- c(f$conc)
    }
  }
  

  if (etype == 'ts') {
    x.lab <- 'Sampling Date'
    xx <- f %>% pull(date)
    f <- f %>% mutate(xx=xx)
  } else if (etype != 'ts' && !is.null(substitute(vx))) {
    x.lab <- toupper(quo_name(enquo(vx)))
    xx <- f %>% pull({{vx}})
    f <- f %>% mutate(xx=xx)
  } else if (etype != 'ts' && is.null(substitute(vx))) {
    x.lab <- 'ALL DATA'
    f <- f %>% mutate(xx='All_Data')
  }
    

  if (!is.null(substitute(vfill))) {
    fill.lab <- toupper(quo_name(enquo(vfill)))
    fill.vec <- f %>% pus({{vfill}})
  }
  
  
# only create plot if at least 3 observations; otherwise return NULL
  n <- nrow(f)
  if (n < 3) return(NULL)
  
  f <- f %>% mutate(nd2= case_when(nd == 0~'Detect',nd==1~'ND'))
  
# adjust plot range as necessary
  any.nd <- sum(f$nd,na.rm=T) > 0
  low.lev <- max(clvec,na.rm=T) < 100
  all.same <- var(clvec,na.rm=T) == 0
  
  if (!low.lev && !all.same) {
    crng <- c(min(clvec,na.rm=T),min(clvec,na.rm=T) + 1.2*diff(c(min(clvec,na.rm=T),max(clvec,na.rm=T))))
  } else if (low.lev && !all.same) {
    crng <- c(min(clvec,na.rm=T),min(clvec,na.rm=T) + 1.5*diff(c(min(clvec,na.rm=T),max(clvec,na.rm=T))))
  } else if (all.same) {
    crng <- c(0,min(clvec,na.rm=T))
  }
  
# plot and label (regulatory) limit if available
  if (!is.null(substitute(vlim))) {
# check for missing limit values even if variable exists
    lim.na <- is.na(lim[1])
# adjust plot y-range to fit missing limit text label at bottom if needed
    if (lim.na) crng[1] <- crng[1] - 0.1*diff(crng)
# configure text label as to whether or not <lim> exists
    text_label <- NA_character_
    if (lim.na) text_label <- 'No Limit Available'
# create label for lim when it exists
    text_limit <- NA_character_
    if (!lim.na) text_limit <- paste0(lim.lab,' = ',as.character(lim[1]),' ',f$units[1])
  }

  

# time series plots
  if (etype == 'ts') {
    
# adjust span for Loess as needed depending on grid facet sample size
    if (grid && loess) {
      f <- f %>% group_by({{vfacet}}) %>% mutate(nn=length(conc)) %>% ungroup() %>% mutate(nspan=case_when(nn < 10~1,nn >= 10~0.75,TRUE~NA_real_))
    }

# adjust plot range for Loess if needed
    if (!loess) {
      p <- ggplot(data=f,aes(date,conc)) + ylim(crng)
    } else if (loess) {
      p <- ggplot(data=f,aes(date,conc))
    }
    
    p <- p + labs(title=thead,x=x.lab,y=paste0('Concentration (',as.character(f$units[1]),')'))
    p <- p + labs(shape='',color='') + theme_bw()

  
    if (showlines) p <- p + geom_line()
    
# check for all NDs or all detects; adjust symbol plotting accordingly
    if (all(f$nd==1)) {
      p <- p + geom_point(data=f,aes(color='ND'),size=pt.size,shape=0) + scale_color_manual('',breaks=c('ND'),values=c('red'))
    } else if (all(f$nd==0)) {
      p <- p + geom_point(data=f,aes(color='Detect'),size=pt.size,shape=16) + scale_color_manual('',breaks=c('Detect'),values=c('black'))
    } else {
      p <- p + geom_point(aes(shape=nd2,color=nd2),data=f,size=pt.size) + scale_shape_manual(values=c(16,0)) + scale_color_manual(values=c(Detect='black',ND='red'))
    }
    
# set up faceting arrangement based on number of plots to display if grid = T
    if (grid && !is.null(substitute(vfacet))) {
      tx <- f %>% pus({{vfacet}})
      ftot <- length(tx)
      nc <- ceiling(sqrt(ftot))
      p <- p + facet_wrap(vars({{vfacet}}),ncol= nc)
    }
    
# box plots and violin plots
  } else if (etype != 'ts') {
    
    p <- ggplot(data=f,aes(x=xx,y=conc.r)) + ylim(crng)
    
    p <- p + labs(title=thead,x=x.lab,y=paste0('Concentration (',as.character(f$units[1]),')'))
    p <- p + labs(shape='',color='') + theme_bw()
    if (x.lab == 'ALL DATA') p <- p + theme(axis.text.x = element_blank())
    
    if (is.null(substitute(vfill))) {
      if (etype == 'box') {
        p <- p + geom_boxplot(fill='lightgreen',...,outlier.size = 0,outlier.stroke = 0)
      } else if (etype == 'violin') {
        p <- p + geom_violin(draw_quantiles = c(0.25,0.5,0.75),fill='lightblue',scale = 'width')
      }
    } else if (!is.null(substitute(vfill))) {
      if (etype == 'box') {
        p <- p + geom_boxplot(aes(fill = {{vfill}}),...,outlier.size = 0,outlier.stroke = 0)
      } else if (etype == 'violin') {
        p <- p + geom_violin(aes(fill= {{vfill}}),draw_quantiles = c(0.25,0.5,0.75),scale = 'width')
      }
      
      
      if (length(fill.vec) > 3) {
        p <- p + labs(fill=fill.lab)
      } else if (length(fill.vec) == 2) {
        p <- p + scale_fill_manual(breaks=rev(fill.vec),values = c('lightblue','wheat')) + labs(fill=fill.lab)
      } else if (length(fill.vec == 3)) {
        p <- p + scale_fill_manual(breaks=rev(fill.vec),values = c('lightgreen','lightblue','wheat')) + labs(fill=fill.lab)
      }
    }
    
# check for all NDs or all detects; adjust symbol plotting accordingly
    if (all(f$nd==1)) {
      p <- p + geom_jitter(data=f,aes(y=conc,color='ND'),size=pt.size,shape=0,width = 0.15,height=0) + scale_color_manual('',breaks=c('ND'),values=c('red'))
    } else if (all(f$nd==0)) {
      p <- p + geom_jitter(data=f,aes(y=conc,color='Detect'),size=pt.size,shape=16,width = 0.15,height=0) + scale_color_manual('',breaks=c('Detect'),values=c('black'))
    } else {
      p <- p + geom_jitter(aes(y=conc,shape=nd2,color=nd2),data=f,size=pt.size,width = 0.15,height=0) + scale_shape_manual(values=c(16,0)) + scale_color_manual(values=c(Detect='black',ND='red'))
    }

  }
    
    
    
    
# if <lim> exists, add to plot
  if (!is.null(substitute(vlim)) && !lim.na) p <- p + geom_hline(aes(yintercept=lim[1]),linewidth=1.25,color='blue')

  
# add Loess trend to time series plots if requested
  if (etype == 'ts') {
    if (loess && grid) {
      p <- p + geom_smooth(aes(span=nspan),color='slateblue3',fill='gray85',alpha=0.5) + theme_bw()
    } else if (loess && !grid) {
      p <- p + geom_smooth(color='slateblue3',fill='gray85',alpha=0.5) + theme_bw()
    }
  }
  
# rotate axis date labels if requested
  if (rotate.lab) p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
# display on log base 10 scale if requested
  if (logscale) p <- p + scale_y_log10()
  

# create tibble to hold labels as needed
  if (!is.null(substitute(vlim))) {
    if (etype == 'ts') {
      f_label <- tibble(x1=min(xx),y1=crng[1],label1=text_label,x2=max(xx),y2=lim[1],label2=text_limit)
    } else if (etype != 'ts') {
      tx <- unique(xx); nn <- length(tx)
      f_label <- tibble(x1=tx[1],y1=crng[1],label1=text_label,x2=tx[nn],y2=lim[1],label2=text_limit)
    }
    if (!is.na(f_label$label1)) p= p + geom_text(data=f_label,aes(x=x1,y=y1,label=label1),hjust=0,size=lab.size)
    if (!is.na(f_label$label2)) p= p + geom_text(data=f_label,aes(x=x2,y=y2,label=label2),hjust=1,vjust=-0.5,size=lab.size)
  }
  

  p
}
