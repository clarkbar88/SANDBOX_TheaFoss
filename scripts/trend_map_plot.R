# Author: Kirk Cameron, PhD, MacStat Consulting, Ltd
# Last Update: December 14, 2022, rv15

# Purpose: plot trend map

# Revisions: rv15 alters input parameters; slight name change from <trendmap_plot>;
#  rv14 updates parameters;
#  rv13 updates code;
#  rv12 uses <htmlwidgets> to prepend plot title in interactive mode;
#  rv11 fixes error in plotly title;
#  rv10 minor cleanup;
#  rv09 widens plotly margins; alters angle of easting labels;
#  rv08 adds slope value to plotly hypertext; changes color of convex hull;
#  rv07 minor tweaks; defaults to user.zone='ALL';
#  rv06 adjusts title placement; adjusts legends when interact=T; 
#  rv05 adds <user> inputs; adds <plotly>
#  rv04 renames function and updates dependencies; also adds expansion factor
#  to computation of convex hull
#  rv03 simplifies GIS layer option; also assumes slope tags have been computed externally;
#  rv02 adds in GIS layering function


trend_map_plot <- function(map,period='hist',locid=locid,hdr,gis=F,gis.lst=NULL,bex=0.01,interact=T,wlabs=F) {
  library(tidyverse)
 # library(vcd) isn't used
  library(plotly)
  library(htmlwidgets)
  library(ggrepel)
  library(rlang)
  library(purrr)
  
##  source('add_gis.R')
##  source('convex_hull.R')
  
  locid <- ensym(locid)
  
  loc <- map %>% group_by(!!locid) %>% summarise(east=mean(east,na.rm=T),north=mean(north,na.rm=T),.groups = 'drop')
  hull <- inject(convex_hull(loc,locid=!!locid,bex=bex))
  xrange <- with(hull,range(east)); yrange <- with(hull,range(north))
  xlim.hull <- xrange + c(-1,1)*diff(xrange)*bex; ylim.hull <- yrange + c(-1,1)*diff(yrange)*bex
  
  map <- map %>% mutate(slope.size= as.character(slope.tag))
  
  p <- ggplot(map) + theme_bw() + geom_polygon(data=hull,aes(x=east,y=north),fill=NA,color='blue',size=0.9) + xlim(xlim.hull) + ylim(ylim.hull)
  p <- p  + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  if (!interact) p <- p + labs(title=hdr) + theme(plot.title = element_text(hjust=0))
  
# add GIS layers if available
  if (gis) p <- add_gis(p,gis.lst)
  
  if (!interact) {
    p <- p + geom_point(data=map,aes(label=!!locid,x=east,y=north,color=trend.tag,size=slope.size,shape=trend.tag)) + scale_color_manual(values= c('INCR'='red','PROB-INCR'='pink1','DECR'='blue','PROB-DECR'='skyblue','ND'='green','FLAT'='tan','N/A'='black'),name='Trend Type') + scale_shape_manual(values= c('INCR'=16,'PROB-INCR'=16,'DECR'=16,'PROB-DECR'=16,'ND'=17,'FLAT'=5,'N/A'=1),name='Trend Type') + scale_size_manual(values= c('1'=2.5,'2'=3,'3'=3.5,'4'=4.5,'5'=5.5),name='Slope Size') + labs(x='Easting',y='Northing') + guides(shape=guide_legend(order=1),color=guide_legend(order=1),size=guide_legend(order=2))
  } else if (interact) {
    p <- p + geom_point(data=map,aes(label=!!locid,x=east,y=north,color=trend.tag,shape=trend.tag,
              text = paste('Locid: ',!!locids,'<br>East: ',east,'<br>North: ',north,'<br>Slope: ',round(slope,6),'<br>Slope Size: ',slope.size)),size=3)
    p <- p + scale_color_manual(values= c('INCR'='red','PROB-INCR'='pink1','DECR'='blue','PROB-DECR'='skyblue','ND'='green','FLAT'='tan','N/A'='black')) + scale_shape_manual(values= c('INCR'=16,'PROB-INCR'=16,'DECR'=16,'PROB-DECR'=16,'ND'=17,'FLAT'=5,'N/A'=1)) + labs(x='Easting',y='Northing',color='Trend Type',shape='Trend Type') + guides(shape=guide_legend(order=1),color=guide_legend(order=1))
  }

  if (wlabs && !interact) p <- p + geom_label_repel(aes(label=!!locid,x=east,y=north),color='blue',size=1.5,label.padding = unit(.15,'lines'))
  

# use <plotly> to make plot interactive if desired
  if (interact) {
    p <- ggplotly(p,tooltip = 'text')
    header <- paste0(header,' (Hover for Slope in Units Per Year)')
    p <- htmlwidgets::prependContent(p,htmltools::strong(header))
    p <- p %>% layout(legend=list(y=0.5))
    p$x$layout$annotations[[1]]$y <- 0.65
    p <- p %>% layout(margin=list(l=75,b=75))
  }
  print(p)
  map
}

