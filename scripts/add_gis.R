# Geostatistical Temporal/Spatial software (GTS), v1.1
# Flowchart or Task: B2.1, COC Analysis/Maps of Exceedance Status
# Author -- Kirk Cameron, MacStat Consulting, Ltd
# Last Update -- September 26, 2017, rv07

# R function for adding specific GIS layers to annotate GTS maps/plots

# Revision Notes: rv06 adds ability to clip shape file layers to bounding box around
#  convex hull of plotted site area;
#  rv05 [mjk] bug fixes to rv04
#  rv04 builds in check to see if there are current plotting limits
#  and whether they should be expanded to accommodate shape file extent;
#  rv03 changes inputs to add GIS parameter dataframe and converted shapefiles;
#  rv02 alters function name and inputs to generalize routine
# External Calls: libraries <tidyverse>, <RColorBrewer>, <raster>, <broom>
# Defaults: none
# Inputs: 
# Outputs: 


add_gis= function(plot,gis.lst) {
  library(tidyverse)
  library(RColorBrewer)
  library(raster)
  library(broom)

# compute plot limits; then combine into 2x2 matrix
  xr= range(plot$data$east)
  yr= range(plot$data$north)
  rmat= matrix(c(xr,yr),nrow=2,ncol=2,byrow = T)
  
# set range of possible colors for shape layers
  bluegreen= brewer.pal(9, 'BuGn')
  redpurple= brewer.pal(9,'PuRd')
  gray= brewer.pal(9,'Greys')
  green= brewer.pal(9,'Greens')
  orange= brewer.pal(9,'Oranges')
  blue= brewer.pal(9,'Blues')
  
  N= gis.lst$n.shp
  for (i in 1:N) {
    if (gis.lst$gis.par$use[i]) {
      color.vec= eval(as.name(gis.lst$gis.par$color[i]))
      tlayer= crop(gis.lst$shp.lst[[i]],rmat)
      tdf= broom::tidy(tlayer)
      if (gis.lst$gis.par$type[i]=='pts') {
        idx= which(names(tdf)==c('coords.x1','coords.x2'))
        names(tdf)[idx]= c('long','lat')
        tdf$group= 1:nrow(tdf)
      }
      if (gis.lst$gis.par$type[i]=='path') {
        plot= plot + geom_path(aes(x=long,y=lat,group=group),data=tdf,color=color.vec[5],alpha=gis.lst$gis.par$alpha[i])
      } else if (gis.lst$gis.par$type[i]=='poly') {
        plot= plot + geom_polygon(aes(x=long,y=lat,group=group),data=tdf,color=color.vec[5],fill=color.vec[3],alpha=gis.lst$gis.par$alpha[i])
      } else if (gis.lst$gis.par$type[i]=='pts') {
        plot= plot + geom_point(aes(x=long,y=lat,group=group),data=tdf,color=color.vec[5],alpha=gis.lst$gis.par$alpha[i])
      }
    }
  }

  plot
}

