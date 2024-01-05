# Author: Kirk Cameron, PhD, MacStat Consulting, Ltd
# Last Update: December 15, 2022, rv01

# Purpose: Plot Tukey HSD contrasts resulting from one-way permutation test

# Revisions: none

oneway_mc_plot <- function(oneway.f,hdr) {
  
  library(tidyverse)
  
# reusable theme
  theme_contrast.plot <- theme_bw() +
    theme(axis.text.y = element_text(size = rel(.75)),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size = rel(.75)),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(linewidth = 0.5,linetype = 'dotted'),
          panel.grid.minor.x = element_blank())
  
  
# create plot
  p <- ggplot(oneway.f,aes(x=contrast,y=diff))+ geom_pointrange(aes(ymin=lwr,ymax=upr),color='blue') + theme_contrast.plot
  p <- p + geom_hline(yintercept = 0,linetype='dotted')
  p <- p + labs(x='',y='Differences in Mean Levels',title = hdr) + coord_flip()
  p
  
  
}