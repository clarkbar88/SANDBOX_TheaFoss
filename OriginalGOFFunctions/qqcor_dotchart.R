# Author: Kirk Cameron, PhD, MacStat Consulting, Ltd
# Last update: March 15, 2023, rv03

# Purpose: Cleveland dotchart of Q-Q plot correlations

# Revisions: rv03 changes name from <gg_dotchart>; updates code;
#  rv02 minor update; also fixes x range



qqcor_dotchart <- function(x,y,hdr) {
  
  library(tidyverse)
  
  if (missing(hdr)) hdr <- 'TEST'
  
  f <- tibble(x,y) |> mutate(y=factor(y,levels = y,ordered = T))
  
## create subset of highest correlation case(s)
  ix <- which(x == max(x,na.rm = T))
  f.tmp <- f |> slice(ix)


## reusable plot theme
theme_dotplot <- theme_bw() +
  theme(axis.text.y = element_text(size = rel(.75)),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = rel(.75)),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.5,linetype = 'dotted'),
        panel.grid.minor.x = element_blank())


## create plot with <ggplot2>
p <- ggplot(f,aes(x,y)) + xlim(0,1) + geom_point(color='blue',size=3,shape=1) + theme_dotplot
p <- p + labs(x='Q-Q Plot Correlation',y='',title = hdr)
p <- p + geom_point(data = f.tmp,aes(x,y),color='purple',size=3,shape=16)
p
}
