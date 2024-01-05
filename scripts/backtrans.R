# UG package v1.0
# Author: Kirk Cameron, PhD, MacStat Consulting, Ltd
# Last revised: February 8, 2022, rv06

# Purpose: Back-transform scaled estimates

# Revisions: rv06 better handles negative values;
#  rv05 adds other edge cases;
#  rv04 updates code to be consistent with related functions;
#  rv03 makes transformations consistent with other functions;
#  rv02 eliminates inverse transform



backtrans <- function(yt,tlab='Normal') {
  library(tidyverse)
  
  tf.df <- tibble(tf=c('log','0.125','0.1428571','0.1666666','0.2','0.25','0.3333333','0.5','1','2','3','4','5','6','7','8'),tf.lab=c('Log','Eighth Root','Seventh Root','Sixth Root','Fifth Root','Fourth Root','Cube Root','Square Root','Normal','Square','Cube','Fourth Power','Fifth Power','Sixth Power','Seventh Power','Eighth Power'))
  
# pull selected transformation
  tf <- tf.df %>% filter(tf.lab == tlab) %>% pull(tf)
  pow <- ifelse(tf=='log',NA_real_,as.numeric(tf))

# back-transform scaled estimate or vector of scaled estimates
  if (tlab=='Log') {
    y <- exp(yt)
  } else if (tlab %in% c('Eighth Power','Sixth Power','Fourth Power','Square')) {
    y <- ifelse(is.na(yt) | is.nan(yt) | yt < 0,NA_real_,exp((1/pow)*log(abs(yt)))*sign(yt)^max(pow,(1/pow)))
  } else {
    y <- ifelse(is.na(yt) | is.nan(yt),NA_real_,exp((1/pow)*log(abs(yt)))*sign(yt)^max(pow,(1/pow)))
  }
  y
}
