# UG package v1.0
# Author: Kirk Cameron, PhD, MacStat Consulting, Ltd
# Last revised: April 22, 2021, rv04

# Purpose: Transform scaled estimates

# Revisions: rv05 better handles negative values;
#  rv04 updates code to be consistent with related functions;
#  rv03 makes transformations consistent with other functions;
#  rv02 eliminates inverse transform



trans <- function(y,tlab='Normal') {
  library(tidyverse)
  
  tf.df <- tibble(tf=c('log','0.125','0.1428571','0.1666666','0.2','0.25','0.3333333','0.5','1','2','3','4','5','6','7','8'),tf.lab=c('Log','Eighth Root','Seventh Root','Sixth Root','Fifth Root','Fourth Root','Cube Root','Square Root','Normal','Square','Cube','Fourth Power','Fifth Power','Sixth Power','Seventh Power','Eighth Power'))
  
# pull selected transformation
  tf <- tf.df %>% filter(tf.lab == tlab) %>% pull(tf)
  pow <- ifelse(tf=='log',NA_real_,as.numeric(tf))

# Transform scaled estimate or vector of scaled estimates
  if (tlab == 'Log') {
    yt <- ifelse(y < 0,NA_real_,log(y))
  } else if (tlab %in% c('Eighth Root','Sixth Root','Fourth Root','Square Root')) {
    yt <- ifelse(y < 0,NA_real_,ifelse(y == 0,0,exp(pow*log(abs(y)))))
  } else {
    yt <- ifelse(y == 0,0,exp(pow*log(abs(y)))*sign(y)^max(pow,(1/pow)))
  }
  yt
}
