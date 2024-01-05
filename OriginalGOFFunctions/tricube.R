# Author: Kirk Cameron, MacStat Consulting, Ltd
# Last Update: November 28, 2023, rv04

# Purpose: compute tri-cube weights for given dataset

# Revisions: rv04 minor update;
#  rv03 updates code;
#  rv02 adds check for vectors of length 1 or less


tricube <- function(x,eps=1e-5) {
  
## if x has length <= 1, simply return weight = 1
	if (length(x) <= 1) return(1)
  
## compute upper bound as either max(x)*(1+eps)
## or as 1 if max(x) < eps; re-scale x by upper bound
	x <- abs(x)
	xmax <- ifelse(max(x) < eps,1,max(x)*(1+eps))
	wt <- (1 - (x/xmax)^3)^3
	wt
}
