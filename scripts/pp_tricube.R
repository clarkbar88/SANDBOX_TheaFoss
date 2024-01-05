# Author -- Kirk Cameron, PhD, MacStat Consulting, Ltd
# Last Revised -- November 08, 2022, rv02

# Purpose: compute tri-cube weights for probability plot residuals

# Revisions: rv02 changes name from <tricube_pp> to <pp_tricube>


pp_tricube <- function(x,mx=1) {
  
# if x has length <= 1, return weight of 1
	if (length(x) <= 1) return(1)
  
# if max(x) > mx; re-scale x by new upper bnd
	x <- abs(x)
	xmax <- ifelse(max(x) <= mx,mx,max(x) + mx)
	wt <- (1 - (x/xmax)^3)^3
	wt
}
