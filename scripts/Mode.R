# Author: Kirk Cameron, PhD, MacStat Consulting, Ltd
# Last Update: November 29, 2023, rv05

# Purpose: compute statistical mode of a vector 
# (i.e., most frequent single value); can be character or numeric

# Revisions: rv05 minor update;
#  rv04 fixes inconsistency in return types;
#  rv03 finds first mode if more than one exist;
#  rv02 alters function to return mode of non-missing values


Mode <- function(x) {
  if (all(is.na(x))) return(x[1])
  x <- x[!is.na(x)]
  ux <- unique(x)
  tab <- tabulate(match(x,ux))
  ux[tab == max(tab)][1]
}
