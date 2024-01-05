# pull-sort-unique function

pus <- function(f,var) {
  library(tidyverse)
  
  f %>% pull({{var}}) %>% unique() %>% sort()
}