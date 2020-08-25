library(dplyr)
nREP = 100
SEP = 3

AA <- matrix(1:nREP, nrow=max(SEP), byrow=TRUE) %>% 
  as.numeric %>% 
  ifelse( duplicated(.), NA, . ) %>% 
  matrix(nrow=max(SEP)) %>% 
  .[sep,] %>% 
  na.omit %>% 
  as.numeric
