library(dplyr)
II <- 1

nREP = 100
SEP = 1:3
SET = expand.grid(SEP=SEP)
sep = SET$SEP[II]

AA <- matrix(1:nREP, nrow=max(SEP), byrow=TRUE) %>% 
  as.numeric %>% 
  ifelse( duplicated(.), NA, . ) %>% 
  matrix(nrow=max(SEP)) %>% 
  .[sep,] %>% 
  na.omit %>% 
  as.numeric
