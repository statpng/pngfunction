
# R base function
pb <- txtProgressBar(min=0, max=nrow(genotype.hapmap[-1,]), style=3)
setTxtProgressBar(pb, i)
