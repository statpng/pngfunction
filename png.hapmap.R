png.hapmap <- function(x){
  # Initial of function -----------------------------------------------------  
  # 1.Sorting -----------------------------------------------------------------
  ord.x <- gtools::mixedorder(colnames(x[,-(1:11)]))
  
  x <- x[,c(1:11, ord.x+11)]
  
  rm.missing <- which( x[,-(1:11)] %>% apply(1, function(x) mean(is.na(x) | x == "NN")) > 0.2 ) # 2804
  rm.missing.sample <- which( x[,-(1:11)] %>% apply(2, function(x) mean(is.na(x) | x == "NN")) > 0.2 ) # 1
  
  set.heterozygote <- apply(expand.grid(c("T","C","A","G"), c("T","C","A","G")) %>% 
                              filter(Var1 != Var2), 1, paste0, collapse="")
  
  png.heterozygousCalls <- function(x){
    x %>% replace( ., .=="NN", NA ) %>% unlist %>% { mean( . %in% set.heterozygote ) }
  }
  
  rm.hetero <- which( x[,-(1:11)] %>% apply(1, function(x) png.heterozygousCalls(x)) > 0.2 ) # 37
  
  print( x[,-(1:11)][ rm.hetero , ] %>% apply(1, table) )
  
  rm.dup.var <- which(duplicated(x)) # There are no duplicates in variants
  rm.dup.sample <- which(duplicated(t(x[,-(1:11)]))) # There are no duplicates in samples
  
  rm.var <- c(rm.missing, rm.hetero, rm.missing.sample, rm.dup.var)
  rm.sample <- rm.dup.sample
  
  x.removed <- x
  if( length(rm.var)>0 ) x.removed <- x[ -unique( rm.var ), ]
  if( length(rm.sample)>0 ) x.removed <- x[, -(unique( rm.sample )+11) ]
  
  print( dim(x) ) # 49683 x 395
  print( dim(x.removed) ) # 46852 x 395
  
  
  myX <- as.matrix( rbind( colnames(x.removed), x.removed ) )
  myGD <- apply(myX[-1,-(1:11)], 1,
                function(one) GAPIT.Numericalization(one, bit=2, impute="None", Major.allele.zero=TRUE)) %>% 
    as.data.frame(stringsAsFactors = FALSE)
  myGM <- myX[,1:4]
  myGT <- myX[,c(12:ncol(myX))]
  
  
  png.snpimpute <- function(xx){
    x.na <- xx[is.na(xx)]
    x.value <- xx[!is.na(xx)]
    
    tb <- table( factor( x.value, levels=0:2 ) )
    
    y <- sum(x.value)
    n <- 2*sum(tb)
    # maf <- (tb[2]+2*tb[3])/(2*sum(tb))
    # pi(p) ~ beta(2, 2)
    # L(y|p) ~ B(p)
    # pi(p|y) \prop pi(p) * L(y|p)
    #         ~ Beta(y+2, n+2-y)
    maf <- rbeta(n=length(x.na), shape1=y+2, n+2-y )/2
    impute.value <- rbinom(n=length(x.na), size = 2, prob = maf)
    xx[is.na(xx)] <- impute.value
    xx
  }
  
  myGD.impute <- apply(myGD, 2, png.snpimpute)
  
  out = list (myX, myGD, myGD.impute, myGM, myGT)
  
  out
}
