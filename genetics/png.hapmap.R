png.hapmap <- function(x, cutoff.hetero=0.2, cutoff.missing=0.2, cutoff.HWE=10e-6){
  print( "If case-control status is available, limit the filtering of cutoff.HWE to control group as a violation in case group may be an indication of association." )
  # Initial of function -----------------------------------------------------  
  # 1.Sorting -----------------------------------------------------------------
  
  x <- as.matrix(x)
  ord.x <- gtools::mixedorder(colnames(x[,-(1:11)]))
  
  
  x <- x[,c(1:11, ord.x+11)]
  
  pvalue.HWE <- sapply( 1:nrow(x), function(indx) {
    input <- strsplit(as.character(unlist(x[indx,-(1:11)]) %>% 
                                     {replace(., list = (. %in% c("NN", "00", "--", "//", "++", "XX") ), NA)} ),"")
    input.genotype <- sapply( input, paste0, collapse="/") %>% 
      genetics::genotype()
    if(length(table(input.genotype))>1 & length(levels(input.genotype))<6){
      input.genotype %>%
        genetics::HWE.exact() %>% .$p.value
    } else {
      1
    }
  })
  rm.HWE <- which( pvalue.HWE <= cutoff.HWE )
  
  rm.missing <- which( x[,-(1:11)] %>% apply(1, function(x) mean(is.na(x) | x %in% c("NN", "00", "--", "//", "++", "XX") )) > cutoff.missing ) # 2804
  rm.missing.sample <- which( x[,-(1:11)] %>% apply(2, function(x) mean(is.na(x) | x %in% c("NN", "00", "--", "//", "++", "XX") )) > cutoff.missing ) # 1
  
  set.heterozygote <- apply(expand.grid(c("T","C","A","G"), c("T","C","A","G")) %>% 
                              filter(Var1 != Var2), 1, paste0, collapse="")
  
  png.heterozygousCalls <- function(x){
    x %>% replace( ., .=="NN", NA ) %>% unlist %>% { mean( . %in% set.heterozygote ) }
  }
  
  rm.hetero <- which( x[,-(1:11)] %>% apply(1, function(x) png.heterozygousCalls(x)) > cutoff.hetero ) # 37
  
  print( x[,-(1:11)][ rm.hetero , ] %>% apply(1, table) )
  
  rm.dup.var <- which(duplicated(x)) # There are no duplicates in variants
  rm.dup.sample <- which(duplicated(t(x[,-(1:11)]))) # There are no duplicates in samples
  
  rm.var <- c(rm.missing, rm.hetero, rm.missing.sample, rm.dup.var, rm.HWE)
  rm.sample <- rm.dup.sample
  
  x.removed <- x
  if( length(rm.var)>0 ) x.removed <- x[ -unique( rm.var ), ]
  if( length(rm.sample)>0 ) x.removed <- x[, -(unique( rm.sample )+11) ]
  
  write.table( x[ -unique( rm.var ), ] )
  write.table( x[, -(unique( rm.sample )+11) ] )
  
  print( dim(x) ) # 49683 x 395
  print( dim(x.removed) ) # 46852 x 395
  
  myX <- as.matrix( rbind( colnames(x.removed), x.removed ) )
  myGD <- apply(myX[-1,-(1:11)], 1,
                function(one) GAPIT.Numericalization(one, bit=2, impute="None", Major.allele.zero=TRUE)) %>% 
    # function(one) GAPIT.Numericalization(one, bit=2, impute="Major", Major.allele.zero=TRUE)) %>% 
    as.data.frame(stringsAsFactors = FALSE)
  myGM <- myX[-1,1:4]
  # myGT <- myX[,c(12:ncol(myX))]
  
  set.seed(120120)
  myX.impute <- rbind( myX[1,], cbind( myX[-1,(1:11)], t( apply(t(myX[-1,-(1:11)]), 2, png.impute.snp) ) ) )
  set.seed(120120)
  myGD.impute <- apply(myGD, 2, png.impute.numeric)
  
  Taxa <- as.character( myX[1, -(1:11)] )
  out = list (myX=as.data.frame(myX), 
              myX.impute=as.data.frame(myX.impute), 
              myGD=data.frame(Taxa=Taxa, myGD), 
              myGD.impute=data.frame(Taxa=Taxa, myGD.impute), 
              myGM=myGM)#, myGT=myGT)
  
  out
}



