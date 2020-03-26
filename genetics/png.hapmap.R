devtools::source_url("https://raw.githubusercontent.com/statpng/pngfunction/master/genetics/png.impute.R")

png.impute.numeric <- function(xx){
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
  maf <- rbeta(n=length(x.na), shape1=y+2, n+2-y )
  impute.value <- rbinom(n=length(x.na), size = 2, prob = maf)
  xx[is.na(xx)] <- impute.value
  xx
}

png.impute.snp <- function(xx){
  
  xx <- ifelse(xx=="NN", NA, xx)
  x.na <- xx[is.na(xx)]
  x.value <- xx[!is.na(xx)]
  
  
  tb <- table( x.value )
  alleles <- unique( unlist( strsplit(names(tb), "") ) )
  major.allele <- unique( unlist( strsplit( names( tb[which.max(tb)] ), "" ) ) )
  minor.allele <- alleles[ !alleles %in% major.allele ]
  
  if(length(minor.allele)==0){
    xx[is.na(xx)] <- names(tb)
    return(xx)
  }
  
  combs <- unique( apply( expand.grid( alleles, alleles ), 1, function(x) paste0(sort(x), collapse="" ) ) )
  
  tb <- table( factor( x.value, levels=combs ) )
  ord.tb <- order( sapply( strsplit( names(tb), "" ), function(x) sum(x==minor.allele) ) )
  tb.new <- tb[ord.tb]
  
  y <- sum( tb.new[2] + 2*tb.new[3] )
  n <- 2*sum(tb.new)
  # maf <- (tb[2]+2*tb[3])/(2*sum(tb))
  # pi(p) ~ beta(0.5, 0.5)
  # L(y|p) ~ B(p)
  # pi(p|y) \prop pi(p) * L(y|p)
  #         ~ Beta(y+0.5, n+0.5-y)
  maf <- rbeta(n=length(x.na), shape1=y+0.5, n+0.5-y )
  # curve(dbeta(x, shape1=y+2, n+2-y ))
  impute.value <- rbinom(n=length(x.na), size = 2, prob = maf)
  xx[is.na(xx)] <- names(tb.new)[ impute.value+1 ]
  xx
}
                          

png.hapmap <- function(x, cutoff.hetero=0.2, cutoff.missing=0.2, cutoff.HWE=10e-6, write=FALSE){
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

  if( write ){
    write.table( x[ -unique( rm.var ), ] )
    write.table( x[, -(unique( rm.sample )+11) ] )
  }  
  
  print( dim(x) ) # 49683 x 395
  print( dim(x.removed) ) # 46852 x 395
  
  myX <- as.matrix( x.removed )
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
