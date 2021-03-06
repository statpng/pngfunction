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
  # pi(p) ~ beta(2, 2)
  # L(y|p) ~ B(p)
  # pi(p|y) \prop pi(p) * L(y|p)
  #         ~ Beta(y+2, n+2-y)
  maf <- rbeta(n=length(x.na), shape1=y+2, n+2-y )
  # curve(dbeta(x, shape1=y+2, n+2-y ))
  impute.value <- rbinom(n=length(x.na), size = 2, prob = maf)
  xx[is.na(xx)] <- names(tb.new)[ impute.value+1 ]
  xx
}
