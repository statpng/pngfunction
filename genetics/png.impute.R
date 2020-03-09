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

                          
                          
                          
                          
                          
## Using "synbreed" package
                          
mat.order <- function(mat, row=TRUE, col=TRUE){
  library(gtools)
  if( row ){
    row.order <- mixedorder(rownames(mat))
  } else {
    row.order <- TRUE
  }
  
  if( col ){
    col.order <- mixedorder(colnames(mat))
  } else {
    col.order <- TRUE
  }
  return( mat[row.order, col.order] )
}


hapmap.impute <- function(hapmap.missing.df){
  library(synbreed)
  
  header <- hapmap.missing.df[1,]
  snp.info <- hapmap.missing.df[-1,c(1:11)]
  
  NN.df <- t( hapmap.missing.df[-1,-c(1:11)] )
  missing.df <- ifelse(NN.df == "NN", NA, NN.df)
  
  gp.created <- create.gpData(geno=missing.df)
  gp.coded <- codeGeno(gp.created, 
                       impute=TRUE, 
                       impute.type = "random", 
                       # label.heter=function(x) substr(x,1,1)!=substr(x,3,3),
                       label.heter=NULL,
                       verbose=TRUE)
  gp.coded$geno
  print( lapply( summary(gp.coded), head ) )
  imputed.df <- mat.order( gp.coded$geno )
  
  removed.snp <- which( apply(missing.df, 2, function(x) all(is.na(x))) )
  
  # print(
  #   missing.df[,-removed.snp][1:10,unique( which( is.na(missing.df[,-removed.snp]), arr.ind=TRUE )[,2] )]
  # )
  # print(
  #   imputed.df[1:10,unique( which( is.na(missing.df[,-removed.snp]), arr.ind=TRUE )[,2] )]
  # )
  
  # data.frame(
  #   missing=missing.df[,-removed.snp][ which( is.na(missing.df[,-removed.snp]) ) ],
  #   imputed=imputed.df[ which( is.na(missing.df[,-removed.snp]) ) ]
  # )
  
  gp.map <- mat.order( gp.coded$alleles, row=TRUE, col=FALSE )
  
  hapmap.df <- imputed.df
  for( h in 1:ncol(imputed.df) ){
    hh <- imputed.df[,h]
    hapmap.df[,h] <- unlist( gp.map[h, hh+2] )
  }
  
  print(
    table( data.frame( numeric=imputed.df[,1], genotype=unlist(gp.map[1, imputed.df[,1]+2]) ) )
  )
  
  h=1
  print(
    table( data.frame( numeric=imputed.df[,h], genotype=hapmap.df[,h] ) )
  )
  
  if(length(removed.snp)>0){
    out <- rbind(header, data.frame(snp.info[-removed.snp,], t(hapmap.df) ))
  } else {
    out <- rbind(header, data.frame(snp.info, t(hapmap.df) ))
  }
  
  print(paste0("There are no missing values(", sum(is.na(out[-1,-c(1:11)])),")."))
  return(out)
}
