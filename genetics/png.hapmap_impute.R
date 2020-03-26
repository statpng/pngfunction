devtools::source_url("https://raw.githubusercontent.com/statpng/pngfunction/master/png.mat_order.R")

png.hapmap_impute <- function(hapmap.missing.df){
  # You must read the genotype data with "header=FALSE"
  # dim(hapmap.missing.df) = (p+1) x (n+11)
  library(synbreed)
  
  header <- hapmap.missing.df[1,]
  snp.info <- hapmap.missing.df[-1,c(1:11)]
  
  NN.df <- t( hapmap.missing.df[-1,-c(1:11)] )
  missing.df <- ifelse(NN.df %in% c("NN", "00", "--", "//", "++", "XX"), NA, NN.df)
  
  gp.created <- create.gpData(geno=missing.df)
  gp.coded <- codeGeno(gp.created, 
                       impute=TRUE, 
                       impute.type = "random", 
                       # label.heter=function(x) substr(x,1,1)!=substr(x,3,3),
                       label.heter=NULL,
                       verbose=TRUE)
  gp.coded$geno
  print( lapply( summary(gp.coded), head ) )
  imputed.df <- png.mat_order( gp.coded$geno )
  
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
  
  gp.map <- png.mat_order( gp.coded$alleles, row=TRUE, col=FALSE )
  
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
