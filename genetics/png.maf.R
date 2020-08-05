png.maf <- function(xx, sep=""){
  
  # xx: vecor with c(GG, GC, CC, CC, CC)
  # This function will return 1e-22 for the SNP with no minor allele
  if(all(is.na(xx))) return(NA)
  if(length(unique(xx))<=1) return(1e-22)
  
  
  xx <- ifelse(xx %in% c("NN", "00", "--", "//", "++", "XX"), NA, xx)
  x.na <- xx[is.na(xx)]
  x.value <- xx[!is.na(xx)]
  
  tb <- table( unlist( unlist( strsplit( x.value, sep ) ) ) )
  maf <- min( prop.table(tb) )
  
  
  if( FALSE ){
    alleles <- names(tb)
    major.allele <- alleles[which.max(tb)]
    minor.allele <- alleles[ !alleles %in% major.allele ]
    
    if(length(minor.allele)==0){
      return(1e-22)
    }
    
    combs <- unique( apply( expand.grid( alleles, alleles ), 1, function(x) paste0(sort(x), collapse=sep ) ) )
    x.value <- sapply( x.value, function(xx) strsplit(as.character(xx), sep) %>% unlist %>% gtools::mixedsort() %>% paste0(collapse=sep) )
    
    tb <- table( factor( x.value, levels=combs ) )
    ord.tb <- order( sapply( strsplit( names(tb), sep ), function(x) sum(x==minor.allele) ) )
    tb.new <- tb[ord.tb]
    
    maf <- sum( tb.new[2] + 2*tb.new[3] ) / sum( 2*tb.new )
  }
  maf
}
