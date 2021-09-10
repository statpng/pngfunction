png.split <- function(nrep, nsplit = 10){
  n <- floor( nrep / nsplit )
  n1 <- nrep %% nsplit
  n2 <- nsplit - n1
  
  n.vec <- rep( c(n+1, n), c(n1, n2) )
  
  start <- end <- 0
  AA <- NULL
  for( i in 1:length(n.vec) ){
    start <- start + c(1,n.vec)[i]
    end <- end + c(1,n.vec)[i+1]
    
    AA[[i]] <- seq(start, end, by=1)
  }
  
  cat( "The cardinal number of each set =", sapply(AA, length), "\n" )
  
  AA
}
