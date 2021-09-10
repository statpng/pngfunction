CTGSA.perm <- function(fit, nperm = 1000, seed = 1){
  n <- nrow(fit$params$x)
  p <- ncol(fit$params$x)
  
  out <- matrix(NA, nperm, 4)
  
  perm.params <- fit$params
  pb <- txtProgressBar(min=0, max=nperm, style=3)
  for( perm.i in 1:nperm ){
    setTxtProgressBar(pb, perm.i)
    
    if( perm.i == 1 ) start <- proc.time()
    # permutation -------------------------------------------------------------
    set.seed( seed*perm.i )
    perm.params$y <- fit$params$y[sample(n)]
    
    fit.perm <- do.call("CTGSA", perm.params)
    
    out[perm.i,] <- unlist( fit.perm$stat )
    
    if( perm.i == 1 ) cat( "The expected time left = ", (proc.time() - start)["elapsed"] * nperm / 60, " (min) \n" )
  }
  
  out.pval <- NULL
  for( k in 1:length(fit$stat) ){
    out.pval[k] <- (sum( out[,k] > fit$stat[[k]] ) + 1) / (nperm+1)
  }
  
  list(pvalue = out.pval, stat = out)
  
}
