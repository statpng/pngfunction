flag_mean <- function(X.list, nrank=NULL, type=c("svd","eigen")){
  
  # X.list <- data[[1]] %>% { lapply( 1:4, function(idx) .$X[,which(.$pvec==idx)] ) }
  
  X.list.orth <- lapply(X.list, function(x) qr.Q(qr(x)))
  # X.list.orth <- lapply(X.list, function(X) svd(X)$u)
  
  X.orth <- do.call("cbind", X.list.orth)
  
  if( type == "eigen" ){
    A <- Reduce("+", lapply( X.list.orth, function(x) tcrossprod( x ) ))
    # r <- sum(svd( X )$d>0)
    eigen.A <- eigen(A)
    U <- eigen.A$vectors
    L <- eigen.A$values
  } else {
    U <- svd(X.orth)$u
  }
  
  if(!is.null(nrank)){
    U <- U[,1:nrank]
  }
  
  U
}
