flag_mean <- function(X.list, type="svd"){
  
  rank <- function(x){
    sum( svd(x)$d > 1e-15 )
  }
  
  X.list.orth <- lapply(X.list, function(X) svd(X)$u[,1:rank(X)])
  # X.list.orth <- lapply(X.list, function(x) qr.Q(qr(x))[,1:rank(X))
  
  if( type == "eigen" ){
    
    A <- Reduce("+", lapply( X.list.orth, function(x) tcrossprod( x ) ))
    # r <- sum(svd( X )$d>0)
    eigen.A <- eigen(A)
    U <- eigen.A$vectors
    L <- eigen.A$values
    
  } else if( type == "svd" ){
    
    X.orth <- do.call("cbind", X.list.orth)
    nrank <- rank( X.orth )
    U <- with(svd(X.orth), u[,1:nrank])
    # if( dim(X.orth)[1] >= dim(X.orth)[2] ){
    #   U <- with(svd(X.orth), u[,d>1e-15])
    # } else {
    #   U <- with(svd(t(X.orth)), v[,d>1e-15])
    # }
    
  }
  
  
  U
  
}
