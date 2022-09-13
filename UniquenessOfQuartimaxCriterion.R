png.quartimax <- function(X){
  
  Q <- diag(ncol(X))
  p <- nrow(X)
  
  Z <- X %*% Q
  dQ <- - Z^3
  G <- crossprod(X, dQ)
  while(TRUE){
    Qold <- Q
    
    Q <- with( svd( Q + 0.1*G ), tcrossprod(u, v) )
    
    Z <- X %*% Q
    dQ <- - Z^3
    # dQ <- 1/p * (Z^3 - Z %*% diag(drop(rep(1, p) %*% Z^2))/p)
    G <- crossprod(X, dQ)
    
    if( norm(Q - Qold, "F") < 1e-11 ) break
  }
  list(loadings = X %*% Q, rotation = Q)
}



A = with( svd(matrix(rnorm(20*5),10,5)), tcrossprod(u,v) );  

crossprod(A)

out = NULL; 
for( i in 1:10 ){ 
	Q = with( svd(matrix(rnorm(5*5),5,5)), tcrossprod(u,v) );  
	out[[i]] = A %*% Q
}


GPArotation::quartimax( A %*% Q )$Th
png.quartimax( A %*% Q )$rot

out[1:5]

lapply(out[1:5], function(L) png.quartimax(L)$load )
lapply(out[1:5], function(L) GPArotation::quartimax(L, maxit=100000, eps=5e-7)$load )

A=matrix(rnorm(6*3),6,3); out=NULL; for(i in 1:100) out[[i]] = GPArotation::Varimax( A %*% with( svd(matrix(rnorm(9),3,3)), tcrossprod(u,v) ))$load


