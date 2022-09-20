png.quartimax <- function(X, lr=0.1, eps=1e-10){
  
  Q <- diag(ncol(X))
  p <- nrow(X)
  
  Z <- X %*% Q
  dQ <- - Z^3
  G <- crossprod(X, dQ)
  while(TRUE){
    Qold <- Q
    
    Q <- with( svd( Q + lr*G ), tcrossprod(u, v) )
    
    Z <- X %*% Q
    dQ <- - Z^3
    # dQ <- 1/p * (Z^3 - Z %*% diag(drop(rep(1, p) %*% Z^2))/p)
    G <- crossprod(X, dQ)
    
    if( norm(Q - Qold, "F") < eps ) break
  }
  
  print(norm(Q - Qold, "F"))
  list(loadings = X %*% Q, rotation = Q)
}






out_total = NULL
for( j in 1:100 ){
	set.seed(j)
	A = matrix(rnorm(10*5),10,5)
	# A = with( svd(matrix(rnorm(20*5),10,5)), tcrossprod(u,v) );  

	out = NULL; 
	for( i in 1:100 ){ 
		set.seed(i+100)
		Q = with( svd(matrix(rnorm(5*5),5,5)), tcrossprod(u,v) );  
		out[[i]] = A %*% Q
	}

	# out2 = lapply(out[1:100], function(L) png.quartimax(L)$load )
	out2 = lapply(out[1:100], function(L) GPArotation::quartimax(L, maxit=1000, eps=1e-3)$load )
	out3 = lapply(out2, function(x) x[,order(apply(x, 2, function(y) abs(y[1])))] )
	out4 = sapply(out3, function(x) abs(x[1,1]) )
	out_total[j] = mean( diff(out4) < 1e-4 ) 
	print(out_total[j])
}

