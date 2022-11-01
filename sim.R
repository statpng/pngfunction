png.sim_lowrank = function(n, p, q, d, r, p0, q0, pvec=NULL, rvec=NULL, rho_X=0.5, rho_E=0.5, snr=0.5, sigma=NULL, type.cov="AR"){
  
  # First-order autoregressive structure
  CorrAR <- function(p, rho){
    rho^outer(1:p, 1:p, function(x, y) abs(x-y))
  }

  # Compound symmetry structure
  CorrCS <- function(p, rho){
    Sigma <- matrix(nrow = p, ncol = p, rho)
    diag(Sigma) <- 1
    Sigma
  }

  Sigma = switch(type.cov, 
                 "AR" = CorrAR, "CS" = CorrCS
                 )
  
  
  if( is.null(pvec) ) pvec <- rep( c(floor(p/d)+1, floor(p/d)), c(p%%d,d-p%%d) )
  if( is.null(rvec) ) rvec <- rep( c(floor(r/d)+1, floor(r/d)), c(r%%d,d-r%%d) )
  
  
  p = sum(pvec)
  r = sum(rvec)
  pvec2 = rep(1:length(pvec), pvec)
  rvec2 = rep( seq_len(length(rvec)), rvec )
  
  
  A1 <- matrix(ncol = r, nrow = q0, rnorm(q0 * r))
  A0 <- matrix(ncol = r, nrow = q - q0, 0)
  A <- rbind(A1, A0)
  A <- with( svd(A), tcrossprod(u,v) )
  
  B1 <- matrix(ncol = r, nrow = p0, rnorm(p0 * r))
  B0 <- matrix(ncol = r, nrow = p - p0, 0)
  B <- rbind(B1, B0)
  
  C <- B %*% t(A)
  
  
  
  X <- MASS::mvrnorm(n, rep(0, p), Sigma(p, rho_X))
  UU <- MASS::mvrnorm(n, rep(0, q), Sigma(q, rho_E))
  
  svdC <- svd(C)
  C3 <- with( svdC, u[, r] %*% t(v[, r]) * d[r] )
  Y3 <- X %*% C3
  
  if (is.null(sigma)) {
    sigma <- sqrt(sum(as.numeric(Y3)^2)/sum(as.numeric(UU)^2)/snr)
  }
  UU <- UU * sigma
  Y <- matrix(nrow = n, ncol = q, NA)
  Y <- X %*% C + UU
  
  params = list(n=n, p=p, q=q, d=d, pvec=pvec, rvec=rvec,
                p0=p0, q0=q0,
                sigma=sigma, snr=snr, rho_X=rho_X, rho_E=rho_E)
  
  list(Y=Y, X=X, A=A, B=B, params=params)
}
