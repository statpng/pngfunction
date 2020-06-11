init <- function(){
  # Variational parameters
  ## s2_jk
  ## mu_jk
  mu <- array(0, c(p, q) )
  s2 <- array(1, c(p, q) )
  
  # Model parameters
  ## theta_inv
  ## a_k
  ## alpha_jk
  ## sigb2_k
  Theta <- matrix(0, q, q)
  diag(Theta) <- diag(2/cov(y))
  alpha <- array(0.2, c(p, q) )
  a <- apply(alpha, 2, mean)
  sigb2 <- diag(cov(y)/2)
  # sigb2 <- rep(1, q)
  beta <- alpha * mu
  
  parameters <- list(
    mu = mu,
    s2 = s2,
    Theta = Theta,
    alpha = alpha,
    a = a,
    sigb2 = sigb2,
    beta = beta
  )
  parameters
}



emstep <- function(x, y, parameters){
  # E-step ------------------------------------------------------------------
  mu <- parameters$mu
  s2 <- parameters$s2
  Theta <- parameters$Theta
  alpha <- parameters$alpha
  a <- parameters$a
  sigb2 <- parameters$sigb2
  beta <- parameters$beta
  
  xx <- t(x) %*% x
  ytilde <- x %*% (alpha * mu)
  for( j in 1:p ){
    for( k in 1:q ){
      beta <- alpha * mu
      y_j <- ytilde - kronecker(beta[j,,drop=F], x[,j,drop=F])
      s2[j,k] <- 1 / ( Theta[k,k] * xx[j,j] + 1/sigb2[k] )
      mu[j,k] <- (sum( sapply(1:q, function(tt) Theta[tt,k] * x[,j] %*% (y[,tt] - y_j[,tt] ) ) ) -
                    sum( sapply((1:q)[-k], function(tt) Theta[tt,k] * xx[j,j] * beta[j,tt] ) ) ) * s2[j,k]
      alpha[j,k] <- 1/(1+exp( -(log(a[k]/(1-a[k])) + 1/2*(mu[j,k]^2 / s2[j,k] + log(s2[j,k]/sigb2[k])) ) ))
      ytilde <- y_j + kronecker(alpha[j,,drop=F] * mu[j,,drop=F], x[,j,drop=F])
    }
  }
  # for( j in 1:p ){
  #   for( k in 1:q ){
  #     beta <- alpha * mu
  #     s2[j,k] <- 1 / ( Theta[k,k] * xx[j,j] + 1/sigb2[k] )
  #     mu[j,k] <- (sum( sapply(1:q, function(tt) Theta[tt,k] * x[,j] %*% (y[,tt] - (x[,-j] %*% beta[-j,tt]) ) ) ) - 
  #                     sum( sapply((1:q)[-k], function(tt) Theta[tt,k] * xx[j,j] * beta[j,tt] ) ) ) * s2[j,k]
  #     alpha[j,k] <- 1/(1+exp( -(log(a[k]/(1-a[k])) + 1/2*(mu[j,k]^2 / s2[j,k] + log(s2[j,k]/sigb2[k])) ) ))
  #   }
  # }
  
  
  # M-step ------------------------------------------------------------------
  
  
  E <- y - ytilde
  
  Theta_inv <- matrix(0, nrow(Theta), ncol(Theta))
  beta <- alpha * mu
  for( kk in 1:q ){
    xamu <- x %*% beta[,kk]
    # Theta_kk_1 <- t(y[,kk] - xamu) %*% (y[,kk] - xamu)
    Theta_kk_1 <- t(E[,kk]) %*% E[,kk]
    # Theta_kk_2 <- sum(xx %*% (alpha[,kk]*(mu[,kk]^2+s2[,kk]) - beta[,kk]^2))
    Theta_kk_2 <- sum( sapply(1:p, function(jj) (xx[j,j] %*% (alpha[jj,kk]*(mu[jj,kk]^2+s2[jj,kk]) - beta[jj,kk]^2)) ) )
    Theta_inv[kk,kk] <- 1/n*( Theta_kk_1 + Theta_kk_2 )
    if( kk <= 4 ){
      for( tt in (kk+1):q){
        Theta_inv[kk,tt] <- 1/n*( t(E[,kk]) %*% E[,tt] )
        Theta_inv[tt,kk] <- Theta_inv[kk,tt]
      }
    }
  }
  
  Theta <- solve(Theta_inv)
  
  # cat("l = ", l, "\n")
  for( k in 1:q ){
    a[k] <- mean(alpha[,k])
    sigb2[k] <- sum(alpha[,k]*(mu[,k]^2 + s2[,k] ))/sum(alpha[,k])
  }

  
  parameters <- list(
    mu = mu,
    s2 = s2,
    Theta = Theta,
    alpha = alpha,
    a = a,
    sigb2 = sigb2,
    beta = beta
  )
  parameters
}



show.paramerters <- function(parameters){
  
  print("mu : j * l")
  print(parameters$mu[1:5,])
  
  print("s2 : j * l")
  print(parameters$s2[1:5,])
  
  print("Theta")
  print(parameters$Theta)
  
  print("alpha")
  print(parameters$alpha[1:20,])
  
  print("a")
  print(parameters$a)
  
  print("sigb2")
  print(parameters$sigb2)
  
}





ELBO <- function(x, y, L, parameters){
  gr <- parameters$gr
  PI <- parameters$PI
  eta <- parameters$eta
  eta0 <- parameters$eta0
  r <- parameters$r
  rho <- parameters$rho
  mu <- parameters$mu
  s2 <- parameters$s2
  Theta <- parameters$Theta
  alpha <- parameters$alpha
  a <- parameters$a
  sigb2 <- parameters$sigb2
  beta <- parameters$beta
  
  # L <- 3
  
  n <- nrow(y)
  p <- ncol(x)
  q <- ncol(y)
  
  elbo <- 0
  term1 <- sum( sapply( 1:q, function(kk) 1/2*alpha[,kk]*( 1 + log(s2[,kk]/sigb2[kk]) - (mu[,kk]^2 + s2[,kk])/(sigb2[kk]) ) ) )
  term2 <- - sum( sapply(1:q, function(kk) alpha[,kk] * log(alpha[,kk]/a[kk]) + (1-alpha[,kk]) * log((1-alpha[,kk])/(1-a[k])) ) )
  term3 <- n/2 * log(det(Theta))
  term5 <- - 1/2 * sum( Theta * t(y - x %*% (alpha * mu) ) %*% 
                        (y - x %*% (alpha * mu) ) )
  term6 <- - 1/2 * sum( diag(Theta) * 
                        apply( 
                          sapply(1:p, function(jj) 
                            sum(x[,jj] * x[,jj]) * 
                              ( alpha[jj,] * (mu[jj,]^2 + s2[jj,]) - alpha[jj,]^2 * mu[jj,]^2 )  
                          ), 1, sum ) )
  
  elbo <- elbo + (term1 + term2 + term3 + term5 + term6)

  elbo
}







# Simulating Dataset
varcov_rho <- function(p, rho){
  out <- matrix(rho, p, p)
  diag(out) <- 1
  return(out)
}

png.varcov <- function(p, rho=0, type=NULL, Omega = 0.001, PropOfNeg = 0.25){
  # Depends : varcov_rho
  
  if( ! type %in% c("arcov", "random_sign", "all.equal") ) stop("type should be one among arcov, random_sign, and all.equal")
  
  if( is.null(type) ) stop("'type' should be entered")
  
  if( type == "arcov" ){
    out <- outer(1:p, 1:p, function(x,y) rho^abs(x-y) ) 
    return(out)
  } 
  
  if( type == "random_sign" ){
    
    if( PropOfNeg<0 | PropOfNeg>0.5 ) stop("PropOfNeg must be in [0,0.5].");
    
    e.rho <- replicate( p*(p-1)/2, runif(1,0,1)*(-1)^rbinom(1,1,PropOfNeg) )
    e.varcov <- matrix(1, p, p)
    e.varcov[upper.tri(e.varcov)] <- e.rho
    e.varcov[lower.tri(e.varcov)] <- t(e.varcov)[lower.tri(e.varcov)]
    
    if( isSymmetric(e.varcov) ){
      x <- (e.varcov %*% e.varcov) + diag(Omega, p)
      e.varcov2 <- (x/sqrt(diag(x)%*%t( diag(x) )))
      return(e.varcov2)
    } else {
      stop("Error")
    }
  }
  
  if ( type == "all.equal" ) {
    # if( rho < 0 ) stop("rho have to be equal or greater than 0")
    e.varcov2 <- varcov_rho(p=p, rho=rho)
    return(e.varcov2)
  }
  
  # e.varcov2 %>% .[upper.tri(.)] %>% hist(main="Error-term correlation", xlab=expression(rho,"e"))
}


png.snp <- function(n, p, rho=0.95, min.maf = 0.05, maf.type = c("unif", "beta", "chisq"), SNP=TRUE, ...){
  
  X <- mnormt::rmnorm( n, varcov=png.varcov(p=p, rho=rho, type="arcov"), ...)
  if( !SNP ) return(X)
  Y <- mnormt::rmnorm( n, varcov=png.varcov(p=p, rho=rho, type="arcov"), ...)
  # Y <- mnormt::rmnorm( n, mean=mu, varcov=png.varcov(p=(p/20), rho=rho, type="arcov") ) ) )
  
  if( SNP ){
    if( maf.type == "unif" ){
      MAF <- runif(p, min.maf, 0.5)
      # ref
      # Dai, M., Ming, J., Cai, M., Liu, J., Yang, C., Wan, X., & Xu, Z. (2017). IGESS: a statistical approach to integrating individual-level genotype data and summary statistics in genome-wide association studies. Bioinformatics, 33(18), 2882-2889.
      # Yang, Y., Shi, X., Jiao, Y., Huang, J., Chen, M., Zhou, X., ... & Liu, J. (2020). CoMM-S2: a collaborative mixed model using summary statistics in transcriptome-wide association studies. Bioinformatics, 36(7), 2009-2016.
      # Jiang, W., & Yu, W. (2017). Controlling the joint local false discovery rate is more powerful than meta-analysis methods in joint analysis of summary statistics from multiple genome-wide association studies. Bioinformatics, 33(4), 500-507.
    } else if ( maf.type == "beta" ){
      MAF <- rbeta(p, 0.14, 0.73) %>% {ifelse(.<0.5, ., 1-.)*(1-min.maf*2) + min.maf }
      # ref
      # Ionita-Laza, I., Lange, C., & Laird, N. M. (2009). Estimating the number of unseen variants in the human genome. Proceedings of the National Academy of Sciences, 106(13), 5008-5013.
    } else if ( maf.type == "chisq" ){
      MAF <- truncnorm::rtruncnorm(p, a=-0.65, b=0.65, mean=0, sd=5)^2
      # MAF <- truncnorm::rtruncnorm(p, a=sqrt(min.maf), b=0.65, mean=0, sd=5)^2
      # No reference
    }
    
    
    for(j in seq_len(p)){    
      perct_X <- rank( (X[,j]) )/length(X[,j])
      perct_Y <- rank( (Y[,j]) )/length(Y[,j])
      
      X[perct_X <= MAF[j], j] <- 1
      X[perct_X >  MAF[j], j] <- 0
      Y[perct_Y <= MAF[j], j] <- 1
      Y[perct_Y >  MAF[j], j] <- 0
    }
  }
  Data <- rbind( X + Y )
  
  return(list(snp=Data, MAF=MAF));
}
###############








library(dplyr)

# Variational Bayesian Variable Selection for finite mixture model of multivariate linear regressions with correlated outcomes
set.seed(1)
i1 <- 2
i2 <- 1

q <- 5
g <- c(1, 3, 5)[i1]
y.rho <- c(0.2, 0.5, 0.8)[i2]


n=100; p=500;
snp.rho = 0.8;  maf.min = 0.05
e.varcov <- matrix(y.rho, q, q)
diag(e.varcov) <- 1

es <- rnorm(1,0,1)
beta <- replicate( q, rep(0, p) )
beta[1:10,1:g] <- es
true <- which( beta != 0, arr.ind=TRUE )

library(mnormt)
SNP <- png.snp(n=n, p=p, rho=snp.rho, min.maf=maf.min, maf.type="unif", SNP=TRUE)

X <- SNP$snp
error <- rmnorm(n=n, varcov=e.varcov)
Z <- X %*% beta + error
Y <- as.data.frame(Z)

x <- scale(X)
y <- scale(Y)

for( REP in 1:500 ){
  cat("---------------------------------------------------------\n")
  cat("REP = ", REP, "\n")
  if( REP == 1 ){
    parameters_list <- NULL
    parameters2 <- init()
  } 
  show.paramerters(parameters2)
  
  parameters_list[[REP]] <- parameters2
  # parameters <- parameters_list[[7]]
  parameters2 <- emstep(x, y, parameters2)
  
  cat("ELBO = ", ELBO(x=x, y=y, L=L, parameters=parameters), "\n")
}
