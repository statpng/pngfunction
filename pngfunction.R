library(mnormt)
library(glmnet)
library(dplyr)


####################################################
ArrayToLong <- function(df){

dm <- dim(df)
expd <- expand.grid( lapply(dm, function(x) 1:x ) )

out <- NULL
for( i in 1:nrow(expd) ){
	out <- rbind( out, df[ as.matrix( expd[i,], nrow=1 ) ] )
}
return( data.frame( expd, value=out ) )
}
####################################################

####################################################
btob <- function(x, name) paste0( deparse(substitute(name)),"=",x[1],"to",x[length(x)] )
####################################################

####################################################
sp.glmnet <- function(x, y, psub=0.5, K=100, seq.alpha=NULL, n.lambda=NULL, family="gaussian", type.mgaussian=NULL, ...){
  if( NROW(y) != nrow(x) ) stop("x and y should be equal length of row")
  if( NCOL(y)>1 & (family!="mgaussian") ) stop("The family should be 'mgaussian'")
  if( NCOL(y)==1 & (family=="mgaussian") ) stop("The family should not be 'mgaussian'")
  if( family=="mgaussian" & is.null(type.mgaussian) ) stop("Type.mgaussian should be entered.")
    
  if(is.null(seq.alpha)) seq.alpha <- 1:9*0.1
  if(is.null(n.lambda)) n.lambda <- 10
  
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  n <- nrow(x);
  p <- ncol(x);
  nsub <- n*psub;
  
  vector.lambda <- NULL
  
  for( j in 1:length(seq.alpha) ){
    for( i in 1:10 ){
      wsub <- sample(n, nsub)
      xsub <- x[wsub,]
      ysub <- y[wsub,]
      fitsub <- glmnet(x=xsub, y=ysub, alpha=seq.alpha[j], family=family, ... )
      vector.lambda <- c( vector.lambda, fitsub$lambda )
    }
  }
  
  lambda.min <- median(vector.lambda)
  lambda.max <- max(vector.lambda)
  seq.lambda <- seq(lambda.min, lambda.max, length.out=n.lambda);
  
  ncol.y <- length( append(fitsub$beta, NULL) )

  out <- array(0, c(ncol(x), length(seq.lambda), length(seq.alpha)) );    
for( j in 1:length(seq.alpha) ){
    for( i in 1:K ){
        set.seed(2018*i)
        wsub <- sample(n, nsub);
        xsub <- x[wsub,];
        ysub <- y[wsub,];  
            
        out.h <- array(0, c(ncol(x), length(seq.lambda), ncol.y) )
        for( h in seq_len(ncol.y) ){
            glmnet.fit <- glmnet(x=xsub, y=ysub, alpha=seq.alpha[j], lambda=seq.lambda, family=family, ... )
            out.h[,,h] <- as.numeric( append( NULL, glmnet.fit$beta )[[h]]!=0 ) ;
        }
        if(is.null(type.mgaussian) ){
            out[,,j] <- out[,,j] + apply(out.h, c(1,2), function(x) ifelse(all(x==1), 1, 0) )
        } else if(type.mgaussian == "any"){
            out[,,j] <- out[,,j] + apply(out.h, c(1,2), function(x) ifelse(any(x==1), 1, 0) )
        } else if(type.mgaussian == "all"){
            out[,,j] <- out[,,j] + apply(out.h, c(1,2), function(x) ifelse(all(x==1), 1, 0) )
        } 
  }
}
  
  return(out, seq.lambda, seq.alpha );
}
####################################################

####################################################
tpr.top <- function(sp, top, true){
    sp.sort <- sort( sp, decreasing=TRUE )
    thr <- sp.sort[top]
    sg <- sum( which(sp>thr) %in% true )
    ng <- length( which(sp>thr) )
    se <- sum( which(sp==thr) %in% true)
    ne <- length( which(sp==thr) )
    (sg+se*(top-ng)/ne)/length(true)
}
####################################################

####################################################
# print.sp <- function(out.sp, true, top=80, alpha=TRUE, lambda=TRUE, K=K){
#     
#     sp <- apply( out.sp[,lambda,alpha,drop=FALSE], 1, max )/K
#     list(
#         TPR = tpr.top(sp, top, true),
#         sp.top = data.frame(order=order( sp, decreasing=TRUE), 
#                             sp=sort( sp, decreasing=TRUE))[1:top, ],
#         nonzero = out.sp[order( sp, decreasing=TRUE )[1:top],lambda,alpha] ,
#         true = out.sp[true,lambda,alpha]
#     )
# }
####################################################

####################################################
png.GenSNP <- function(n, p, rho, threads=threads, display_progress=FALSE ){
  
  X <- do.call("cbind", lapply( 1:20, function(x) mnormt::rmnorm( n, varcov=ARCOV_C(p=(p/20), rho=rho, threads=1, display_progress=display_progress) ) ) )
  Y <- do.call("cbind", lapply( 1:20, function(x) mnormt::rmnorm( n, varcov=ARCOV_C(p=(p/20), rho=rho, threads=1, display_progress=display_progress) ) ) )
  
  MAF <- truncnorm::rtruncnorm(p, a=-0.65, b=0.65, mean=0, sd=5)^2
  
  for(j in seq_len(p)){
    
    perct_X <- rank( (X[,j]) )/length(X[,j])
    perct_Y <- rank( (Y[,j]) )/length(Y[,j])
    
    X[perct_X <= MAF[j], j] <- 1
    X[perct_X >  MAF[j], j] <- 0
    Y[perct_Y <= MAF[j], j] <- 1
    Y[perct_Y >  MAF[j], j] <- 0
    
  }
  
  Data <- rbind( X + Y )
  
  return(list(snp=Data, MAF=MAF));
}
####################################################



# ภฺมุ
####################################################
th <- function(q, alpha, p) q^2/(2*alpha*p)+1/2
####################################################

####################################################
gen.net<-function(beta,cvm,n=100,p=1000,g=100,s0=1) {
  k<-p/g
  bt<-matrix(beta,g,k)
  for (i in 1:k) {
    cx0<-rmnorm(n,bt[,i],cvm)
    if (i==1) cx<-cx0
    else cx<-cbind(cx,cx0)
  }
  for (j in 1:k) {
    tx0<-rmnorm(n*s0,rep(0,g),cvm)
    if (j==1) tx<-tx0
    else tx<-cbind(tx,tx0)
  }
  x<-rbind(cx,tx)
  y<-c(rep(1,n),rep(0,n*s0))
  return(list(x=x,y=y))
}
####################################################

####################################################
binom.glmnet.sp <- function(x, y, alpha, lambda=NULL, K=100, psub=0.5, nlamb=10) {
  w <- which(!is.na(y))
  x <- x[w,]
  y <- y[w]
  wc <- which(y==1)
  wt <- which(y==0)
  nc <- floor(length(wc) * psub)
  nt <- floor(length(wt) * psub)
  if(is.null(lambda)) {
    lam <- NULL
    for (i in 1:10) {
      for (a in 1:length(alpha)) {
        ss <- c(sample(wc, nc), sample(wt, nt))
        g0 <- glmnet(x[ss,], y[ss], alpha=alpha[a], family="binomial")
        lam <- c(lam, g0$lambda)
      }
    }
    lambda <- seq(median(lam), max(lam), length.out=nlamb)
    lambda <- sort(lambda, decreasing=T)
  }
  sp <- array(0, c(ncol(x),length(lambda), length(alpha)))
  N <- min(K, choose(length(wc), nc) * choose(length(wt), nt))
  OUT <- as.list(1:N)
  for (i in 1:N) {
    print(i)
    ss <- c(sample(wc, nc), sample(wt, nt))
    out <- as.list(1:length(alpha))
    for (a in 1:length(alpha)) {
      g <- glmnet(x[ss,], y[ss], alpha=alpha[a], lambda=lambda, family="binomial")
      sp[,,a] <- sp[,,a]+as.numeric(g$beta!=0)
#      out[[a]] <- g
    }
#   OUT[[i]] <- out
  }
#  av <- apply(sp, 1, max)/N
#  qhat <- sum(sp)/(length(alpha)*length(lambda)*N)
  return(
    list(
      sp = sp,
#      av = av,
 #     q = qhat,
      lambda = lambda,
	K=K,
	alpha=alpha,
	psub=psub
#      out = OUT
    )
  )
}
####################################################

####################################################
gaussian.glmnet.sp <- function(x, y, alpha, lambda=NULL, K=100, psub=0.5, nlamb=10) {
  w <- which(!is.na(y))
  x <- x[w,]
  y <- y[w]
  nsub <- floor(length(y) * psub)
  
  if(is.null(lambda)) {
    lam <- NULL
    for (i in 1:10) {
      for (a in 1:length(alpha)) {
        ss <- sample(length(y), nsub);
        g0 <- glmnet(x[ss,], y[ss], alpha=alpha[a], family="gaussian")
        lam <- c(lam, g0$lambda)
      }
    }
    lambda <- seq(median(lam), max(lam), length.out=nlamb)
    lambda <- sort(lambda, decreasing=T)
  }
  sp <- array(0, c(ncol(x),length(lambda), length(alpha)))
  N <- K
  for (i in 1:N) {
    print(paste0("Resampling for selection probability = ", i))
    ss <- sample(length(y), nsub)
    for (a in 1:length(alpha)) {
      g <- glmnet(x[ss,], y[ss], alpha=alpha[a], lambda=lambda, family="gaussian")
      sp[,,a] <- sp[,,a]+as.numeric(g$beta!=0)
    }
  }
  av <- apply(sp, 1, max)/N
  qhat <- sum(sp)/(length(alpha)*length(lambda)*N)
  return(
    list(
      sp = sp,
      av = av,
      q = qhat,
      lambda = lambda
    )
  )
}
####################################################
