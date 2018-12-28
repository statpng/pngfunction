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
mglmnet.selprob <- function(SNP, Phenotype, idx=NULL, Alpha=NULL, nlambda=10, psub=0.5, K=100, cov=c("none","remainder","pc")){
  if( NCOL(Phenotype)>1 & is.null(idx) ) stop("Include an integer of idx")
  if( NCOL(Phenotype)==1 & !is.null(idx) ) stop("Include a data.frame or matrix of Phenotype")
  
  M <- ncol(Phenotype)
  
  n <- NROW(Phenotype);
  nsub <- n*psub;
  if(is.null(Alpha)) Alpha <- 1:9*0.1
  

  vector.lambda <- NULL
  for( i in 1:10 ){
    for( j in 1:length(Alpha) ){
      wsub <- sample(n, nsub)
      SNPsub <- SNP[wsub,]
      Phenotypesub <- Phenotype[wsub, ]
      fitsub <- mglmnet(x=SNPsub, y=Phenotypesub, idx=idx, alpha=Alpha[j], family="gaussian", cov=cov)
      vector.lambda <- c( vector.lambda, fitsub$lambda )
    }
  }
  
  lambda.min <- min(vector.lambda)
  lambda.max <- max(vector.lambda)
  
  seq.lambda <- seq(lambda.min, lambda.max, length.out=nlambda);
  out <- array(0, c(ncol(SNP), length(seq.lambda), length(Alpha)) );
  for( i in 1:K ){
    for( j in 1:length(Alpha) ){
      wsub <- sample(n, nsub);
      SNPsub <- SNP[wsub,];
      Phenotypesub <- Phenotype[wsub, ];
      mglmnet.fit <- mglmnet(x=SNPsub, y=Phenotypesub, idx=idx, family="gaussian", alpha=Alpha[j], lambda=seq.lambda, cov=cov)
      out[,,j] <- out[,,j] + as.numeric( mglmnet.fit$beta!=0 ) ;
    }
  }
  
  return(out);
  
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
  
#  X <- do.call("cbind", lapply( 1:20, function(x) mnormt::rmnorm( n, # varcov=ARCOV_C(p=(p/20), rho=rho, threads=1,  # display_progress=display_progress) ) ) )

#  Y <- do.call("cbind", lapply( 1:20, function(x) mnormt::rmnorm( n, # varcov=ARCOV_C(p=(p/20), rho=rho, threads=1, # display_progress=display_progress) ) ) )

  X <- do.call("cbind", lapply( 1:20, function(x) mnormt::rmnorm( n, varcov=ARCOV_R(p=(p/20), rho=rho, threads=1) ) ) )

  Y <- do.call("cbind", lapply( 1:20, function(x) mnormt::rmnorm( n, varcov=ARCOV_R(p=(p/20), rho=rho, threads=1) ) ) )
  
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



####################################################
GenSimulation.png <- function(n, p, M=25, snp.rho=0.6, ptrue=c(0.001, 0.005, 0.01), Beta=1, Beta.pdec=0.1, pcross=0.2,
                              y.pNeg=0.2, Omega=0, rho=0.3){
  
  # snp <- do.call( "cbind", lapply(MAF, function(x) rbinom(n, 2, x) ) )
  GenSNP <- png.GenSNP(n=n, p=p, rho=snp.rho)
  snp <- GenSNP$snp
  MAF <- GenSNP$MAF
  dimnames(snp) <- list(paste0("N",1:n), paste0("snp", 1:p))
  ##  range( cor(snp)[upper.tri(matrix( 0, ncol(snp), ncol(snp) ) )] )
  
  K <- p*ptrue
  ponly <- 1-pcross
  Swhole <- 1:floor(pcross*K)
  Sonly <- outer( 100*seq_len(M), (seq_len( ceiling(ponly*K)) ), "+" )
  wbeta <- as.matrix( cbind(matrix(replicate(nrow(Sonly), Swhole), nrow=nrow(Sonly), byrow=TRUE), Sonly)) 
  NameWhole <- paste0(Swhole, "(Whole)")
  if( pcross==1 ){
    NameOnly <- NA
  } else{ 
    NameOnly <- paste0(length(Swhole)+seq_len(ncol(Sonly)), "(Only)")
  }
  Colname_wbeta <- as.character( na.omit( c( NameWhole, NameOnly ) ) )
  dimnames(wbeta) <- list(paste0("Phenotype.",1:M), Colname_wbeta )
  beta <- matrix(0, nrow=p, ncol=M)
  for( m in 1:M ){
    beta[Swhole, m]    <- Beta*(--1)^(m+1)
    if( pcross != 1 ){
      beta[Sonly[m,], m] <- (--1)^(m)*(Beta-Beta.pdec*(sapply(Sonly[m,], digits, which=1)-1))
    }
  }
  
  ##  beta[1:10,1:m]
  ##  beta[101:110, 1:2]
  ##  beta[201:210, 1:2]
  ##  hist( MAF[as.numeric(wbeta)] )
  
  e.varcov <- GenVarcov(M, PropOfNeg = y.pNeg, Omega=Omega, rho=rho)
  e.varcov[upper.tri(e.varcov)] %>% sign %>% table
  y <- snp %*% beta + mnormt::rmnorm(n,varcov=e.varcov)
#  corrplot::corrplot(cor(y), tl.pos = "n")
  ValuesOfArguments = list(n=n, p=p, M=M, snp.rho=snp.rho, ptrue=ptrue, Beta=Beta, pcross=pcross,
                           y.pNeg=y.pNeg, Omega=Omega, rho.y = rho)
  return(list(snp=snp, y=y, MAF=MAF, beta=beta, wbeta=wbeta, e.varcov=e.varcov, args = ValuesOfArguments));
}
####################################################

####################################################
mglmnet <- function(x=x, y=y, idx=NULL, cov=c("none", "remainder", "pc"), ...){
  yidx <- y[,idx]
  if(cov == "none"){
    return( glmnet(x=x, y=yidx, ...) )
  } else{
    if( NCOL(y)<2 ) stop("y should be a data.frame or matrix")
    if(is.null(idx)) stop("idx must be needed");
    yidx <- y[,idx]
    Remainder <- y[,-idx]
    if(cov == "remainder") {
      PF <- c( rep( 0, NCOL(Remainder) ), rep( 1, ncol(x) ) )
      out <- glmnet(x=cbind(Remainder, x), y=yidx, penalty.factor=PF, ...)
      out$dim[1] <- ncol(x)
      out$beta <- out$beta[-c(1:NCOL(Remainder)),]
      
      
      return( out )
    } else if(cov == "pc"){
      PC <- prcomp(Remainder, scale.=TRUE)
      PCcov <- PC$x[,which( summary( PC )$sdev > 1), drop=FALSE]
      PF <- c( rep( 0, NCOL(PCcov) ), rep( 1, ncol(x) ) )
      out <- glmnet(x=cbind(PCcov, x), y=yidx, penalty.factor=PF, ...)
      out$dim[1] <- ncol(x)
      out$beta <- out$beta[-c(1:NCOL(PCcov)),]
      return( out )
    }
  }
}
####################################################

####################################################
Method.lm <- function(SNP, Phenotype){
pvalue.lm <- sapply( 1:ncol(SNP), function(j) {
  fit_coef <- summary( lm( Phenotype ~ SNP[,j] ) )$coef
  ifelse( nrow(fit_coef)==2, fit_coef[2,4], 9999 )
} )
return(-log10(pvalue.lm))
}
####################################################

####################################################
mcc <- function(true, pred, p){
  data.frame(
    Actual = ifelse((1:p)%in% true, "TRUE", "FALSE"),
    Condition = ifelse((1:p)%in% pred, "+", "-") ) %>% table -> TABLE
  TP <- TABLE[2,2]
  FP <- TABLE[1,2]
  FN <- TABLE[2,1]
  TN <- TABLE[1,1]
  return( (TP*TN-FP*FN)/(sqrt(TP+FP)*sqrt(TP+FN)*sqrt(TN+FP)*sqrt(TN+FN)) );
}
####################################################

####################################################
glmnet_ycov <- function( SNP, Phenotype, idxY, Alpha, method ){
  
  Y <- Phenotype[,idxY]
  subY <- Phenotype[,-idxY]
  
  if( method == "glmnet" ){
    cvglmnet.fit <- cv.glmnet(SNP, Y, family="gaussian", alpha=Alpha, type.measure="mse" )
    glmnet.fit <- glmnet(SNP, Y, family="gaussian", alpha=Alpha, lambda=cvglmnet.fit$lambda.1se )
    return( which( glmnet.fit$beta != 0 ) );
    
  } else if( method == "ycov" ){
    X <- cbind(subY, SNP)
    cvglmnet.fit <- cv.glmnet(X, Y, family="gaussian", alpha=Alpha, type.measure="mse", penalty.factor = c( rep(0,ncol(subY)), rep(1,ncol(X)-ncol(subY)) ) ) 
    glmnet.fit <- glmnet(X, Y, family="gaussian", alpha=Alpha, lambda=cvglmnet.fit$lambda.1se, penalty.factor = c( rep(0,ncol(subY)), rep(1,ncol(X)-ncol(subY)) ) ) 
    return( which( glmnet.fit$beta != 0 ) );
    
  } else if( method == "pccov"){
    X <- cbind(SNP, Z <- prcomp(subY, scale=TRUE)$x[,1])
    cvglmnet.fit <- cv.glmnet(X, Y, family="gaussian", alpha=Alpha, type.measure="mse", penalty.factor=c(rep(0,1),rep(1,ncol(SNP))) ) 
    glmnet.fit <- glmnet(X, Y, family="gaussian", alpha=Alpha, lambda=cvglmnet.fit$lambda.1se, penalty.factor=c(rep(0,1),rep(1,ncol(SNP))) )
    return( which( glmnet.fit$beta != 0) );
    
  } else if( method == "lm" ){
    pvalue.lm <- sapply( 1:ncol(SNP), function(j) summary( lm( Y ~ SNP[,j] ) )$coef[2,4] )
    return( pvalue.lm );
  }
}
####################################################

####################################################
ARCOV_R <- function(p, rho, threads=1) outer(1:p, 1:p, function(x,y) rho^abs(x-y) )

minmax <- function(x, lag){ 
  out <- ( x-min(x) )/(max(x)-min(x) )
  out 
}
####################################################

####################################################
SampleSign <- function(Matrix, PropOfNeg=PropOfNeg){
  Matrix[ lower.tri(Matrix) ] <- sapply( Matrix[ lower.tri(Matrix) ], function(x) (-1)^rbinom(1,1,PropOfNeg)*x )
  Matrix[ upper.tri(Matrix) ] <- t(Matrix)[ upper.tri(Matrix) ]
  Matrix
}
digits <- function(x, which=1e1){
  as.numeric( unlist(strsplit(as.character(x),""))[nchar(x)+1-log10(10*which)] )
}
####################################################

####################################################
Varcov_rho <- function(p, rho){
  out <- matrix(rho, p, p)
  diag(out) <- 1
  return(out)
}
####################################################

####################################################
GenVarcov <- function(m, Omega = 0, PropOfNeg = 0.25, rho=0){
  if( PropOfNeg<0 | PropOfNeg>0.5 ) stop("PropOfNeg must be in [0,0.5].");
  
  if( rho == 0 ){
  e.rho <- replicate( m*(m-1)/2, runif(1,0,1) )
  e.varcov <- matrix(1, m, m)
  e.varcov[upper.tri(e.varcov)] <- e.rho
  e.varcov[lower.tri(e.varcov)] <- t(e.varcov)[lower.tri(e.varcov)]
  e.varcov <- SampleSign(e.varcov, PropOfNeg = PropOfNeg)
  isSymmetric(e.varcov)
  
  x <- (e.varcov %*% e.varcov) + diag(Omega, m)
  e.varcov2 <- (x/sqrt(diag(x)%*%t( diag(x) )))
  
  } else if ( rho>0 ) {
    e.varcov2 <- Varcov_rho(p=m, rho=rho)
  }
  
  # e.varcov2 %>% .[upper.tri(.)] %>% hist(main="Error-term correlation", xlab=expression(rho,"e"))
  
  invisible(e.varcov2)
}
####################################################

