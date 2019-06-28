library(mnormt)
library(glmnet)
library(dplyr)


png.snpimpute <- function(xx){
  x.na <- xx[is.na(xx)]
  x.value <- xx[!is.na(xx)]
  
  tb <- table( factor( x.value, levels=0:2 ) )
  
  y <- sum(x.value)
  n <- 2*sum(tb)
  # maf <- (tb[2]+2*tb[3])/(2*sum(tb))
  # pi(p) ~ beta(2, 2)
  # L(y|p) ~ B(p)
  # pi(p|y) \prop pi(p) * L(y|p)
  #         ~ Beta(y+2, n+2-y)
  maf <- rbeta(n=length(x.na), shape1=y+2, n+2-y )/2
  impute.value <- rbinom(n=length(x.na), size = 2, prob = maf)
  xx[is.na(xx)] <- impute.value
  xx
}



sparse.dim <- function(sparsematrix){
  rows <- sparsematrix@i + 1
  cols <- findInterval(seq(sparsematrix@x)-1,sparsematrix@p[-1])+1
  list( rows=rows, cols=cols )
}

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
digits <- function(x, which=1e1){
  as.numeric( unlist(strsplit(as.character(x),""))[nchar(x)+1-log10(10*which)] )
}
####################################################

####################################################
minmax <- function(x, lag){ 
  out <- ( x-min(x) )/(max(x)-min(x) )
  out 
}
####################################################


####################################################
# Replace this by plyr::adplyr( array, c(1,2,3) )
# ArrayToLong <- function(df){
# 
# dm <- dim(df)
# expd <- expand.grid( lapply(dm, function(x) 1:x ) )
# 
# out <- NULL
# for( i in 1:nrow(expd) ){
# 	out <- rbind( out, df[ as.matrix( expd[i,], nrow=1 ) ] )
# }
# return( data.frame( expd, value=out ) )
# }
####################################################

####################################################
btob <- function(x, name) paste0( deparse(substitute(name)),"=",x[1],"to",x[length(x)] )
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
th <- function(q, alpha, p) q^2/(2*alpha*p)+1/2
####################################################

####################################################
varcov_rho <- function(p, rho){
  out <- matrix(rho, p, p)
  diag(out) <- 1
  return(out)
}
####################################################


####################################################
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
####################################################

####################################################
png.snp <- function(n, p, rho, min.maf = 0.05){
  
  X <- do.call("cbind", lapply( 1:20, function(x) mnormt::rmnorm( n, varcov=png.varcov(p=(p/20), rho=rho, type="arcov") ) ) )
  Y <- do.call("cbind", lapply( 1:20, function(x) mnormt::rmnorm( n, varcov=png.varcov(p=(p/20), rho=rho, type="arcov") ) ) )
  
  MAF <- truncnorm::rtruncnorm(p, a=sqrt(min.maf), b=0.65, mean=0, sd=5)^2
  # MAF <- truncnorm::rtruncnorm(p, a=-0.65, b=0.65, mean=0, sd=5)^2
  # MAF <- rbeta(n=p, shape1=2, shape2=5)
  
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

####################################################
png.msnp <- function(n, p, M=25, snp.rho=0.6, y.rho=0.3, maf.min=0.01, 
					sp.mu, cp.mu, ncp, nsp, ncp.size, 
					cov.type="all.equal", y.p.neg=NULL, omega=NULL,
					standardization=TRUE){
	# Depends: png.snp, png.varcov, mnormt
  
 # sp.mu <- rnorm(1, 0.3, 0.05)
 # cp.mu <- rnorm(1, 0.3, 0.05)
  
  PNG.SNP <- png.snp(n=n, p=p, rho=snp.rho, min.maf=maf.min)
  if( is.matrix(y.rho) ){
	VAR <- y.rho
  } else {
	VAR <- png.varcov(M, rho=y.rho, type=cov.type)
  }
  SNP <- PNG.SNP$snp
  dimnames(SNP) <- list(paste0("N",1:n), paste0("snp", 1:p))
  MAF <- PNG.SNP$MAF
  
  nsig <- ncp+nsp
  # SNP <- do.call( "cbind", lapply(MAF, function(x) rbinom(n, 2, x) ) )
  ##  range( cor(SNP)[upper.tri(matrix( 0, ncol(SNP), ncol(SNP) ) )] )  
     
  if(nsp %% M != 0) warnings("The number of single-phenotypic variants is not proportional to M!")
  if(ncp>0)  var.cp <- (1:nsig)[  (1:nsig)%in%seq_len(ncp)]
  if(nsp>0)  var.sp <- (1:nsig)[!(1:nsig)%in%seq_len(ncp)]

beta <- replicate( M, rep(0, p) )
if(ncp>0){
	for( m in 1:ncp.size ){
		beta[ var.cp*10, m ] <- cp.mu
	}
}
if(nsp>0)  {
	var.sp.list <- tapply(var.sp, rep(1:M, each=ceiling(nsp/M))[1:nsp], list)
	for( m in 1:length(var.sp.list) ){
		beta[ var.sp.list[[m]]*10, m] <- sp.mu
	}
}


#  for( ss in 1:nsig ) {
#	if( ss %in% var.cp ){
#		for( m in 1:ncp.size ){
#			beta[ var.cp[ss]*10, m] <- cp.mu
#		}
#	} else if ( ss %in% var.sp ){
#		for( m in 1:length(var.sp) ){
#			beta[ var.sp[m]*10, ((m-1)%%M+1)] <- sp.mu
#		}
#	}
# }
  
  
  true <- which( beta!=0, arr.ind = TRUE )
if(ncp>0)  {
	true.cp <- unique( true[duplicated(true[,1]), 1, drop=T] )
} else {
	true.cp <- NULL
}
  if(nsp>0)  {
	  if( is.null(true.cp) ) {
			true.sp <- true[,1]
		} else {
			true.sp <- true[ ! true[,1] %in% true.cp, 1 ]
		}
		true.sp <- tapply(true.sp , rep(1:M, each=ceiling(nsp/M))[1:nsp], list )
  } else {
	true.sp <- NULL
  }

  
  
  Y <- SNP %*% beta + rmnorm(n=n, varcov=VAR)
  #  cor(Y)
  #  corrplot::corrplot(cor(y), tl.pos = "n")
        
  ValuesOfArguments = 
	list(
	n=n, 
	p=p, 
	M=M, 
	snp.rho=snp.rho, 
	y.rho=y.rho, 
	maf.min=maf.min, 
	sp.mu=sp.mu, 
	cp.mu=cp.mu, 
	ncp=ncp, 
	nsp=nsp,
	cov.type=cov.type
	)
  
  if(standardization) SNP <- scale(SNP)
	
	Data <- list(snp=SNP, y=Y, maf=MAF, beta=beta, true=true, true.cp=true.cp, true.sp=true.sp, args=ValuesOfArguments)
  return(Data);
}
####################################################

