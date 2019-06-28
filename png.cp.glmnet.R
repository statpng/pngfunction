png.cp.glmnet <- function(x, y, family, seq.alpha=NULL, seq.lambda=NULL, K=100, setseed, psub=0.5, ...){
  library(mnormt)
  library(glmnet)
  library(dplyr)
  
  # x=Data$snp
  # y=Data$y
  # family="mixed"
  # psub=0.5
  # seq.alpha=seqalpha
  # n.lambda=nlambda
  # setseed=1129
  
  if( NROW(y) != nrow(x) ) stop("x and y should be equal length of row")
  if(is.null(seq.alpha)) seq.alpha <- 1:9*0.1
  if(family=="mgaussian") standardize.response <- TRUE
  
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  n <- nrow(x);
  p <- ncol(x);
  nsub <- n*psub;

  
  n.alpha <- length(seq.alpha)
  n.lambda <- unique( sapply(seq.lambda, length) )
  
  beta.array <- array(NA,
                      dim = c(p, n.alpha, n.lambda, ncol(y), K), 
                      dimnames = list( paste0("", 1:p),  # paste0("v", 1:p), 
                                       paste0("",seq.alpha),   # paste0("Alpha=",seq.alpha),  
                                       paste0("", seq_len(n.lambda)), # paste0("Lambda=", seq_len(n.lambda)),
                                       paste0("", seq_len(ncol(y))), # paste0("Phenotype=", seq_len(ncol(y))),
                                       paste0("", seq_len(K)) ) ) # paste0("Rep=", seq_len(K)) ) )
  
  names(attributes(beta.array)$dimnames) <- c("Variable", "Alpha", "Lambda", "Phenotype", "Replications")
  
  for( kk in 1:K ){
    print(kk)
    start <- proc.time();
    set.seed( setseed*kk )
    wsub <- sample(n, nsub)
    xsub <- x[wsub,,drop=F];
    ysub <- y[wsub,,drop=F];

    
    for( aa in 1:length(seq.alpha) ){
      if( family=="mgaussian" ){
      mgaussian.fit <- append(NULL, glmnet(x=xsub, 
                                           y=ysub, 
                                           alpha=seq.alpha[aa], 
                                           lambda=unlist(seq.lambda[1]), 
                                           family=family, 
                                           standardize.response=TRUE, ... )$beta )
      }
      
      for( colcol in 1:ncol(y) ){
        if( family=="mgaussian" ){
          beta.array[,aa,,colcol,kk] <- as.numeric( mgaussian.fit[[colcol]] != 0 )
        } else {
          beta.array[,aa,,colcol,kk] <- as.numeric( glmnet(x=xsub, 
                                                       y=ysub[,colcol,drop=F], 
                                                       alpha=seq.alpha[aa], 
                                                       lambda=seq.lambda[[colcol]], 
                                                       family=family, ... )$beta != 0 )
        }
      }
      
    }
    end <- proc.time();
    cat("In this iteration, the elaspsed time=", (end-start)[3], "\n");
    cat("The remained time is", (end-start)[3]*(K-kk), "\n");
  }
  
  return(beta.array);
  
}



# gaussian.union <- apply( gaussian.beta.array, c(1,2), function(nonzero) ifelse(sum(nonzero) > 0, 1, 0) )
# mgaussian.union <- apply( mgaussian.beta.array, c(1,2), function(nonzero) ifelse(sum(nonzero)> 0, 1, 0) )
png.get_sp <- function(array){
  Margin.rep <- which( !names(dimnames(array)) %in% c("Replications") )
  count.array <- apply( array, Margin.rep, mean )
  Margin.tuning <- which( !names(dimnames(count.array)) %in% c("Alpha", "Lambda") )
  as.matrix(apply( count.array, Margin.tuning, max ))
}

























png.cp.sgl <- function(x, y, type, seq.alpha=NULL, seq.lambda=NULL, K=100, setseed, psub=0.5, ...){
  library(mnormt)
  library(glmnet)
  library(dplyr)
  
  # x=Data$snp+1e5
  # y=Data$y
  # type="linear"
  # psub=0.5
  # seq.alpha=seqalpha
  # seq.lambda=Lambda.sgl.y
  # setseed=1129
  
  if( NROW(y) != nrow(x) ) stop("x and y should be equal length of row")
  if(is.null(seq.alpha)) seq.alpha <- 1:9*0.1
  
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  n <- nrow(x);
  p <- ncol(x);
  q <- ncol(y);
  nsub <- n*psub;
  
  
  n.alpha <- length(seq.alpha)
  n.lambda <- unique( sapply(seq.lambda, length) )
  
  
  kronecker <- function(X, Y){
    x <- matrix(0, nrow(Y)*ncol(Y), ncol(X)*ncol(Y) )
    for( i in 1:ncol(Y)){
      n <- nrow(X)
      p <- ncol(X)
      x[ (n*(i-1)+1):(n*i), (p*(i-1)+1):(p*i) ] <- X
    }
    list(x=x, y=as.matrix(as.vector(Y)), 
         group.x=rep(1:ncol(Y), each=ncol(X)), group.y=rep(1:ncol(Y), each=nrow(Y)),
         group.xy=rep(1:ncol(X), ncol(Y)) )
  }
  
  fit.kronecker <- kronecker(x, y)
  x <- fit.kronecker$x
  y <- fit.kronecker$y
  index <- fit.kronecker$group.xy
  
  
  n.new <- nrow(x);
  p.new <- ncol(x);
  q.new <- ncol(y)
  
  
  beta.array <- array(NA,
                      dim = c(p, n.alpha, n.lambda, q, K), 
                      dimnames = list( paste0("", 1:p),  # paste0("v", 1:p), 
                                       paste0("",seq.alpha),   # paste0("Alpha=",seq.alpha),  
                                       paste0("", seq_len(n.lambda)), # paste0("Lambda=", seq_len(n.lambda)),
                                       paste0("", seq_len(q)), # paste0("Phenotype=", seq_len(ncol(y))),
                                       paste0("", seq_len(K)) ) ) # paste0("Rep=", seq_len(K)) ) )
  
  names(attributes(beta.array)$dimnames) <- c("Variable", "Alpha", "Lambda", "Phenotype", "Replications")
  
  for( kk in 1:K ){
    print(kk)
    start <- proc.time();
    set.seed( setseed*kk )
    wsub <- sample(n, nsub)
    wsub.new <- unlist( lapply( (seq_len(q)-1)*n, function(x) x + wsub ) )
    xsub <- x[wsub.new,,drop=F];
    ysub <- y[wsub.new,,drop=F];
    
    
    for( aa in 1:length(seq.alpha) ){
      for( colcol in 1:length(seq.lambda) ){
        beta.sgl <- SGL(list(x=xsub, y=ysub), 
                        index=index,
                        alpha=seq.alpha[aa], 
                        lambda=seq.lambda[[colcol]], 
                        type=type, 
                        verbose=FALSE, ...)$beta
        for( lamlam in 1:ncol(beta.sgl) ){
          beta.array[,aa,lamlam,,kk] <- do.call("rbind", tapply( as.numeric(beta.sgl[,lamlam]!= 0), index, c ) )
        }
        
      }
    }
      
    end <- proc.time();
    cat("In this iteration, the elaspsed time=", (end-start)[3], "\n");
    cat("The remained time is", (end-start)[3]*(K-kk), "\n");
  }
  
  return(beta.array);
  
}
