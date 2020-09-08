# available for multinomial type and with penalty factor

# x <- rbind( mnormt::rmnorm(100, rep(0, 1000), diag(1000)), mnormt::rmnorm(100, rep(1, 1000), diag(1000)) )
# y <- data.frame(y1=rep(0:1, each=100), y2=x%*%rep(0.1, 1000)+rnorm(200), y3=rbinom(200, 5, prob=rep(c(0.4,0.6), each=100)))
# # family <- rep("gaussian", ncol(y))
# family <- c("binomial", "gaussian", "multinomial")
# type.multivariate <- "unet"
# seq.alpha <- 1*0.1
# seq.lambda <- replicate(3, seq(0.5, 20.0, length.out=50), simplify = F)
# K <- 10
# psub <- 0.5
# 
# kk <- aa <- colcol <- 1
# penalty.factor <- c(0,0,0,0, rep(1, ncol(x)-4))
# ... <- NULL

# -------------------------------------------------------------------------

png.cp_glmnet2 <-
  function(x,
           y,
           family,
           type.multivariate = c("unet", "mnet"),
           seq.alpha = NULL,
           seq.lambda = NULL,
           K = 100,
           psub = 0.5,
           penalty.factor = NULL,
           ...) {
    
    library(mnormt)
    library(glmnet)
    library(dplyr)
    
    
    
    if( ! type.multivariate %in% c("enet", "unet", "mnet", "wnet") ) stop("type.multivariate %in% ('enet', 'unet', 'mnet', 'wnet') ")
    if( length(family) != ncol(y) ) stop("The length of family should be equal to ncol(y).")
    if( is.list(seq.lambda) & length(seq.lambda) != ncol(y) ) stop("The length of seq.lambda should be equal to ncol(y).")
    if( any(family %in% c("binomial", "multinomial")) & (type.multivariate %in% c("mnet")) ) stop("type.multivariate cannot be 'mnet' in family of binomial and multinomial")
    if( any(family %in% c("mgaussian")) ) stop("'mgaussian' should be used in type.multivariate")
    
    if( NROW(y) != nrow(x) ) stop("x and y should be equal length of row")
    if( is.null(seq.alpha)) seq.alpha <- 1:9*0.1
    if( type.multivariate=="mnet") standardize.response <- TRUE
    
    x <- as.matrix(x)
    y <- as.data.frame(y, stringsAsFactors = FALSE)
    
    n.pf <- sum(penalty.factor==0)
    if( n.pf > 0 ){
      wh.var <- (1:ncol(x))[-(1:n.pf)]
    } else {
      wh.var <- (1:ncol(x))
    }
    
    n <- nrow(x);
    p <- ncol(x)-n.pf;
    nsub <- n*psub;
    
    n.alpha <- length(seq.alpha)
    n.lambda <- unique( sapply(seq.lambda, length) )
    
    
    # IF sp -------------------------------------------------------------------
    beta.array.enet <- array(0, dim = c(p, n.alpha, n.lambda, ncol(y)), 
                             dimnames = list( paste0("", 1:p),  # paste0("v", 1:p), 
                                              paste0("",seq.alpha),   # paste0("Alpha=",seq.alpha),  
                                              paste0("", seq_len(n.lambda)), # paste0("Lambda=", seq_len(n.lambda)),
                                              paste0("", seq_len(ncol(y))) ) )
    
    names(attributes(beta.array.enet)$dimnames) <- c("Variable", "Alpha", "Lambda", "Phenotype")
    
    
    beta.array.unet <- 
      beta.array.unet_sum <- beta.array.unet_sum_max <- 
      # beta.array.unet_wsum <- beta.array.unet_wsum_max <- 
      beta.array.mnet <- beta.array.enet
    
    beta.array.unet <- slam::as.simple_sparse_array(beta.array.unet)
    beta.array.unet_sum <- slam::as.simple_sparse_array(beta.array.unet_sum)
    beta.array.unet_sum_max <- slam::as.simple_sparse_array(beta.array.unet_sum_max)
    # beta.array.unet_wsum <- slam::as.simple_sparse_array(beta.array.unet_wsum)
    # beta.array.unet_wsum_max <- slam::as.simple_sparse_array(beta.array.unet_wsum_max)
    beta.array.mnet <- slam::as.simple_sparse_array(beta.array.mnet)
    beta.array.enet <- slam::as.simple_sparse_array(beta.array.enet)
    
    
    for( kk in 1:K ){
      print(kk)
      start <- proc.time();
      wsub <- sample(n, nsub)
      xsub <- x[wsub,,drop=F];
      ysub <- y[wsub,,drop=F];
      
      for( aa in 1:length(seq.alpha) ){
        
        if( type.multivariate=="mnet" ){
          mgaussian.fit <- glmnet(x=xsub,
                                   y=ysub,
                                   alpha=seq.alpha[aa],
                                   lambda=unlist(seq.lambda[1]), 
                                   family="mgaussian", 
                                   standardize.response=TRUE, 
                                   penalty.factor = penalty.factor, ... )$beta[[1]][wh.var,]
          for( colcol in 1:ncol(y) ){
            beta.array.mnet[,aa,,colcol] <- as.array(beta.array.mnet[,aa,,colcol]) + as.numeric( mgaussian.fit != 0 )/K
          }
        }
        
        if ( type.multivariate=="unet" ) {
          # begin colcol ------------------------------------------------------------
          unet.matrix <- NULL
          for( colcol in 1:ncol(y) ){
            FAMILY <- family[colcol]
            
            if(FAMILY=="multinomial"){
              tmp.y <- model.matrix( ~.-1, data=as.data.frame(as.factor(ysub[,colcol,drop=T])%>% droplevels()) )
              type.multinomial <- "grouped"
            } else {
              tmp.y <- ysub[,colcol,drop=T]
              type.multinomial <- NULL
            }
            
            fit <- glmnet( x=xsub, y=tmp.y, 
                           alpha=seq.alpha[aa], lambda=seq.lambda[[colcol]], 
                           family=FAMILY, type.multinomial = type.multinomial, 
                           penalty.factor = penalty.factor, ...)
            
            fitted.beta <- switch(as.character(is.list(fit$beta)), "TRUE"=fit$beta[[1]], "FALSE"=fit$beta)
            fitted.beta <- fitted.beta[wh.var,]
            
            unet.matrix <- cbind( unet.matrix, as.numeric( fitted.beta != 0 ) )
            
            beta.array.enet[,aa,,colcol] <- as.array(beta.array.enet[,aa,,colcol]) + as.numeric( fitted.beta != 0 )/K
            
          }
          
          # end colcol --------------------------------------------------------------
          
          unet.sum.vec <- apply( unet.matrix, 2, sum )
          
          uwnet.matrix <- NULL
          for( colcol in 1:ncol(y) ){
            w.vec <- (prod(unet.sum.vec)/sum(unet.sum.vec)/unet.sum.vec) %>% {./sum(.)}
            uwnet.matrix <- cbind( uwnet.matrix, unet.matrix[,colcol]*w.vec[colcol] )
          }
          
          for( colcol in 1:ncol(y) ){
            beta.array.unet[,aa,,colcol] <- as.array(beta.array.unet[,aa,,colcol]) + apply( unet.matrix, 1, function(x) ifelse( sum(x)>0, 1, 0 ) )/K
            
            beta.array.unet_sum[,aa,,colcol] <- as.array(beta.array.unet_sum[,aa,,colcol]) + apply( unet.matrix, 1, function(x) sum(x) )/K
            # beta.array.unet_wsum[,aa,,colcol] <- as.array(beta.array.unet_wsum[,aa,,colcol]) + apply( uwnet.matrix, 1, function(x) sum(x) )/K
          }
          
          
          beta.array.unet_sum_max[,aa,,] <- beta.array.unet_sum
          # beta.array.unet_wsum_max[,aa,,] <- beta.array.unet_wsum
          
          # PI <- c(1, 5, 10, 50, 100)
          # threshold <- matrix(NA, length(PI), ncol(y))
          # for( colcol in 1:ncol(y) ){
          #   qhat <- ( sum( beta.array[,,,colcol,1] ) * K ) / (length(seq.alpha) * length(seq.lambda) * K )
          #   for( pipi in 1:length(PI) ){
          #     threshold[pipi, colcol] <-  qhat^2/(2*pipi*dim(beta.array)[1])+1/2
          #   }
          # }
          
        }
      }
      # end alpha ---------------------------------------------------------------
      
      end <- proc.time();
      cat("In this iteration, the elaspsed time=", (end-start)[3], "\n");
      cat("The remained time is", (end-start)[3]*(K-kk), "\n");
      
    }
    
    # end REP -----------------------------------------------------------------
    
    if( type.multivariate == "unet" ){
      GRID <- expand.grid( lapply( dim(beta.array.unet_sum)[-1], function(x) 1:x ) )
      tmp.denom1 <- tmp.denom2 <- array(0, c(max(GRID$Var1), max(GRID$Var2), max(GRID$Var3)))
      for( JJ in 1:nrow(GRID) ){
        if(JJ %% 10 == 0) cat("JJ is ", JJ, "\n")
        j1 <- GRID[JJ,1]
        j2 <- GRID[JJ,2]
        j3 <- GRID[JJ,3]
        tmp.denom1[j1,j2,j3] <- ifelse( max(beta.array.unet_sum[,j1,j2,j3]) == 0, 1, max(beta.array.unet_sum[,j1,j2,j3]) )
        # tmp.denom2[j1,j2,j3] <- ifelse( max(beta.array.unet_wsum[,j1,j2,j3]) == 0, 1, max(beta.array.unet_wsum[,j1,j2,j3]) )
        
        beta.array.unet_sum_max[,j1,j2,j3] <- as.array(beta.array.unet_sum[,j1,j2,j3]) / tmp.denom1[j1,j2,j3]
        # beta.array.unet_wsum_max[,j1,j2,j3] <- as.array(beta.array.unet_wsum[,j1,j2,j3]) / tmp.denom2[j1,j2,j3]
      }
      
      if( max(beta.array.unet) != 0 ) beta.array.unet_sum <- as.array(beta.array.unet_sum)/max(beta.array.unet_sum)
      beta.array.unet_sum <- slam::as.simple_sparse_array(beta.array.unet_sum)
      # if( max(beta.array.unet_wsum) != 0 ) beta.array.unet_wsum <- as.array(beta.array.unet_wsum)/max(beta.array.unet_wsum)
      
      
      
      out <- list(enet = beta.array.enet, 
                  unet = beta.array.unet,
                  unet_sum = beta.array.unet_sum,
                  unet_sum_max = beta.array.unet_sum_max
                  # unet_wsum = beta.array.unet_wsum,
                  # unet_wsum_max = beta.array.unet_wsum_max 
      )
    } else if (type.multivariate == "mnet"){
      out <- list(mnet = beta.array.mnet)
    }
    
    # beta.array %>% png.get_sp %>% head(15)
    # return( list(beta.array=beta.array, threshold=threshold) )
    
    # END sp ------------------------------------------------------------------
    
    return(out)
    
  }





# x <- rbind( mnormt::rmnorm(100, rep(0, 1000), diag(1000)), mnormt::rmnorm(100, rep(1, 1000), diag(1000)) )
# y <- data.frame(y1=rep(0:1, each=100), y2=x%*%rep(0.1, 1000)+rnorm(200), y3=rbinom(200, 5, prob=rep(c(0.4,0.6), each=100))) %>% .[,3]
# family <- c("binomial", "gaussian", "multinomial") %>% .[3]
# if(family=="multinomial"){
#   y <- model.matrix(~.-1, data=data.frame(as.factor(y)) )
# }
# alpha <- 1*0.1
# nlambda <- nrow(x)*0.5
# type.multinomial <- "grouped"
# DFDF <- 1
# DFDF <- nrow(x)
# df.seq <- floor(seq(0.1*nrow(y), 0.9*nrow(y), length.out=5))
# tol <- 1
# penalty.factor <- c(0,0,0,0, rep(1, ncol(x)-4))
# ... <- NULL


glmnet.df_to_lambda_old <- function(x, y, 
                                    alpha, nlambda,
                                    family, type.multinomial, DFDF, tol, 
                                    penalty.factor = penalty.factor, ...){
  
  fit.lambda_range <- glmnet( x=x, y=y, 
                              alpha=alpha, nlambda=nlambda,
                              family=family, type.multinomial = type.multinomial, 
                              penalty.factor = penalty.factor, ...)
  if(any(penalty.factor==0)) DFDF <- DFDF + sum(penalty.factor==0)
  
  df.loss <- abs(fit.lambda_range$df - DFDF)
  df_min <- fit.lambda_range$df[which.min(df.loss)]
  df_upper <- ifelse( all(!fit.lambda_range$df > df_min), 
                      max(fit.lambda_range$df), 
                      ifelse( which.min(df.loss)==length(df.loss), 
                              length(df.loss), 
                              min( fit.lambda_range$df[ fit.lambda_range$df > df_min ] ) ) )
  df_lower <- ifelse( which.min(df.loss)==1, min(fit.lambda_range$df), max( fit.lambda_range$df[ fit.lambda_range$df < df_min ] ) )
  
  lambda_min <- fit.lambda_range$lambda[which.min(df.loss)]
  lambda_upper <- min( fit.lambda_range$lambda[fit.lambda_range$df == df_upper] )
  lambda_lower <- max( fit.lambda_range$lambda[fit.lambda_range$df == df_lower] )
  
  fit.opt <- optimize( function(LAMBDA, DF, ...){
    fit2 <- glmnet( ..., lambda=LAMBDA ); 
    abs(fit2$df - DF)}, 
    interval = c(lambda_lower, lambda_upper), 
    DF=DFDF, x=x, y=y, alpha=alpha, family=family, type.multinomial=type.multinomial,
    tol=tol, 
    penalty.factor = penalty.factor, ... )
  
  glmnet( x=x, y=y,
          alpha=alpha, lambda=fit.opt$minimum,
          family=family, type.multinomial = type.multinomial,
          penalty.factor = penalty.factor, ...) %>% print
  
  opt.lambda=fit.opt$minimum
  opt.obj=fit.opt$objective
  list(opt.lambda = opt.lambda, 
       opt.obj = opt.obj)
  
}





library(glmnet)

glmnet.df_to_lambda <- function(x, y, 
                                alpha, 
                                family, type.multinomial, df.seq, 
                                penalty.factor = penalty.factor, verbose=TRUE){
  
  if( sum( penalty.factor==0 ) > min(df.seq) ) df.seq[which.min(df.seq)] <- df.seq[which.min(df.seq)] + sum( penalty.factor==0 )
  
  get_range_glmnet <- function(glmnet_fit, DFDF){
    
    df.loss <- abs(glmnet_fit$df - DFDF)
    out.df <- glmnet_fit$df[ which.min(df.loss) ]
    out.lambda <- glmnet_fit$lambda[ which.min(df.loss) ]
    
    if( sum( df.loss == 0 ) > 0 ){
      lambda_upper <- glmnet_fit$lambda[ median( which( df.loss == 0 ) ) ]
      lambda_lower <- glmnet_fit$lambda[ median( which( df.loss == 0 ) ) ]
    } else {
      lambda_minmax <- glmnet_fit$lambda[ order(df.loss, decreasing=FALSE)[1:2] ]
      lambda_upper <- max( lambda_minmax )
      lambda_lower <- min( lambda_minmax )
    }
    
    # 
    # df_min <- glmnet_fit$df[which.min(df.loss)]
    # df_upper <- ifelse( all(!glmnet_fit$df > df_min), 
    #                     max(glmnet_fit$df), 
    #                     ifelse( which.min(df.loss)==length(df.loss), 
    #                             length(df.loss), 
    #                             min( glmnet_fit$df[ glmnet_fit$df > df_min ] ) ) )
    # df_lower <- ifelse( which.min(df.loss)==1, min(glmnet_fit$df), max( glmnet_fit$df[ glmnet_fit$df < df_min ] ) )
    # 
    # lambda_min <- glmnet_fit$lambda[which.min(df.loss)]
    # lambda_lower <- min( glmnet_fit$lambda[glmnet_fit$df == df_upper] )
    # lambda_upper <- max( glmnet_fit$lambda[glmnet_fit$df == df_lower] )
    
    list( df = out.df, lambda = out.lambda, minmax=c(lower=lambda_lower, upper=lambda_upper) )
  }
  
  
  
  fit.lambda_range <- glmnet( x=x, y=y, 
                              alpha=alpha, nlambda=nrow(x)*5,
                              family=family, type.multinomial = type.multinomial, dfmax = 2*nrow(x),
                              penalty.factor = penalty.factor)
  
  lambda.sequence <- fit.lambda_range$lambda
  df.sequence <- fit.lambda_range$df
  
  
  out.lambda.seq <- NULL
  for( DFDF in df.seq ){
    fit.lambda_minmax <- get_range_glmnet(fit.lambda_range, DFDF)
    avg.lambda.minmax <- mean(fit.lambda_minmax$minmax)
    
    out.lambda.seq <- c(out.lambda.seq, avg.lambda.minmax)
  }
  
  
  
  if( verbose ){
    print(
      tmp_glmnet_fit <- 
        glmnet( x=x,
                y=y,
                alpha=alpha,
                lambda=out.lambda.seq,
                family=family,
                type.multinomial = type.multinomial,
                penalty.factor = penalty.factor) )
    
    cat( paste0("Target d.f.= ", df.seq, ", Estimated d.f.= ", tmp_glmnet_fit$df, collapse = "\n") )
    
  }
  
  list( DF=df.seq, df=tmp_glmnet_fit$df, lambda=out.lambda.seq )
  
  # count <- 0
  # while( TRUE ){
  #   
  #   count <- count + 1
  #   fit.lambda_minmax <- get_range_glmnet(fit.lambda_range)
  #   if( fit.lambda_minmax$df == DFDF ){
  #     if( verbose ){
  #       print(
  #         glmnet( x=x, 
  #                 y=y, 
  #                 alpha=alpha, 
  #                 lambda=fit.lambda_minmax$minmax[1],
  #                 family=family, 
  #                 type.multinomial = type.multinomial, 
  #                 penalty.factor = penalty.factor,
  #                 ...) )
  #     }
  #     
  #     return(fit.lambda_minmax)
  #     
  #     break
  #     
  #   } else {
  #     lambda.sequence <- fit.lambda_minmax$minmax %>% { seq( .[1], .[2], length.out = nlambda ) }
  #   }
  # 
  #   fit.lambda_range <- glmnet( x=x, 
  #                               y=y, 
  #                               alpha=alpha, lambda=lambda.sequence,
  #                               family=family, 
  #                               type.multinomial = type.multinomial, 
  #                               penalty.factor = penalty.factor,
  #                               ...)
  #  
  #   if( count > 20 | abs( fit.lambda_minmax$df - DFDF ) <= tol ) {
  #     if( verbose ){
  #       print(
  #         glmnet( x=x, 
  #                 y=y, 
  #                 alpha=alpha, 
  #                 lambda=fit.lambda_minmax$minmax[1],
  #                 family=family, 
  #                 type.multinomial = type.multinomial, 
  #                 penalty.factor = penalty.factor,
  #                 ...) )
  #     }
  #     
  #     return( fit.lambda_minmax )
  #     
  #     break
  #   }
  #      
  # }
  #  
  
}

# TIME.list <- NULL
# count <- 0
# for( lamlam in 10:50){
#   count <- count+1
#   TIME.list[[count]] <- system.time(
#     print(glmnet.df_to_lambda(x, 
#                                y, 
#                                nlambda = lamlam,
#                                alpha, 
#                                family, 
#                                type.multinomial, 
#                                1, 
#                                tol,
#                                penalty.factor))
#   )
#   print(TIME.list[[count]])
# }


# print(glmnet.df_to_lambda( x,
#                            y,
#                            alpha,
#                            family,
#                            type.multinomial,
#                            DFDF = 89,
#                            tol,
#                            penalty.factor,
#                            verbose = TRUE))
#




# x <- rbind( mnormt::rmnorm(100, rep(0, 1000), diag(1000)), mnormt::rmnorm(100, rep(1, 1000), diag(1000)) )
# y <- data.frame(y1=rep(0:1, each=100), y2=x%*%rep(0.1, 1000)+rnorm(200), y3=rbinom(200, 5, prob=rep(c(0.4,0.6), each=100)))
# # family <- rep("gaussian", ncol(y))
# family <- c("binomial", "gaussian", "multinomial")
# type.multivariate <- "unet"
# seq.alpha <- 1*0.1
# df.seq <- ceiling( seq(0.1*nrow(x), 0.9*nrow(x), length.out=10) )
# K <- 10
# psub <- 0.5
# 
# kk <- aa <- colcol <- 1
# penalty.factor <- c(0,0,0,0, rep(1, ncol(x)-4))
# ... <- NULL


png.cp_glmnet.df <-
  function(x,
           y,
           family,
           type.multivariate = c("unet", "mnet"),
           seq.alpha = NULL,
           df.seq = NULL,
           K = 100,
           psub = 0.5,
           penalty.factor = NULL,
           verbose = FALSE,
           tol) {
    
    library(mnormt)
    library(glmnet)
    library(dplyr)
    
    
    
    if( ! type.multivariate %in% c("enet", "unet", "mnet", "wnet") ) stop("type.multivariate %in% ('enet', 'unet', 'mnet', 'wnet') ")
    if( length(family) != ncol(y) ) stop("The length of family should be equal to ncol(y).")
    if( any(family %in% c("binomial", "multinomial")) & (type.multivariate %in% c("mnet")) ) stop("type.multivariate cannot be 'mnet' in family of binomial and multinomial")
    if( any(family %in% c("mgaussian")) ) stop("'mgaussian' should be used in type.multivariate")
    
    if( NROW(y) != nrow(x) ) stop("x and y should be equal length of row")
    if( is.null(seq.alpha)) seq.alpha <- 1:9*0.1
    if( type.multivariate=="mnet") standardize.response <- TRUE
    
    
    x <- as.matrix(x)
    y <- as.data.frame(y, stringsAsFactors = FALSE)
    
    n.pf <- sum(penalty.factor==0)
    if( n.pf > 0 ){
      wh.var <- (1:ncol(x))[-(1:n.pf)]
    } else {
      wh.var <- (1:ncol(x))
    }
    
    n <- nrow(x);
    p <- ncol(x)-n.pf;
    nsub <- n * psub
    
    n.alpha <- length(seq.alpha)
    df.seq <- sort(df.seq, decreasing = FALSE)
    n.df <- length(df.seq)
    
    # IF sp -------------------------------------------------------------------
    beta.array.enet <- array( 0, dim = c(p, n.alpha, n.df, ncol(y)),
                              dimnames = list(
                                paste0("", 1:p),
                                paste0("", seq.alpha),
                                paste0("", df.seq),
                                paste0("", seq_len(ncol(y)))
                              )
    )
    
    names(attributes(beta.array.enet)$dimnames) <-
      c("Variable", "Alpha", "Df", "Phenotype")
    
    beta.array.unet <-
      beta.array.unet_sum <- 
      beta.array.unet_sum_max <-
      beta.array.mnet <- 
      beta.array.enet
    
    beta.array.unet <- slam::as.simple_sparse_array(beta.array.unet)
    beta.array.unet_sum <- slam::as.simple_sparse_array(beta.array.unet_sum)
    beta.array.unet_sum_max <- slam::as.simple_sparse_array(beta.array.unet_sum_max)
    beta.array.mnet <- slam::as.simple_sparse_array(beta.array.mnet)
    beta.array.enet <- slam::as.simple_sparse_array(beta.array.enet)
    
    
    
    for (kk in 1:K) {
      print(kk)
      start <- proc.time()
      
      wsub <- sample(n, nsub)
      xsub <- x[wsub, , drop = F]
      ysub <- y[wsub, , drop = F]
      
      
      
      for (aa in 1:length(seq.alpha)) {
        if (type.multivariate == "mnet") {
          mgaussian.fit <- glmnet(x=xsub, y=ysub, 
                                  alpha=seq.alpha[aa], 
                                  lambda=unlist(seq.lambda[1]), 
                                  family="mgaussian",
                                  standardize.response=TRUE, 
                                  penalty.factor = penalty.factor)$beta[[1]][wh.var,]
          
          for (colcol in 1:ncol(y)) {
            beta.array.mnet[, aa, , colcol] <-
              beta.array.mnet[, aa, , colcol] + as.numeric(mgaussian.fit[[colcol]] != 0) / K
          }
          
        }
        
        
        
        if (type.multivariate == "unet") {
          # begin colcol ------------------------------------------------------------
          unet.matrix <- NULL
          for (colcol in 1:ncol(y)) {
            FAMILY <- family[colcol]
            
            if (FAMILY == "multinomial") {
              tmp.y <- model.matrix(~ . - 1, 
                                    data=as.data.frame(as.factor(ysub[, colcol, drop= T]) %>% droplevels()))
              type.multinomial <- "grouped"
            } else {
              tmp.y <- ysub[, colcol, drop = T]
              type.multinomial <- NULL
            }
            
            
            # start <- proc.time()
              fit_df_to_lambda <- glmnet.df_to_lambda( x = xsub, y = tmp.y,
                                                       alpha = seq.alpha[aa], 
                                                       # nlambda = 0.2*floor(nrow(xsub)),
                                                       # nlambda = 100,
                                                       family = FAMILY,
                                                       type.multinomial = type.multinomial,
                                                       df.seq = df.seq,
                                                       penalty.factor = penalty.factor, 
                                                       verbose = verbose )
              
              seq.lambda <- fit_df_to_lambda$lambda
              
            # end <- proc.time(); end-start
            
            
            fit <- glmnet( x = xsub, y = tmp.y,
                           alpha = seq.alpha[aa],
                           lambda = seq.lambda,
                           family = FAMILY,
                           type.multinomial = type.multinomial, penalty.factor = penalty.factor )
            
            
            fitted.beta <-
              switch(as.character(is.list(fit$beta)),
                     "TRUE" = fit$beta[[1]],
                     "FALSE" = fit$beta)
            
            fitted.beta <- fitted.beta[wh.var,]
              
            
            unet.matrix <- cbind(unet.matrix, as.numeric(fitted.beta != 0))
            
            beta.array.enet[, aa, , colcol] <- 
              as.array(beta.array.enet[, aa, , colcol]) + as.numeric(fitted.beta != 0) / K
          }
          # end colcol --------------------------------------------------------------
          
          for (colcol in 1:ncol(y)) {
            beta.array.unet[, aa, , colcol] <- as.array( beta.array.unet[, aa, , colcol] ) + apply(unet.matrix, 1, function(x) ifelse(sum(x) > 0, 1, 0)) / K
            beta.array.unet_sum[, aa, , colcol] <- as.array( beta.array.unet_sum[, aa, , colcol] ) + apply(unet.matrix, 1, function(x) sum(x)) / K
          }
          
          
          # PI <- c(1, 5, 10, 50, 100)
          # threshold <- matrix(NA, length(PI), ncol(y))
          # for( colcol in 1:ncol(y) ){
          #   qhat <- ( sum( beta.array[,,,colcol,1] ) * K ) / (length(seq.alpha) * length(seq.lambda) * K )
          #   for( pipi in 1:length(PI) ){
          #     threshold[pipi, colcol] <-  qhat^2/(2*pipi*dim(beta.array)[1])+1/2
          #   }
          # }
          
        }
      }
      
      # end alpha ---------------------------------------------------------------
      
      end <- proc.time()
      
      cat("In this iteration, the elaspsed time=", (end -
                                                      start)[3], "\n")
      
      cat("The remained time is", (end - start)[3] *
            (K - kk), "\n")
      
      
    }
    
    # end REP -----------------------------------------------------------------
    
    if (type.multivariate == "unet") {
      GRID <- expand.grid( lapply(dim(beta.array.unet_sum)[-1], function(x) 1:x) )
      tmp.denom1 <- array(0, c(max(GRID$Var1), max(GRID$Var2), max(GRID$Var3)))
      for (JJ in 1:nrow(GRID)) {
        if (JJ %% 10 == 0) cat("JJ is ", JJ, "\n")
        j1 <- GRID[JJ, 1]
        j2 <- GRID[JJ, 2]
        j3 <- GRID[JJ, 3]
        
        tmp.denom1[j1, j2, j3] <- ifelse(max(beta.array.unet_sum[, j1, j2, j3]) == 0, 1, max(beta.array.unet_sum[, j1, j2, j3]))
        beta.array.unet_sum_max[, j1, j2, j3] <- as.array(beta.array.unet_sum[, j1, j2, j3]) / tmp.denom1[j1, j2, j3]
      }
      
      beta.array.unet_sum <- as.array(beta.array.unet_sum) / max(beta.array.unet_sum)
      beta.array.unet_sum <- slam::as.simple_sparse_array(beta.array.unet_sum)
      
      out <- list( enet = beta.array.enet,
                   unet = beta.array.unet,
                   unet_sum = beta.array.unet_sum,
                   unet_sum_max = beta.array.unet_sum_max#,
                   # unet_2toq = beta.array.unet_2toq,
                   # unet_3toq = beta.array.unet_3toq,
                   # unet_2 = beta.array.unet_2,
                   # unet_3 = beta.array.unet_3,
                   # unet_4 = beta.array.unet_4
      )
      
    } else if (type.multivariate == "mnet") {
      
      out <- list(mnet = beta.array.mnet)
      
    }
    
    # beta.array %>% png.get_sp %>% head(15)
    # return( list(beta.array=beta.array, threshold=threshold) )
    
    # END sp ------------------------------------------------------------------
    
    return(out)
    
  }

# beta.array.unet_sum_max %>% apply(c(1,4), max) %>% .[,1] %>% order(decreasing=TRUE) %>% head(20)
# beta.array.unet %>% apply(c(1,4), max) %>% .[,1] %>% order(decreasing=TRUE) %>% head(20)
# beta.array.mnet %>% apply(c(1,4), max) %>% .[,1] %>% order(decreasing=TRUE) %>% head(20)
# beta.array.unet_sum_max[1:500,1,,1] %>%
#   gplots::heatmap.2(dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none', col=rev(heat.colors(16)))


png.get_lambda <-
  function(x, y,
           family,
           iter = 10,
           seq.alpha = NULL,
           n.lambda = NULL,
           psub = 0.5,
           ...) {
    # x=Data$snp
    # y=Data$y
    # family="mixed"
    # psub=0.5
    # seq.alpha=seqalpha
    # n.lambda=nlambda
    # setseed=1129
    
    if (NROW(y) != nrow(x))
      stop("x and y should be equal length of row")
    
    if (is.null(seq.alpha)) seq.alpha <- 1:9 * 0.1
    if (is.null(n.lambda)) n.lambda <- 10
    
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    nsub <- n * psub
    
    y.column.set <- seq_len(ncol(y))
    
    if (length(family) > 1) {
      lambda.array <- array( NA, dim = c(iter, n.lambda, length(seq.alpha), ncol(y), 2),
                             dimnames = list(
                               paste0("iter=", 1:iter),
                               paste0("lambda", 1:n.lambda),
                               paste0("alpha=", seq.alpha),
                               paste0("Y", 1:ncol(y)),
                               c("lambda", "df")
                             )
      )
      
      seq.lambda <- as.list(1:length(y.column.set))
      
      for (colcol in 1:length(y.column.set)) {
        FAMILY <- family[colcol]
        h <- unlist(y.column.set[colcol])
        
        lambda.vec <- NULL
        for (i in 1:iter) {
          for (j in 1:length(seq.alpha)) {
            wsub <- sample(n, nsub)
            xsub <- x[wsub, , drop = F]
            ysub <- y[wsub, , drop = F]
            if (FAMILY == "multinomial") {
              fitsub <- glmnet(
                x = xsub,
                y = model.matrix( ~ . - 1, data = as.data.frame(as.factor(ysub[, h]) %>% droplevels()) ),
                alpha = seq.alpha[j],
                family = FAMILY,
                nlambda = n.lambda,
                type.multinomial = "grouped",
                ... )
            } else {
              fitsub <- glmnet(
                x = xsub,
                y = ysub[, h],
                alpha = seq.alpha[j],
                family = FAMILY,
                nlambda = n.lambda,
                ...
              )
            }
            
            
            lambda.array[i, 1:length(fitsub$df), j, colcol, 1] <- fitsub$lambda
            lambda.array[i, 1:length(fitsub$df), j, colcol, 2] <- fitsub$df
            
            lambda.vec <- c(lambda.vec, fitsub$lambda)
          }
        }
        
        seq.lambda[[colcol]] <- lambda.vec
        
      }
      
    } else if (family == "mgaussian") {
      seq.lambda <- as.list(1:length(y.column.set))
      
      lambda.vec <- NULL
      for (i in 1:iter) {
        for (j in 1:length(seq.alpha)) {
          wsub <- sample(n, nsub)
          xsub <- x[wsub, , drop = F]
          ysub <- y[wsub, , drop = F]
          fitsub <-
            glmnet(
              x = xsub,
              y = ysub,
              alpha = seq.alpha[j],
              family = family,
              nlambda = n.lambda,
              ...
            )
          
          lambda.vec <- c(lambda.vec, fitsub$lambda)
        }
      }
      
      for (colcol in 1:length(y.column.set)) {
        seq.lambda[[colcol]] <- lambda.vec
        
      }
      
    } else {
      seq.lambda <- as.list(1:length(y.column.set))
      for (colcol in 1:length(y.column.set)) {
        h <- unlist(y.column.set[colcol])
        
        lambda.vec <- NULL
        for (i in 1:iter) {
          for (j in 1:length(seq.alpha)) {
            wsub <- sample(n, nsub)
            xsub <- x[wsub, , drop = F]
            ysub <- y[wsub, , drop = F]
            fitsub <- glmnet( x = xsub, y = ysub[, h],
                              alpha = seq.alpha[j],
                              family = family,
                              nlambda = n.lambda,
                              ... )
            lambda.vec <- c(lambda.vec, fitsub$lambda)
          }
        }
        seq.lambda[[colcol]] <- lambda.vec
        
      }
    }
    
    return(seq.lambda)
    
  }


png.get_sp <- function(array, K = 100) {
  Margin.rep <- which(!names(dimnames(array)) %in% c("Replications"))
  count.array <- apply(array, Margin.rep, mean)
  Margin.tuning <- which(!names(dimnames(count.array)) %in% c("Alpha", "Lambda"))
  out.sp <- as.matrix(apply(count.array, Margin.tuning, max))
  
  c("Variable", "Alpha", "Lambda", "Phenotype", "Replications")
  
  Margin.phenotype <- which(names(dimnames(count.array)) %in% c("Phenotype"))
  
  PI.seq <- c(1, 5, 10, 50, 100)
  threshold <- matrix(NA, length(PI.seq), dim(count.array)[Margin.phenotype])
  dimnames(threshold) <- list(PI.seq, paste0("y", seq_len(ncol(threshold))))
  for (colcol in 1:ncol(threshold)) {
    qhat <- (sum(count.array[, , , colcol]) * K) / (prod(dim(count.array)[-Margin.tuning]) * K)
    for (pipi in 1:length(PI.seq)) {
      PI <- PI.seq[pipi]
      threshold[pipi, colcol] <- min(qhat ^ 2 / (2 * PI * dim(count.array)[1]) + 1 / 2, 1)
    }
  }
  
  list(sp = out.sp, threshold = threshold)
}



































png.venn2 = function(x,y,duplicated=FALSE){
  inner = intersect(x,y)
  xnyc = x[!x%in%inner]
  ynxc = y[!y%in%inner]
  if(!duplicated){
    inner = unique(inner)
    xnyc = unique(xnyc)
    ynxc = unique(ynxc)
  }
  list( x = xnyc , y = ynxc , inner = inner )
}

