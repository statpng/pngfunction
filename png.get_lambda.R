library(mnormt)
library(glmnet)
library(dplyr)

# png.get_lambda <- function(x, y, family, iter=10, seq.alpha=NULL, n.lambda=NULL, psub=0.5, ...){
#         # x=Data$snp
#         # y=Data$y
#         # family="mixed"
#         # psub=0.5
#         # seq.alpha=seqalpha
#         # n.lambda=nlambda
#         # setseed=1129
#         
#         if( NROW(y) != nrow(x) ) stop("x and y should be equal length of row")
#         
#         if(is.null(seq.alpha)) seq.alpha <- 1:9*0.1
#         if(is.null(n.lambda)) n.lambda <- 10
#         
#         x <- as.matrix(x)
#         y <- as.matrix(y)
#         
#         n <- nrow(x);
#         p <- ncol(x);
#         nsub <- n*psub;
#         
#         if( family == "mgaussian" ){
#                 y.column.set <- list(seq_len(ncol(y)))
#         } else {
#                 y.column.set <- seq_len(ncol(y))
#         }
#         
#         seq.lambda <- as.list(1:length(y.column.set))
#         for( colcol in 1:length(y.column.set) ){
#                 h <- unlist( y.column.set[colcol] )
#                 
#                 lambda.vec <- NULL
#                 for( i in 1:iter ){
#                         for( j in 1:length(seq.alpha) ){
#                                 wsub <- sample(n, nsub)
#                                 xsub <- x[wsub,,drop=F]
#                                 ysub <- y[wsub,,drop=F]
#                                 fitsub <- glmnet(x=xsub, y=ysub[,h], alpha=seq.alpha[j], family=family, nlambda = n.lambda, ...)
#                                 lambda.vec <- c( lambda.vec, fitsub$lambda )
#                         }
#                 }
#                 seq.lambda[[colcol]] <- lambda.vec;
#         }
#         
#         
#         return(seq.lambda)
#         
# }


png.scale <- function(mat){
        mat <- as.matrix(mat)
        out <- mat
        sd.vec <- apply(mat, 2, function(x) sqrt(mean((x-mean(x))^2)) )
        
        out[,sd.vec != 0] <- t(t(scale( mat[,sd.vec != 0], scale=FALSE ))/(sd.vec[sd.vec != 0]))
        out[,sd.vec == 0] <- scale( mat[,sd.vec == 0], scale=FALSE, center=FALSE )
        out
}


png.lambda <- function(x, y, seq.alpha, family){
        x.std <- png.scale(x)
        y.std <- as.matrix(y) # png.scale(y)
        if(family=="multinomial") y.std <- model.matrix(~-1+as.factor(y.std))
        # In glmnet, you should set the type.multinomial to multinomial.
        
        n <- nrow(y.std)
        p <- ncol(x)
        lambda.min.ratio <- ifelse( n < p, 0.01, 0.0001 )
        K <- 100
        
        out <- as.list(1:length(seq.alpha))
        names(out) <- seq.alpha
        for( h in 1:length(seq.alpha) ){
                alpha = seq.alpha
                lambda.max <- max( apply( crossprod(y.std, x.std), 2, 
                                          function(z) sqrt(sum(z^2))/(alpha*n) )) 
                lambda.min <- lambda.max*lambda.min.ratio
                out[[h]] <- rev( seq( lambda.min, lambda.max, length.out = K ) )
        }
        
        out
                
        
        # png.lambda(x = iris[,1:2], y = iris$Species, family = "multinomial", seq.alpha = .1)
        # glmnet(x = as.matrix(iris[,1:2]), y = iris$Species, family = "multinomial", type.multinomial = "grouped", alpha = .1)$lambda
        
        # png.lambda(x = iris[,1:2], y = iris$Petal.Width, family = "gaussian", seq.alpha = .1)
        # glmnet(x = as.matrix(iris[,1:2]), y = iris$Petal.Width, family = "gaussian", type.multinomial = "grouped", alpha = .1)$lambda
        
        # png.lambda(x = iris[,1:2], y = iris[,3:4], family = "mgaussian", seq.alpha = 1:9*0.1)
        # glmnet(x = as.matrix(iris[,1:2]), y = iris[,3:4], family = "mgaussian", type.multinomial = "grouped", alpha = .1)$lambda
}


png.sp.lambda <- function(x, y, family, iter=10, seq.alpha=NULL, n.lambda=NULL, psub=0.5){
        # x=Data$snp
        # y=Data$y
        # family="mixed"
        # psub=0.5
        # seq.alpha=seqalpha
        # n.lambda=nlambda
        # setseed=1129

        if( NROW(y) != nrow(x) ) stop("x and y should be equal length of row")

        if(is.null(seq.alpha)) seq.alpha <- 1:9*0.1
        if(is.null(n.lambda)) n.lambda <- 10

        x <- as.matrix(x)
        y <- as.matrix(y)

        n <- nrow(x);
        p <- ncol(x);
        nsub <- n*psub;

        y.column.set <- seq_len(ncol(y))

        seq.lambda <- as.list(1:length(y.column.set))
        for( colcol in 1:length(y.column.set) ){
                h <- unlist( y.column.set[colcol] )

                lambda.vec <- NULL
                for( i in 1:iter ){
                        
                        wsub <- sample(n, nsub)
                        xsub <- x[wsub,,drop=F]
                        ysub <- y[wsub,,drop=F]
                        lambda.list <- png.lambda(x=xsub, y=ysub[,h], seq.alpha=seq.alpha, family=family)
                        lambda.vec <- c( lambda.vec, unlist(lambda.list) )
                        
                }
                seq.lambda[[colcol]] <- lambda.vec;
        }


        return(seq.lambda)

}







png.get_lambda.sgl <- function(x, y, index=NULL, type, iter=10, seq.alpha=NULL, n.lambda=NULL, psub=0.5, ...){
        # x=Data$snp
        # y=Data$y
        # type="linear"
        # psub=0.5
        # seq.alpha=seqalpha
        # n.lambda=nlambda
        # setseed=1129
        
        if( NROW(y) != nrow(x) ) stop("x and y should be equal length of row")
        
        if(is.null(seq.alpha)) seq.alpha <- 1:9*0.1
        if(is.null(n.lambda)) n.lambda <- 10
        
        if(type=="mgaussian") type <- "linear"
        if(type=="gaussian") type <- "linear"
        if(type=="binomial") type <- "logit"
        
        x <- as.matrix(x)
        y <- as.matrix(y)
        
        n <- nrow(x);
        p <- ncol(x);
        q <- ncol(y)
        nsub <- n*psub;
        
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
        
        if(q>1) index <- fit.kronecker$group.xy
        
        n.new <- nrow(x);
        p.new <- ncol(x);
        q.new <- ncol(y)
        
        seq.lambda <- as.list(1:ncol(y))
        lambda.vec <- NULL
        for( i in 1:iter ){
                for( j in 1:length(seq.alpha) ){
                        wsub <- sample(n, nsub)
                        wsub.new <- unlist( lapply( (seq_len(q)-1)*n, function(x) x + wsub ) )
                        xsub <- x[wsub.new,,drop=F]
                        ysub <- y[wsub.new,,drop=F]
                        
                        fitsub <- SGL(data=list(x=xsub, y=ysub), 
                                      index=index, 
                                      type=type, 
                                      alpha=seq.alpha[j], 
                                      nlam=n.lambda, verbose=TRUE)
                        
                        lambda.vec <- c( lambda.vec, fitsub$lambdas )
                }
        }
        
        seq.lambda[[1]] <- lambda.vec
        
        return(seq.lambda)
        
}
