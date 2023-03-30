png.cmdscale <- function(dist, dim, eig=TRUE, ...){
    dist <- as.matrix(dist)
    
    if( length(rownames(dist))>0 ){
        Label <- rownames(dist)
    } else {
        Label <- seq_len(nrow(dist))
    }
    
    cmdscale( dist, k = max(dim), eig=eig, ... ) %$% 
        {
            .$points %>% 
                magrittr::set_colnames(paste0("Dim", seq_len(max(dim)) )) %>% 
                data.frame(Label=Label, .) %>% 
                as.data.frame %>% 
                ggplot( aes_string(x=paste0("Dim", dim[1]), y=paste0("Dim", dim[2])) )+
                xlab(paste0("Dim", dim[1], " (", round(.$eig[dim[1]]/sum(.$eig),3)*100 ,"%)"))+
                ylab(paste0("Dim", dim[2], " (", round(.$eig[dim[2]]/sum(.$eig),3)*100 ,"%)"))+
                geom_point(size=1, shape=3)+
                geom_text(aes(label=Label), size=5,
                          check_overlap = FALSE,
                          # hjust = "center", vjust = "bottom",
                          nudge_x = 0, nudge_y = 0.05)+
                geom_hline(yintercept=0, lty=2)+
                geom_vline(xintercept=0, lty=2)+
                theme_bw()
        }
}


png.isoMDS <- function(dist, dim, ...){
    library(MASS)
    dist <- as.matrix(dist)
    
    if( length(rownames(dist))>0 ){
        Label <- rownames(dist)
    } else {
        Label <- seq_len(nrow(dist))
    }
    
    isoMDS( dist+1e-5, k = ceiling(nrow(dist)/2), ... ) %$% 
        {
            .$points %>% 
                magrittr::set_colnames(paste0("Dim", seq_len(ncol(.)) )) %>% 
                data.frame(Label=Label, .) %>% 
                as.data.frame %>% 
                ggplot( aes_string(x=paste0("Dim", dim[1]), y=paste0("Dim", dim[2])) )+
                xlab(paste0("Dim", dim[1], " (", round(apply(.$points, 2, function(x) sd(x)^2 )[dim[1]]/sum(apply(.$points, 2, function(x) sd(x)^2 )),3)*100 ,"%)"))+
                ylab(paste0("Dim", dim[2], " (", round(apply(.$points, 2, function(x) sd(x)^2 )[dim[2]]/sum(apply(.$points, 2, function(x) sd(x)^2 )),3)*100 ,"%)"))+
                geom_point(size=1, shape=3)+
                geom_text(aes(label=Label), size=5,
                          check_overlap = FALSE,
                          # hjust = "center", vjust = "bottom",
                          nudge_x = 0, nudge_y = 0.05)+
                geom_hline(yintercept=0, lty=2)+
                geom_vline(xintercept=0, lty=2)+
                theme_bw()
        }
}



png.sammon <- function(dist, dim, ...){
    
    library(MASS)
    dist <- as.matrix(dist)
    
    if( length(rownames(dist))>0 ){
        Label <- rownames(dist)
    } else {
        Label <- seq_len(nrow(dist))
    }
    
    con <- sammon(dist, k=max(dim), ...) %$% 
        {
            .$points %>% 
                magrittr::set_colnames(paste0("Dim", seq_len(max(dim)) )) %>% 
                data.frame(Label=Label, .) %>% 
                as.data.frame %>% 
                ggplot( aes_string(x=paste0("Dim", dim[1]), y=paste0("Dim", dim[2])) )+
                xlab(paste0("Dim", dim[1], " (", round(.$eig[dim[1]]/sum(.$eig),3)*100 ,"%)"))+
                ylab(paste0("Dim", dim[2], " (", round(.$eig[dim[2]]/sum(.$eig),3)*100 ,"%)"))+
                geom_point(size=1, shape=3)+
                geom_text(aes(label=Label), size=5,
                          check_overlap = FALSE,
                          # hjust = "center", vjust = "bottom",
                          nudge_x = 0, nudge_y = 0.05)+
                geom_hline(yintercept=0, lty=2)+
                geom_vline(xintercept=0, lty=2)+
                theme_bw()
        }
}



png.get_dim <- function(Vector){
  
  equation <- function(n){
    abs( n*(n-1)/2 - length(Vector) )
  }
  
  round(optimize(equation, c(0,100))$minimum, 0)
}

png.trans_to_dist <- function(mat){
  J <- matrix(1, nrow(mat), nrow(mat))
  sqrt( diag(mat)*J - 2*mat + t( diag(mat)*J ) )
}

png.get_symmetric <- function(x, diag=FALSE){
  
  if( diag ){
    nc <- nr <- png.get_dim(x)-1
    diags <- x[ cumsum(cumsum(rep(1,nc))) ]
    x <- x[ -cumsum(cumsum(rep(1,nc))) ]
  } else {
    nc <- nr <- png.get_dim(x)
  }
  
  out <- matrix(0, nrow=nr, ncol=nc)
  out[upper.tri(out)] <- x
  out[lower.tri(out)] <- t(out)[lower.tri(out)]
  
  if( diag ) diag(out) <- diags
  
  return(out)
}


png.get_mds <- function(mat, k=2, method="cmdscale", labels=NULL, print.plot=TRUE, ...){
  library(MASS)
  con <- switch( method, 
                 cmdscale = cmdscale( mat, k=k, eig=TRUE, ...),
                 isoMDS = isoMDS(mat, k=k, ...),
                 sammon = sammon(mat, k=k, ...)
  )
  if( print.plot ){
    x<-con$points[,1]
    y<-con$points[,2]
    lim<-c(-max(abs(con$points))*1.1,
           max(abs(con$points))*1.1)
    plot(x, y, xlab="Dim1", ylab="Dim2", xlim=lim, ylim=lim)
    text(x, y, labels, cex=0.8,pos=1)
    abline(v=0, h=0)
  }
  
  return(append(list(dist=mat), con) )
}

png.plot_mds <- function( df, labels, ... ){
    x<-df[,1]
    y<-df[,2]
    lim<-c(-max(abs(df))*1.1,
           max(abs(df))*1.1)
    # plot(x,y)
    plot(x, y, xlim=lim, ylim=lim, ...)
    text(x, y, labels, cex=0.8,pos=1, ...)
    abline(v=0, h=0)
}

png.get_stress <- function(mat, method = "cmdscale", labels = NULL, k=ncol(mat)-1, ...){
  stress_vs_dim <- NULL
  for( kk in 1:k ){
    res <- png.get_mds(mat=mat, k=kk, method = method, labels = method, print.plot=FALSE, ...)
    delta <- as.matrix( dist(x = res$points, method = "euclidean", diag=TRUE, upper=TRUE) )
    stress <- sqrt(sum((( res$dist - delta )[lower.tri(res$dist)])^2)/
                     sum(res$dist[lower.tri(res$dist)]^2))
    stress_vs_dim[[kk]] <- stress
  }
  plot( 1:k, stress_vs_dim, type="b" )
  
  return(Kruskal_stress = stress_vs_dim)
}

png.rotation <- function(mat, theta){
  rot <-
    rbind(
      c( cos(theta), -sin(theta) ),
      c( sin(theta), cos(theta) )
    )
  t( rot %*% t(mat) )
}



png.similarity <- function(df, method=c("simple", "RR", "jacarrd")){
  p <- ncol(df)
  n <- nrow(df)
  J <- matrix(1, n, p)
  
  if( method == "simple" ) out <- ( df %*% t(df) + (J-df)%*%t(J-df) ) / p
  if( method == "RR" ) out <- df %*% t(df) / p
  if( method == "jaccard" ) out <- (df %*% t(df)) / (p-(df==0) %*% t(df==0))
  
  return(out)
}


png.mahalanobis <- function(df){

N <- nrow(df)
COV <- cov(df)
out <- matrix( NA, N, N )

for(i in 1:N){
  for(j in 1:N){
    out[i,j]<-sqrt(t(df[i,]-df[j,])%*% COV %*% (df[i,]-df[j,]))
  }}
rownames(out)<-colnames(out)<-rownames(df)
out
}




library(dplyr)


Torgerson <- function(D, add=FALSE){
  D <- as.matrix(D)
    n <- nrow(D)
    Ones <- matrix(1,n,n)
    H <- diag(n) - Ones/n
    B1 <- -1/2 * (H %*% D %*% H)
    B2 <- -1/2 * (H %*% D^2 %*% H)
    
  if( add ){
    tmp.mat <- rbind(  cbind(matrix(0,n,n), 2*B2),
                       cbind(-diag(n), -4*B1)  )
    ca <- max( Re( eigen(tmp.mat)$values ) )
    PSD.mat <- -1/2 * ( D + ca*(1-diag(n)) )^2
    PSD.cnt.mat <- H %*% PSD.mat %*% H
    eigen.PSD <- eigen( PSD.cnt.mat, symmetric = TRUE )
    positive.ev.bool <- eigen.PSD$values > 0
    
    ev <- eigen.PSD$values[positive.ev.bool]
    evec <- eigen.PSD$vectors[ , positive.ev.bool]
    out <- evec %*% diag( sqrt(ev) )
    dimnames(out) <- list(rownames(D), NULL)
  } else {
    out <- B2
  }
  return(out)
}

# Torgerson(dist(skull[1:30,-1]), add=TRUE)[,1:10] %>% head
# cmdscale(dist(skull[1:30,-1]), k=10,add=TRUE)$points[,1:10] %>% head

