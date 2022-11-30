
# vec <- pseq@otu_table@.Data[,1]
AlphaDiversity <- function(vec){
  vec <- vec[vec>0]
  p <- length(vec)
  pj <- vec / sum(vec)
  richness = sum(vec>0)
  shannon = -sum(pj * log(pj))
  pielou = shannon / log(p)
  
  list( richness=richness, shannon=shannon, pielou=pielou )
}

png.AlphaDiversity <- function(mat){
  n <- ncol(mat)
  sapply(1:n, function(j) AlphaDiversity(abundances(pseq)[,j])) %>% t
}

# x1 <- pseq@otu_table@.Data[,1]
# x2 <- pseq@otu_table@.Data[,2]
BetaDiversity <- function(x1, x2){
  wh1 <- which(x1>0)
  wh2 <- which(x2>0)
  
  Jaccard = 1 - length(intersect(wh1,wh2))/length(unique(c(wh1,wh2)))
  BrayCurtis = sum(abs(x1-x2))/(sum(x1)+sum(x2))
  
  list( jaccard=Jaccard, braycurtis=BrayCurtis )
}

png.BetaDiversity <- function(mat, type="braycurtis"){
  # mat: p x n
  n <- ncol(mat)
  
  BetaDivMat <- matrix(0,n,n)
  for(u in 1:(n-1)){
    for(v in (u+1):n){
      BetaDivMat[u,v] <- BetaDiversity(mat[,u], mat[,v])[[type]]
    }
  }
  BetaDivMat[lower.tri(BetaDivMat)] <- t(BetaDivMat) %>% .[lower.tri(.)]
  cat( "IsSymmetric:", Matrix::isSymmetric(BetaDivMat), "\n" )
  
  BetaDivMat
}


if(FALSE){
## Validation
# BiocManager::install("microbiome")
library(microbiome)
data(dietswap)
pseq <- dietswap

# Package
microbiome::alpha(pseq, index = c("diversity_shannon","evenness_pielou"))[1:10,]
mapply(function(u,v) abundances(pseq) %>% {microbiome::divergence(.[,u], .[,v])}, u=1:4, v=2:5)

# My Function
sapply(1:10, function(j) AlphaDiversity(abundances(pseq)[,j])) %>% t
mapply(function(u,v) abundances(pseq) %>% {BetaDiversity(.[,u], .[,v])}, u=1:4, v=2:5)



## Matrix form
Alpha <- png.AlphaDiversity(abundances(pseq))
Beta <- png.BetaDiversity(abundances(pseq))

## PCoA
fit.pcoa <- ape::pcoa(Beta)
plot( fit.pcoa$vectors[,1:2], pch=18 )
}
