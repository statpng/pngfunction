
detach("package:png.compositionalPCA",unload=TRUE)
library(png.compositionalPCA)
library(parallel)
library(purrr)
library(tidyverse)
dirmult
Ternary
c("reshape2", "Rtsne", "mnormt", "quadprog")



{
  

  {
    n.seq <- c(50,100,500,1000)
    p.seq <- c(100,500)
    eta.seq <- c(0, 0.1)
    seed.seq <- c(1:10)
    param.list <- list(n.seq=n.seq,
                       p.seq=p.seq,
                       eta.seq=eta.seq,
                       seed.seq=seed.seq)

    # Create the indices
    grid <- expand.grid(param.list)
  }
}



# consistency
function(){
  {
    #1: eta=0
    eta=0

    Fnorm_V <- function(V,vhat){
      if(FALSE){
        V=data$V; vhat=fit$vhat
      }
      r=NCOL(V)
      out <- NULL
      for( k in 1:r ){
        out[k] <- min(sum((V[,k]-vhat[,k])^2),
                      sum((V[,k]+vhat[,k])^2))
      }
      sqrt(sum(out))
    }



    n.seq <- c(10,100,1000,5000)
    eta.seq <- c(-0.1, 0, 0.1)
    seed.seq <- 1:10

    out <- array(NA, dim=c(4,3,2,3,10), dimnames=list( paste0("n=",c(10,100,1000,5000)), c("V vs PC", "V vs GPC", "PC vs GPC"), c("angle", "fnorm"), paste0("eta=", c(-0.1, 0, 0.1)), paste0("seed=", seed.seq) ))

    for( k in 1:length(seed.seq) ){
      seed <- seed.seq[k]
      for( j in 1:length(eta.seq) ){
        eta=eta.seq[j]

        for( i in 1:length(n.seq) ){
          n=n.seq[i]

          p=10; r=2; d=1; d0=0.1; snr=5; seed=seed
          data <- sim.simplex(n=n,p=p,r=r,snr=snr,d=d,d0=d0,seed=seed,eta=eta )
          X <- data$X2
          fit <- png.gppca_qp(X, nrank=r, kappa=1e-4)

          out[i,1,1,j,k] <- png.utils::png.angle(data$V, prcomp(X)$rot[,1:2])$max
          out[i,2,1,j,k] <- png.utils::png.angle(data$V, fit$vhat)$max
          out[i,3,1,j,k] <- png.utils::png.angle(prcomp(X)$rot[,1:2], fit$vhat)$max

          out[i,1,2,j,k] <- Fnorm_V(data$V, prcomp(X)$rot[,1:2])
          out[i,2,2,j,k] <- Fnorm_V(data$V, fit$vhat)
          out[i,3,2,j,k] <- Fnorm_V(prcomp(X)$rot[,1:2], fit$vhat)
        }
      }
    }

    apply(out,1:4,mean)
    #
    #

  }

}








function(){

  {


    {
      n=100; p=4; r=2; snr=5; seed=1; eta=0
      X <- sim.simplex(n=n,p=p,r=r,snr=snr,d=10,d0=0.01,seed=seed,eta=eta)$X2
      fit1_a <- png.lrpca(X, nrank=r, zero.replace="simple")
      fit1_b <- png.lrpca(X, nrank=r, zero.replace="additive")
      fit1_c <- png.lrpca(X, nrank=r, zero.replace="multiplicative")
      # png.quaternary3d(fit1$Xnew, vhat=fit1$vhat, xhat=fit1$xhat)
      fit2 <- png.ppca(X, nrank=r)
      fit3 <- png.gppca(X, nrank=r)
      kappa=1e-4; gamma=0.2
      fit4 <- png.ppca_qp(X, nrank=r, kappa=kappa, gamma=gamma)
      fit5 <- png.gppca_qp(X, nrank=r, kappa=kappa, gamma=gamma)

      png.ternary(fit1_a$Xnew, vhat=fit1_a$vhat)
      png.ternary(X, vhat=fit4$vhat)
      png.ternary(X, vhat=fit5$vhat)

      png.ternary(fit1_a$Xnew, vhat=fit1_a$vhat)
      png.ternary(X, vhat=fit4$vhat)
      png.quaternary3d(X, vhat=fit4$vhat)
      png.quaternary3d(X, vhat=fit5$vhat)
    }




    # main function
    sim1_convergence <- function(n, p, r, snr, eta, seed){
      if(FALSE){
        n=100; p=4; r=2; snr=2; eta=0.1/log(p); seed=1
        n=100; p=50; r=5; snr=5; eta=0.0; seed=1
      }

      data <- sim.simplex(n=n,p=p,r=r,snr=snr,d=10,d0=0.01,seed=seed,eta=eta/log(p))
      X <- data$X2

      kappa=1e-6; gamma=0.5
      fit4 <- png.ppca_qp(X, nrank=r, kappa=kappa, maxit=500, eps=1e-6, gamma=gamma)
      fit5 <- png.gppca_qp(X, nrank=r, kappa=kappa, maxit=500, eps=1e-6, gamma=gamma)

      # png.pca.plot_convergence(fit4)
      # png.pca.plot_convergence(fit5)

      return( list(fit4, fit5) %>% purrr::map(~png.pca.convergence(.x)) )

    }




    res1_convergence <- function(){
      # n.seq eta.seq seed.seq
      # 1     10     0.0        1
      # 2    100     0.0        1
      # 3    500     0.0        1
      # 4   1000     0.0        1
      # 5     10     0.1        1
      # 6    100     0.1        1
      # 7    500     0.1        1
      # 8   1000     0.1        1
      # 9     10     0.0        2
      # 10   100     0.0        2
      # 11   500     0.0        2
      # 12  1000     0.0        2
      # 13    10     0.1        2
      # 14   100     0.1        2
      # 15   500     0.1        2
      # 16  1000     0.1        2

      {
        n.seq <- c(50,100,500)
        p.seq <- c(100,200)
        eta.seq <- c(0, 0.1)
        seed.seq <- c(1:10)
        param.list <- list(n.seq=n.seq,
                           p.seq=p.seq,
                           eta.seq=eta.seq,
                           seed.seq=seed.seq)
        grid <- expand.grid(param.list)

        res1 <- mcmapply(function(n, p, r, snr, eta, seed) sim1_convergence(n=n, p=p, r=r, snr=snr, eta=eta/log(p), seed=seed),
                         n=grid$n.seq, p=grid$p.seq, eta=grid$eta.seq, seed=grid$seed.seq, r=5, snr=2,
                         mc.cores = 1)

        save(res1, file="./res1.RData")
      }

      load("./res1.RData")


      # Reshape the res1 to a matrix form
      n.out <- length(res1) / prod(sapply(param.list,length))
      res1_array <- array(res1,
                             dim=c(n.out, sapply(param.list,length)),
                             dimnames=append(list(out=paste0("V", 1:n.out)), param.list) )

      dimnames(res1_array)

      library(ggsci)
      # eta=0
      # ppca
      ## n=10
      res1_array[1,1,2,1,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)
      ## n=100
      res1_array[1,2,2,1,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)
      ## n=500
      res1_array[1,3,2,1,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)
      
      # gppca
      ## n=10
      res1_array[2,1,2,1,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)
      ## n=100
      res1_array[2,2,2,1,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)
      ## n=500
      res1_array[2,3,2,1,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)
      
      # eta=0
      # ppca
      ## n=10
      res1_array[1,1,2,2,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)
      ## n=100
      res1_array[1,2,2,2,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)
      ## n=500
      res1_array[1,3,2,2,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)
      
      # gppca
      ## n=10
      res1_array[2,1,2,2,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)
      ## n=100
      res1_array[2,2,2,2,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)
      ## n=500
      res1_array[2,3,2,2,1][[1]] %>% filter(!is.na(value)) %>% png.pca.plot_convergence(maxit=100)
      
    }



