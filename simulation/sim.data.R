sim.data <- function(n, p, q, snp.rho, y.rho, maf.min, 
                     scenario = c("scenarioA", "scenarioB", "scenarioC"),
                     response.type = c("continuous", "binary", "mixed"), mu, wh.true = NULL){
  
  
  # n=400; p=40000;
  # snp.rho = 0.95;  maf.min = 0.05
  # y.rho = runif(1, 0.0, 0.2)
  # mu <- c(0.5, 0.75, 1.0, 1.25, 1.5) #1:5*0.1
  
  if( is.null(wh.true) ){
    wh.true <- rep(1:2, gamma/2) + rep( 1:(q/2)-1, each=gamma/q*2 )*2
  }
  
  
  
  
  
  GetBeta <- function(sbeta, p, q, type, true.y, scenario){
    
    beta <- matrix(0, nrow=p, ncol=q)
    if ( scenario == "scenarioB" ){ ## Setting46-2
      for( mm in true.y ){
        if( mm %% 2 == 0 ){
          beta[ 1100+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 1200+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 1300+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 1400+seq_len(length(sbeta)), mm] <- sbeta
          
          beta[ 2100+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 2200+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 2300+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 2400+seq_len(length(sbeta)), mm] <- sbeta
          
          beta[ 3100+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 3200+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 3300+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 3400+seq_len(length(sbeta)), mm] <- sbeta
          
          beta[ 4100+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 4200+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 4300+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 4400+seq_len(length(sbeta)), mm] <- sbeta
        } else {
          beta[ 1100+seq_len(length(sbeta)), mm] <- rev(sbeta)
          beta[ 1200+seq_len(length(sbeta)), mm] <- rev(sbeta)
          beta[ 1300+seq_len(length(sbeta)), mm] <- rev(sbeta)
          beta[ 1400+seq_len(length(sbeta)), mm] <- rev(sbeta)
          
          beta[ 2100+seq_len(length(sbeta)), mm] <- rev(sbeta)
          beta[ 2200+seq_len(length(sbeta)), mm] <- rev(sbeta)
          beta[ 2300+seq_len(length(sbeta)), mm] <- rev(sbeta)
          beta[ 2400+seq_len(length(sbeta)), mm] <- rev(sbeta)
          
          beta[ 3100+seq_len(length(sbeta)), mm] <- rev(sbeta)
          beta[ 3200+seq_len(length(sbeta)), mm] <- rev(sbeta)
          beta[ 3300+seq_len(length(sbeta)), mm] <- rev(sbeta)
          beta[ 3400+seq_len(length(sbeta)), mm] <- rev(sbeta)
          
          beta[ 4100+seq_len(length(sbeta)), mm] <- rev(sbeta)
          beta[ 4200+seq_len(length(sbeta)), mm] <- rev(sbeta)
          beta[ 4300+seq_len(length(sbeta)), mm] <- rev(sbeta)
          beta[ 4400+seq_len(length(sbeta)), mm] <- rev(sbeta)
        }
      }
      
    } else if ( scenario == "scenarioA" ){ ## Setting46-4
      for( mm in true.y ){
        beta[ 1100+seq_len(length(sbeta)), mm] <- sbeta
        beta[ 1200+seq_len(length(sbeta)), mm] <- sbeta
        beta[ 1300+seq_len(length(sbeta)), mm] <- sbeta
        beta[ 1400+seq_len(length(sbeta)), mm] <- sbeta
        
        beta[ 2100+seq_len(length(sbeta)), mm] <- sbeta
        beta[ 2200+seq_len(length(sbeta)), mm] <- sbeta
        beta[ 2300+seq_len(length(sbeta)), mm] <- sbeta
        beta[ 2400+seq_len(length(sbeta)), mm] <- sbeta
        
        beta[ 3100+seq_len(length(sbeta)), mm] <- sbeta
        beta[ 3200+seq_len(length(sbeta)), mm] <- sbeta
        beta[ 3300+seq_len(length(sbeta)), mm] <- sbeta
        beta[ 3400+seq_len(length(sbeta)), mm] <- sbeta
        
        beta[ 4100+seq_len(length(sbeta)), mm] <- sbeta
        beta[ 4200+seq_len(length(sbeta)), mm] <- sbeta
        beta[ 4300+seq_len(length(sbeta)), mm] <- sbeta
        beta[ 4400+seq_len(length(sbeta)), mm] <- sbeta
      }
      
    } else if ( scenario == "scenarioC" ){ ## Setting46-5_sbetaUp
      for( mm in true.y ){
        if( mm %% 2 == 0 ){
          beta[ 1100+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 1200+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 1300+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 1400+seq_len(length(sbeta)), mm] <- sbeta
          
          beta[ 2100+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 2200+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 2300+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 2400+seq_len(length(sbeta)), mm] <- sbeta
          
          beta[ 3100+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 3200+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 3300+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 3400+seq_len(length(sbeta)), mm] <- sbeta
          
          beta[ 4100+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 4200+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 4300+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 4400+seq_len(length(sbeta)), mm] <- sbeta
        } else {
          beta[ 1100+seq_len(length(sbeta)), mm] <- rev(sbeta)
          beta[ 1200+seq_len(length(sbeta)), mm] <- rev(sbeta)
          beta[ 1300+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 1400+seq_len(length(sbeta)), mm] <- sbeta
          
          beta[ 2100+seq_len(length(sbeta)), mm] <- rev(sbeta)
          beta[ 2200+seq_len(length(sbeta)), mm] <- rev(sbeta)
          beta[ 2300+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 2400+seq_len(length(sbeta)), mm] <- sbeta
          
          beta[ 3100+seq_len(length(sbeta)), mm] <- rev(sbeta)
          beta[ 3200+seq_len(length(sbeta)), mm] <- rev(sbeta)
          beta[ 3300+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 3400+seq_len(length(sbeta)), mm] <- sbeta
          
          beta[ 4100+seq_len(length(sbeta)), mm] <- rev(sbeta)
          beta[ 4200+seq_len(length(sbeta)), mm] <- rev(sbeta)
          beta[ 4300+seq_len(length(sbeta)), mm] <- sbeta
          beta[ 4400+seq_len(length(sbeta)), mm] <- sbeta
        }
      }
      
    }
    
    beta
  }
  
  png.varcov <-
    function(p,
             rho = 0,
             type = NULL,
             Omega = 0.001,
             PropOfNeg = 0.25) {
      
      if (!type %in% c("arcov", "random_sign", "all.equal"))
        stop("type should be one among arcov, random_sign, and all.equal")
      
      if (is.null(type))
        stop("'type' should be entered")
      
      if (type == "arcov") {
        out <- outer(1:p, 1:p, function(x, y)
          rho ^ abs(x - y))
        return(out)
      }
      
      if (type == "random_sign") {
        if (PropOfNeg < 0 |
            PropOfNeg > 0.5)
          stop("PropOfNeg must be in [0,0.5].")
        
        
        e.rho <-
          replicate(p * (p - 1) / 2, runif(1, 0, 1) * (-1) ^ rbinom(1, 1, PropOfNeg))
        e.varcov <- matrix(1, p, p)
        e.varcov[upper.tri(e.varcov)] <- e.rho
        e.varcov[lower.tri(e.varcov)] <-
          t(e.varcov)[lower.tri(e.varcov)]
        
        if (isSymmetric(e.varcov)) {
          x <- (e.varcov %*% e.varcov) + diag(Omega, p)
          e.varcov2 <- (x / sqrt(diag(x) %*% t(diag(x))))
          return(e.varcov2)
        } else {
          stop("Error")
        }
      }
      
      if (type == "all.equal") {
        e.varcov2 <- matrix(rho, p, p)
        diag(e.varcov2) <- 1
        
        return(e.varcov2)
      }
      
    }
  
  png.msnp <-
    function(n, p, M = 25, snp.rho = 0.6, y.rho = 0.3,
             maf.min = 0.01, sp.mu, cp.mu, ncp, nsp, ncp.size,
             cov.type = "all.equal", y.p.neg = NULL, omega = NULL, standardization = TRUE) {
      # Depends: png.snp, png.varcov, mnormt
      
      # sp.mu <- rnorm(1, 0.3, 0.05)
      # cp.mu <- rnorm(1, 0.3, 0.05)
      
      PNG.SNP <- png.snp( n = n, p = p, rho = snp.rho, min.maf = maf.min )
      if (is.matrix(y.rho)) {
        VAR <- y.rho
      } else {
        VAR <- png.varcov(M, rho = y.rho, type = cov.type)
      }
      SNP <- PNG.SNP$snp
      dimnames(SNP) <- list(paste0("N", 1:n), paste0("snp", 1:p))
      MAF <- PNG.SNP$MAF
      
      nsig <- ncp + nsp
      
      if (nsp %% M != 0)
        warnings("The number of single-phenotypic variants is not proportional to M!")
      if (ncp > 0)
        var.cp <- (1:nsig)[(1:nsig) %in% seq_len(ncp)]
      if (nsp > 0)
        var.sp <- (1:nsig)[!(1:nsig) %in% seq_len(ncp)]
      
      sbeta <- replicate(M, rep(0, p))
      if (ncp > 0) {
        for (m in 1:ncp.size) {
          sbeta[var.cp * 10, m] <- cp.mu
        }
      }
      if (nsp > 0)  {
        var.sp.list <-
          tapply(var.sp, rep(1:M, each = ceiling(nsp / M))[1:nsp], list)
        for (m in 1:length(var.sp.list)) {
          sbeta[var.sp.list[[m]] * 10, m] <- sp.mu
        }
      }
      
      
      true <- which(sbeta != 0, arr.ind = TRUE)
      if (ncp > 0)  {
        true.cp <- unique(true[duplicated(true[, 1]), 1, drop = T])
      } else {
        true.cp <- NULL
      }
      if (nsp > 0)  {
        if (is.null(true.cp)) {
          true.sp <- true[, 1]
        } else {
          true.sp <- true[!true[, 1] %in% true.cp, 1]
        }
        true.sp <-
          tapply(true.sp , rep(1:M, each = ceiling(nsp / M))[1:nsp], list)
      } else {
        true.sp <- NULL
      }
      
      Y <- SNP %*% sbeta + rmnorm(n = n, varcov = VAR)
      
      ValuesOfArguments =
        list(
          n = n, p = p, M = M, snp.rho = snp.rho, y.rho = y.rho, maf.min = maf.min, 
          sp.mu = sp.mu, cp.mu = cp.mu, ncp = ncp, nsp = nsp,
          cov.type = cov.type
        )
      
      if (standardization)
        SNP <- scale(SNP)
      
      Data <-
        list(
          snp = SNP, y = Y, maf = MAF,
          sbeta = sbeta, true = true, true.cp = true.cp, true.sp = true.sp,
          args = ValuesOfArguments
        )
      return(Data)
    }
  
  
  png.snp <-
    function(n, p, rho = 0.95, min.maf = 0.05, maf.type = c("unif", "sbeta", "chisq")) {
      X <- do.call("cbind", lapply(1:20, function(x)
        mnormt::rmnorm(n, varcov = png.varcov(
          p = (p / 20),
          rho = rho,
          type = "arcov"
        ))))
      Y <- do.call("cbind", lapply(1:20, function(x)
        mnormt::rmnorm(n, varcov = png.varcov(
          p = (p / 20),
          rho = rho,
          type = "arcov"
        ))))
      
      if (maf.type == "unif") {
        MAF <- runif(p, min.maf, 0.5)
        # ref
        # Dai, M., Ming, J., Cai, M., Liu, J., Yang, C., Wan, X., & Xu, Z. (2017). IGESS: a statistical approach to integrating individual-level genotype data and summary statistics in genome-wide association studies. Bioinformatics, 33(18), 2882-2889.
        # Yang, Y., Shi, X., Jiao, Y., Huang, J., Chen, M., Zhou, X., ... & Liu, J. (2020). CoMM-S2: a collaborative mixed model using summary statistics in transcriptome-wide association studies. Bioinformatics, 36(7), 2009-2016.
        # Jiang, W., & Yu, W. (2017). Controlling the joint local false discovery rate is more powerful than meta-analysis methods in joint analysis of summary statistics from multiple genome-wide association studies. Bioinformatics, 33(4), 500-507.
      } else if (maf.type == "sbeta") {
        MAF <-
          rsbeta(p, 0.14, 0.73) %>% {
            ifelse(. < 0.5, ., 1 - .) * (1 - min.maf * 2) + min.maf
          }
        # ref
        # Ionita-Laza, I., Lange, C., & Laird, N. M. (2009). Estimating the number of unseen variants in the human genome. Proceedings of the National Academy of Sciences, 106(13), 5008-5013.
      } else if (maf.type == "chisq") {
        MAF <- truncnorm::rtruncnorm(
          p,
          a = -0.65,
          b = 0.65,
          mean = 0,
          sd = 5
        ) ^ 2
        # MAF <- truncnorm::rtruncnorm(p, a=sqrt(min.maf), b=0.65, mean=0, sd=5)^2
        # No reference
      }
      
      
      for (j in seq_len(p)) {
        perct_X <- rank((X[, j])) / length(X[, j])
        perct_Y <- rank((Y[, j])) / length(Y[, j])
        
        X[perct_X <= MAF[j], j] <- 1
        X[perct_X >  MAF[j], j] <- 0
        Y[perct_Y <= MAF[j], j] <- 1
        Y[perct_Y >  MAF[j], j] <- 0
      }
      
      Data <- rbind(X + Y)
      
      return(list(snp = Data, MAF = MAF))
    }
  
  
  
  
  
  beta <- GetBeta(mu, p, q, scenario, wh.true, scenario = scenario)
  
  tmp = unique( which( beta != 0, arr.ind=TRUE )[,2] )
  cat("Causal variants are associated with", length(tmp), " phenotypes: (", paste0(tmp, collapse=", "), ")\n")
  cat("The effect sizes are (", paste0( sort(unique( beta[which( beta != 0, arr.ind=TRUE )] )), collapse=", "), ") \n")
  
  SNP <- png.snp(n=n, p=p, rho=snp.rho, min.maf=maf.min, maf.type="unif")
  snp <- SNP$snp
  E <- matrix(y.rho, q, q)
  diag(E) <- 1
  error <- rmnorm(n=n, varcov=E)
  Z <- snp %*% beta + error
  Y <- as.data.frame(Z)
  
  if( response.type == "mixed" ){
    for( colcol in 1:n.binary ){
      Y[,colcol] <- as.factor(ifelse( Z[,colcol] > quantile(Z[,colcol], 0.5), 1, 0 ))
      
      # IL <- exp(Z[,colcol]) / (1 + exp(Z[,colcol]))
      # Y[,colcol] <- as.factor( rbinom(length(IL), 1, IL) )
    }
  } else if( response.type == "binary" ){
    for( colcol in 1:ncol(Z) ){
      Y[,colcol] <- as.factor(ifelse( Z[,colcol] > quantile(Z[,colcol], 0.5), 1, 0 ))
      
      # IL <- exp(Z[,colcol]) / (1 + exp(Z[,colcol]))
      # Y[,colcol] <- as.factor( rbinom(length(IL), 1, IL) )
    }
  }
  
  if( response.type == "mixed" ){
    Family <- rep(c("binomial", "gaussian"), c(n.binary, q-n.binary))
  } else if( response.type == "binary" ){
    Family <- rep("binomial", ncol(Y))
  } else if( response.type == "continuous" ){
    Family <- rep("gaussian", ncol(Y))
  }
  
  Data <- list(snp=snp, y=Y, family=Family, maf=SNP$MAF, beta=beta, penalty.factor=NULL)
  
  Data
}



# Data <- sim.data(100, 6000, 8, 0.95, 0.1, 0.05, 
#                  scenario = "scenarioC", response.type = "mixed", 2:6*0.25)

