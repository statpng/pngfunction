TATES <- function(pv, cor_mat, adjustment=TRUE){
      # pv: a q-dimensional vector of p-value
      # cor_mat: a q x q correlation matrix between responses, where q is the number of responses.
      # adjustment: If TRUE, betaT will be used. betaT values were calculated in the following reference.
      # Reference: Sluis, Sophie van der, Danielle Posthuma, and Conor V. Dolan. 2013. “TATES: Efficient Multivariate Genotype-Phenotype Analysis for Genome-Wide Association Studies.” PLOS Genetics 9(1):e1003235.

      pv <- as.numeric(pv)
      pv_sort <- sort(pv, index.return=TRUE)
      pv <- pv_sort$x
      betaT=c(-0.0007609278,-0.0022965148,0.6226249243,0.0148755138,0.1095155903,-0.0218930325,0.2178970393)
      cor2 <- cor_mat[pv_sort$ix, pv_sort$ix]
      
      if( adjustment ){
        for(i in 1:ncol(cor_mat)){
          for(j in 1:i){
            if (i>j) {
              er=cor2[i,j]
              cor2[i,j] <- cor2[j,i] <- betaT[7]*er^6+betaT[6]*er^5+betaT[5]*er^4+betaT[4]*er^3+betaT[3]*er^2+betaT[2]*er+betaT[1]
            }
          }
        }
        diag(cor2) <- 1
      }
      
      
      mej <- NULL
      for( j in 1:length(pv) ){
        eig.values <- eigen(cor2[1:j,1:j])$values
        mej[j] <- j - sum( ifelse( eig.values[1:j]>1, 1, 0 ) * (eig.values[1:j]-1) )
      }
      
      me <- max(mej)
      min(me/mej*pv)
    }
