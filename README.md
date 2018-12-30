# pngfunction
bundle of R functions by statpng to analyze the high-dimensional genomic data
 
 
 ## List
 
- **ArrayToLong** : (Replaced by plyr::adply) This function converts an array into a data.frame with long form.
- **btob** : For writing a title with vectors.
- **sp.glmnet** : Calculates the selection probabilities for several families (gaussian, binomial, mgaussian). In the case of mgaussian, 
              type of hypothesis should be selected from "any", "all". "any": A variant is associated with at least one of multiple traits. "all": A variant is associated with all multiple traits.
- **tpr.top** : Calculating the true positive rate for variants with the top N ranked selection probabilities.
- **png.snp** : Generating the simulated SNP with *n* samples and *p* SNPs correlated withing each other(*rho*).
- **th** : Calculating the threshold according to the paper of "Stability selection" (Meinshausen and Buhlmann; 2010).
- **gen.net** : Generating the simulated dataset with *n* samples and *p* variables with *beta*, *cvm*.
- **binom.glmnet.sp** : The case of *family=binomial* in sp.glmnet.
- **gaussian.glmnet.sp** : The case of *family=binomial* in sp.glmnet.
- **png.msnp** : Generating the simulated dataset with *M* traits *p* SNPs that has caucal variants as *ptrue*x*p*.
- **mglmnet** : To be deleted.
- **glmnet_ycov** : To be deleted.
- **mcc** : Calculating the Mattew's Correlation Coefficients.
- **minmax** : Transforming the range of data from 0 to 1.
- **SampleSign** : Randomizing the sign of a matrix at a proportion of *ProbOfNeg*.
- **digits** : Extracting the k-th digits of a number.
- **Varcov_rho** : Generating the variance-covariance matrix with *rho* off-diagonal elements.
- **png.varcov** : Generating the variance-covariance matrix with three types of "equal *rho* cov", "arcov", and "covariance with randomized signs".
