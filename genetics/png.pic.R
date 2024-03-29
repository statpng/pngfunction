png.pic <- function(allele_frequency, method=c("homo", "hetero", "pic")){
  AF <- allele_frequency
  # png.pic(c(0.01, 0.12, 0.2, 0.67), "homo") # 0.5034
  # png.pic(c(0.01, 0.12, 0.2, 0.67), "hetero") # 0.4966
  # png.pic(c(0.01, 0.12, 0.2, 0.67), "pic") # 0.446507
  # PIC == 0 -> Marker has only one allele
  # PIC == 1 -> Marker would have an infinite number of allele
  # A PIC greather than 0.7 is considered to be highly informative
  
  # c.f. Expected Heterozygosity: 1-sum(c(maf, 1-maf)^2)
  
  switch(method,
         "pic"=sum(
           apply( subset( expand.grid(seq_len(length(AF)), seq_len(length(AF))), Var1!=Var2 ), 1, function(y){
             AF[y[1]]*AF[y[2]] * ( 1 - AF[y[1]]*AF[y[2]] )
           }) ),
         
         "hetero"=sum(
           apply( subset( expand.grid(seq_len(length(AF)), seq_len(length(AF))), Var1!=Var2 ), 1, function(y){
             AF[y[1]]*AF[y[2]]
           }) ),
         
         "homo"=sum(
           apply( subset( expand.grid(seq_len(length(AF)), seq_len(length(AF))), Var1==Var2 ), 1, function(y){
             AF[y[1]]*AF[y[2]]
           }) )
  )
  
}


# install.packages("genetics")
# library(genetics)
# Geno <- genotype(c("AA", "AG", "GG", "AA", "AG", "AG", "AG"), sep="")
# summary(Geno)
