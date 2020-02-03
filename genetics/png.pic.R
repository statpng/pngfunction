png.pic <- function(maf, method=c("homo", "hetero", "pic")){
  
  # png.pic(c(0.01, 0.12, 0.2, 0.67), "homo") # 0.5034
  # png.pic(c(0.01, 0.12, 0.2, 0.67), "hetero") # 0.4966
  # png.pic(c(0.01, 0.12, 0.2, 0.67), "pic") # 0.446507
  # PIC == 0 -> Marker has only one allele
  # PIC == 1 -> Marker would have an infinite number of allele
  # A PIC greather than 0.7 is considered to be highly informative
  
  switch(method,
         "pic"=sum(
           apply( subset( expand.grid(seq_len(length(maf)), seq_len(length(maf))), Var1!=Var2 ), 1, function(y){
             maf[y[1]]*maf[y[2]] * ( 1 - maf[y[1]]*maf[y[2]] )
           }) ),
         
         "hetero"=sum(
           apply( subset( expand.grid(seq_len(length(maf)), seq_len(length(maf))), Var1!=Var2 ), 1, function(y){
             maf[y[1]]*maf[y[2]]
           }) ),
         
         "homo"=sum(
           apply( subset( expand.grid(seq_len(length(maf)), seq_len(length(maf))), Var1==Var2 ), 1, function(y){
             maf[y[1]]*maf[y[2]]
           }) )
  )
  
}
