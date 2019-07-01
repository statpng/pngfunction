png.hapmap <- function(x){
  # Initial of function -----------------------------------------------------  
  # 1.Sorting -----------------------------------------------------------------
  ord.x <- gtools::mixedorder(colnames(x[,-(1:11)]))
  
  x <- x[,c(1:11, ord.x+11)]
  
  rm.missing <- which( x[,-(1:11)] %>% apply(1, function(x) mean(is.na(x) | x == "NN")) > 0.2 ) # 2804
  rm.missing.sample <- which( x[,-(1:11)] %>% apply(2, function(x) mean(is.na(x) | x == "NN")) > 0.2 ) # 1
  
  set.heterozygote <- apply(expand.grid(c("T","C","A","G"), c("T","C","A","G")) %>% 
                              filter(Var1 != Var2), 1, paste0, collapse="")
  
  png.heterozygousCalls <- function(x){
    x %>% replace( ., .=="NN", NA ) %>% unlist %>% { mean( . %in% set.heterozygote ) }
  }
  
  rm.hetero <- which( x[,-(1:11)] %>% apply(1, function(x) png.heterozygousCalls(x)) > 0.2 ) # 37
  
  print( x[,-(1:11)][ rm.hetero , ] %>% apply(1, table) )
  
  rm.dup.var <- which(duplicated(x)) # There are no duplicates in variants
  rm.dup.sample <- which(duplicated(t(x[,-(1:11)]))) # There are no duplicates in samples
  
  rm.var <- c(rm.missing, rm.hetero, rm.missing.sample, rm.dup.var)
  rm.sample <- rm.dup.sample
  
  x.removed <- x
  if( length(rm.var)>0 ) x.removed <- x[ -unique( rm.var ), ]
  if( length(rm.sample)>0 ) x.removed <- x[, -(unique( rm.sample )+11) ]
  
  print( dim(x) ) # 49683 x 395
  print( dim(x.removed) ) # 46852 x 395
  
  myX <- as.matrix( rbind( colnames(x.removed), x.removed ) )
  myGD <- apply(myX[-1,-(1:11)], 1,
                function(one) GAPIT.Numericalization(one, bit=2, impute="None", Major.allele.zero=TRUE)) %>% 
    as.data.frame(stringsAsFactors = FALSE)
  myGM <- myX[-1,1:4]
  # myGT <- myX[,c(12:ncol(myX))]

  myX.impute <- rbind( myX[1,], cbind( myX[-1,(1:11)], t( apply(t(myX[-1,-(1:11)]), 2, png.impute.snp) ) ) )
  myGD.impute <- apply(myGD, 2, png.impute.numeric)
                
  Taxa <- as.character( myX[1, -(1:11)] )
  out = list (myX=myX, myGD=data.frame(Taxa=Taxa, myGD), myGD.impute=data.frame(Taxa=Taxa, myGD.impute), myGM=myGM)#, myGT=myGT)
  
  out
}
                

                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
`GAPIT.Numericalization` <-
  function(x,bit=2,effect="Add",impute="None", Create.indicator = FALSE, Major.allele.zero = FALSE, byRow=TRUE){
    #Object: To convert character SNP genotpe to numerical
    #Output: Coresponding numerical value
    #Authors: Feng Tian and Zhiwu Zhang
    # Last update: May 30, 2011
    ##############################################################################################
    if(bit==1)  {
      x[x=="X"]="N"
      x[x=="-"]="N"
      x[x=="+"]="N"
      x[x=="/"]="N"
      x[x=="K"]="Z" #K (for GT genotype)is replaced by Z to ensure heterozygose has the largest value
    }
    
    if(bit==2)  {
      x[x=="XX"]="N"
      x[x=="--"]="N"
      x[x=="++"]="N"
      x[x=="//"]="N"
      x[x=="NN"]="N"
      x[x=="00"]="N"
      
    }
    
    n=length(x)
    lev=levels(as.factor(x))
    lev=setdiff(lev,"N")
    #print(lev)
    len=length(lev)
    #print(len)
    #Jiabo creat this code to convert AT TT to 1 and 2. 2018.5.29
    if(bit==2)
    {
      inter_store=c("AT","AG","AC","TA","GA","CA","GT","TG","GC","CG","CT","TC")
      inter=intersect(lev,inter_store)
      if(length(inter)>1)
      {
        x[x==inter[2]]=inter[1]
        n=length(x)
        lev=levels(as.factor(x))
        lev=setdiff(lev,"N")
        #print(lev)
        len=length(lev)
      }
      if(len==2&bit==2)
      { #inter=intersect(lev,inter_store)
        if(!is.na(inter[1]))
        {
          lev=union(lev,"UU")
          len=len+1
          
        }
      }
      if(len==3&bit==2)
      {
        inter=intersect(lev,inter_store)
      }
      
    }
    #print(lev)
    #print(len)
    #Jiabo code is end here
    
    #Genotype counts
    count=1:len
    for(i in 1:len){
      count[i]=length(x[(x==lev[i])])
    }
    
    #print(count)
    
    if(Major.allele.zero){
      if(len>1 & len<=3){
        #One bit: Make sure that the SNP with the major allele is on the top, and the SNP with the minor allele is on the second position
        if(bit==1){
          count.temp = cbind(count, seq(1:len))
          if(len==3) count.temp = count.temp[-3,]
          count.temp <- count.temp[order(count.temp[,1], decreasing = TRUE),]
          if(len==3)order =  c(count.temp[,2],3)else order = count.temp[,2]
        }
        
        #Two bit: Make sure that the SNP with the major allele is on the top, and the SNP with the minor allele is on the third position
        if(bit==2){
          count.temp = cbind(count, seq(1:len))
          if(len==3) count.temp = count.temp[-2,]
          count.temp <- count.temp[order(count.temp[,1], decreasing = TRUE),]
          if(len==3) order =  c(count.temp[1,2],2,count.temp[2,2])else order = count.temp[,2]
        }
        
        count = count[order]
        lev = lev[order]
        
      }   #End  if(len<=1 | len> 3)
    } #End  if(Major.allele.zero)
    
    #print(x)
    
    #make two  bit order genotype as AA,AT and TT, one bit as A(AA),T(TT) and X(AT)
    if(bit==1 & len==3){
      temp=count[2]
      count[2]=count[3]
      count[3]=temp
    }
    
    #print(lev)
    #print(count)
    position=order(count)
    
    #Jiabo creat this code to convert AT TT to 1 and 2.2018.5.29
    
    lev1=lev
    if(bit==2&len==3)
    {
      lev1[1]=lev[count==sort(count)[1]]
      lev1[2]=lev[count==sort(count)[2]]
      lev1[3]=lev[count==sort(count)[3]]
      position=c(1:3)
      lev=lev1
    }
    #print(lev)
    #print(position)
    #print(inter)
    #Jiabo code is end here
    
    
    #1status other than 2 or 3
    if(len<=1 | len> 3)x=0
    
    #2 status
    if(len==2)x=ifelse(x=="N",NA,ifelse(x==lev[1],0,2))
    
    #3 status
    if(bit==1){
      if(len==3)x=ifelse(x=="N",NA,ifelse(x==lev[1],0,ifelse(x==lev[3],1,2)))
    }else{
      if(len==3)x=ifelse(x=="N",NA,ifelse(x==lev[lev!=inter][1],0,ifelse(x==inter,1,2)))
    }
    
    #print(paste(lev,len,sep=" "))
    #print(position)
    
    #missing data imputation
    if(impute=="Middle") {x[is.na(x)]=1 }
    
    if(len==3){
      if(impute=="Minor")  {x[is.na(x)]=position[1]  -1}
      if(impute=="Major")  {x[is.na(x)]=position[len]-1}
      
    }else{
      if(impute=="Minor")  {x[is.na(x)]=2*(position[1]  -1)}
      if(impute=="Major")  {x[is.na(x)]=2*(position[len]-1)}
    }
    
    #alternative genetic models
    if(effect=="Dom") x=ifelse(x==1,1,0)
    if(effect=="Left") x[x==1]=0
    if(effect=="Right") x[x==1]=2
    
    if(byRow) {
      result=matrix(x,n,1)
    }else{
      result=matrix(x,1,n)
    }
    
    return(result)
  }#end of GAPIT.Numericalization function
