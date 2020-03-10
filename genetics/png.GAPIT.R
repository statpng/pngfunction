## Must be revised!!
if(FALSE){
png.GAPIT <- function(name.x, name.y){
  source("./[CODE]GAPIT_source.R")
  library(dplyr)
  
  if( any( is.na(y) ) ) stop("The process stops because there are missing values in phenotype data (Y).")
  
  
  devtools::source_url("https://raw.githubusercontent.com/statpng/pngfunction/master/png.etc.R")
  
  id.x <- as.character(x[1,-(1:11)] )
  id.y <- as.character(y[,1])
  id.xy <- png.venn2(id.x, id.y)$inner
  X <- x[, c(1:11, 11 + which(x[1,-(1:11)] %in% id.xy))]
  Y <- y %>% subset(ID %in% id.xy)
  
  png.venn2(X[1,-(1:11)] %>% unlist, Y$ID %>% as.character) %>% sapply(length)
  
  #
  
  # source("http://www.zzlab.net/GAPIT/GAPIT.library.R")
  # source("http://www.zzlab.net/GAPIT/gapit_functions.txt")
  
  system(paste0("mkdir ", strsplit(name.x,"_")[[1]][[1]], "_", stringr::str_match(name.y,"[0-9]{1,5}") ))
  gwd <- getwd()
  setwd( paste0("./", paste0(strsplit(name.x,"_")[[1]][[1]], "_", stringr::str_match(name.y,"[0-9]{1,5}")) ) )
  
  
  
  if( !file.exists("./ModelSelection") ) system("mkdir ./ModelSelection")
  setwd("./ModelSelection")
  
  ModelSelection <- GAPIT(Y=Y, 
                          G=X,
                          kinship.algorithm = "VanRaden",
                          # GD=myGD_final,
                          # GM=myGM_final,
                          SNP.MAF=0.01,
                          PCA.total=20,
                          PCA.View.output=FALSE,
                          Model.selection=TRUE)
  
  setwd("..")
  
  
  
  Phenotype.name <- colnames(Y)[2]
  
  ModelSelection <- data.table::fread(paste0("./ModelSelection/GAPIT.MLM.",Phenotype.name,".BIC.Model.Selection.Results.csv"))
  
  PCAtotal <- ModelSelection$"Number of PCs/Covariates"[which.max(ModelSelection$"BIC (larger is better) - Schwarz 1978")]
  
  if( !file.exists("./ECMLM") ) system("mkdir ./ECMLM")
  setwd("./ECMLM")
  
  ECMLM <- GAPIT(Y=Y,
                 G=X,
                 # GD=myGD_final,
                 # GM=myGM_final,
                 SNP.MAF = 0.01,
                 # kinship.cluster = c("average", "complete", "ward.D"),
                 # kinship.group = c("Mean", "Max"),
                 kinship.algorithm = "VanRaden",
                 PCA.total = PCAtotal,
                 PCA.View.output=FALSE,
                 group.from = 10,
                 group.to = 50,
                 group.by = 10,
                 memo = "ECMLM"
  )
  
  setwd("..")
  
  save.image(file=paste0(strsplit(name.x,"_")[[1]][[1]], "_", stringr::str_match(name.y,"[0-9]{1,5}"),".RData"))
  
  setwd(gwd)
  
  
}
}
