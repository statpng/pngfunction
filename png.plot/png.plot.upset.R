png.plot.upset <- function(List){
  
  df_total <- cbind.data.frame( ID=List[[1]], " "=TRUE ) %>% {.[!duplicated(.),]}
  
  for( i in 2:length(List) ){
    df <- cbind.data.frame( ID=List[[i]], " "=TRUE ) %>% {.[!duplicated(.),]}
    df_total <- merge(df_total, df, by="ID", all=TRUE)
  }
  
  colnames(df_total) <- c("ID", names(List))
  
  head(df_total)
  
  library(ggplot2)
  library(ComplexUpset)
  
  # pdf(file=paste0("./Figure - SampleOverlap.pdf"), width=9, height=6)
  print( ComplexUpset::upset(df_total, colnames(df_total)[-1], name="", min_size=1, width_ratio=0.15) )
  # + ggtitle("TITLE")
  # dev.off()
  
}
