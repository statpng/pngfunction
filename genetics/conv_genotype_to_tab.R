conv_genotype_to_tab <- function(str){
  str_splitted <- strsplit(str, "\t")
  tab <- table( do.call("c", sapply(str_splitted, function(x) ifelse(is.na(x)|x=="NA", NA, strsplit(x, "")))), useNA = "always" )
  
  print("Allele frequencies")
  print(tab)
  
  print("MAF")
  tab_without_na <- tab[!is.na(names(tab))]
  print(min(tab_without_na)/sum(tab_without_na))
  
  tab
}
