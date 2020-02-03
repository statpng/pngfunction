png.heterozygousCalls <- function(x){
    x %>% replace( ., .=="NN", NA ) %>% unlist %>% { mean( . %in% set.heterozygote ) }
}
