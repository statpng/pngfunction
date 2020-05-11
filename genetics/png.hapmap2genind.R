png.hapmap2genind <- function(hapmap, ...){
  # hapmap:: (p+1) x (n+11) hapmap-formatted genotype data used in GAPIT
  colnames(df) <- hapmap[1,-(1:11)]
  hapmap <- df2genind(as.data.frame(t(df)), type="codom", ncode=1, ploidy=2, ...)
  hapmap
}
