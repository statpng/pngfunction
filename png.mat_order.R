png.mat_order <- function(mat, row=TRUE, col=TRUE){
  library(gtools)
  if( row ){
    row.order <- mixedorder(rownames(mat))
  } else {
    row.order <- TRUE
  }
  
  if( col ){
    col.order <- mixedorder(colnames(mat))
  } else {
    col.order <- TRUE
  }
  return( mat[row.order, col.order] )
}
