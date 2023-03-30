png.SparseDim <- function(sparsematrix) {
  rows <- sparsematrix@i + 1
  cols <- findInterval(seq(sparsematrix@x) - 1, sparsematrix@p[-1]) + 1
  list(rows = rows, cols = cols)
}
