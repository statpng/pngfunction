png.parLapply <- function(cl, X, FUN, ...) {
    library(doSNOW)
    registerDoSNOW(cl)
    pb <- txtProgressBar(max=length(X))
    on.exit(close(pb))
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    
    foreach(i=X, .combine='rbind', .options.snow=opts) %dopar% {
      FUN(i, ...)
    }
  }
