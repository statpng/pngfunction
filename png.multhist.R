# x = list(runif(100), runif(100))
# beside = TRUE
# freq = FALSE
# probability = !freq
# plot.it = TRUE
# col=c("red", "blue", paste0("gray", floor(seq(1, 90, length.out=5))))
# xlab = expression(paste("LD (",r^2,")"))
# ylab = "Density"
# nclass = NULL
# args = list(col, xlab, ylab)

# breaks <- c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1.0)
# breaks = c(0, 2^(seq(-5, 0)))

png.multhist <- function (x, beside = TRUE, freq = NULL, probability = !freq, 
                          nclass=NULL, breaks=NULL, log2=FALSE, plot.it = TRUE, ...){
  hist.args <- formals(hist.default)
  args <- list(...)
  hargs <- names(args)[names(args) %in% names(hist.args)]
  hist.args[hargs] <- args[hargs]
  hist.args$plot <- FALSE
  hist.args$nclass <- nclass
  if( !is.null(nclass) ){
    allhist <- hist(unlist(x), nclass=hist.args$nclass, plot = FALSE)
  } else if( !is.null(breaks) ) {
    hist.args$breaks <- breaks
    allhist <- hist(unlist(x), hist.args$breaks, plot = FALSE)
  } else {
    allhist <- hist(unlist(x), hist.args$breaks, plot = FALSE)
  }
  
  if (plot.it) {
    barplot.args <- formals(barplot.default)
    bargs <- names(args)[names(args) %in% names(barplot.args)]
    barplot.args[bargs] <- args[bargs]
    barplot.args$beside <- beside
    if ("ann" %in% names(barplot.args)) 
      barplot.args$ann <- eval(barplot.args$ann, envir = barplot.args)
    barplot.args$... <- barplot.args$inside <- NULL
    if (!"names.arg" %in% bargs) 
      barplot.args$names.arg <- signif(allhist$mids, 2)
    if (is.null(freq)) {
      freq <- if (!missing(probability)) 
        !as.logical(probability)
      else TRUE
    }
    comp <- if (freq) {
      "counts"
    } else comp <- "density"
    
    combhist <- t(sapply(x, function(z) hist(z, breaks = allhist$breaks, 
                                             plot = FALSE)[[comp]]))
    
    if( log2 ){
      barplot.args$names.arg <- paste0( "[", log2( allhist$breaks ) %>% {.[-length(.)]}, ", ", log2( allhist$breaks )[-1], ")" )
    }
    
    
    if (plot.it){
      # barplot.args$axisnames = FALSE
      bar <- do.call("barplot", c(list(combhist), barplot.args))
      # barplot.args
      # axis( 1, allhist$breaks )
      # if( !log ){
      #   axis( 1, get(paste0("log", log))(allhist$breaks) )
      # }
      
    }
      
    invisible(list(allhist, combhist))
  }
}
