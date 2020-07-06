library(dplyr)
library(ggplot2)
library(tidyr)
library(magrittr)
library(R.utils)

# save.figure  -----------------------------------------------------------------------

png.save.figure <- function(expr, file="./Figure/figure.jpeg", print=T, save=T, ...){
  
  extension <- strsplit( file, "\\." )[[1]] %>% .[length(.)]
  if(save){
    switch(extension,
           png = {
             png(file, ...)
             eval(parse(text=expr))
             dev.off()
           },
           jpeg = {
             jpeg(file, ...)
             eval(parse(text=expr))
             dev.off()
           },
           pdf = {
             pdf(file, ...)
             eval(parse(text=expr))
             dev.off()
           },
           eps = {
             setEPS()
             postscript(file, ...)
             # postscript("./Figure/Manhattan.eps", width=10, height=5)
             eval(parse(text=expr))
             dev.off()
           },
           
           {stop("filetype not recognized")}
    )
  }
  
  #print?
  if (print) eval(parse(text=expr))
  invisible(NULL)
  
}	


# bundle for pairs  -----------------------------------------------------------------------

png.panel.lm =  function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                      cex = 1, col.smooth = "black", span = 2/3, iter = 3, ...)  {
  reg <- function(x, y, col) abline(lm(y~x), col=col, cex=1.2) 
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) reg(x[ok], y[ok], col.smooth)
}

png.panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  text(0.5, 0.5, txt, cex = 1.7, font = 4)
}

# pairs(iris[,1:3], panel = panel.lm, lower.panel = panel.cor,
#       cex = 1.5, pch = 19, col = iris$Species, # adjustcolor(4, .4), 
#       cex.labels = 2, font.labels = 2)




png.venn2 = function(x,y,duplicated=FALSE){
  inner = intersect(x,y)
  xnyc = x[!x%in%inner]
  ynxc = y[!y%in%inner]
  if(!duplicated){
    inner = unique(inner)
    xnyc = unique(xnyc)
    ynxc = unique(ynxc)
  }
  list( x = xnyc , y = ynxc , inner = inner )
}

png.venn3 = function(x,y,z){
  
  n1 = unique( x )
  n2 = unique( y )
  n3 = unique( z )
  
  n12 = intersect(n1, n2)
  n13 = intersect(n1, n3)
  n23 = intersect(n2, n3)
  
  n123 = intersect( intersect(n12,n13), n23 )
  
  n1 = n1[ ! n1 %in% union( n12, n13) ]
  n2 = n2[ ! n2 %in% union( n12, n23) ]
  n3 = n3[ ! n3 %in% union( n13, n23) ]
  
  n12 = n12[ ! n12 %in% n123 ]
  n13 = n13[ ! n13 %in% n123 ]
  n23 = n23[ ! n23 %in% n123 ]
  
  list( n1=n1, n2=n2, n3=n3, n12=n12, n13=n13, n23=n23, n123=n123 )  
  
}













png.nonamed = function(x){
  res = x
  names(res) = NULL
  res
}

# png.natochr = function(x){
#   x = as.vector(x)
#   rep( x[!is.na(x)], diff( c( which( !is.na(x) ), length(x)+1 ) , 1 ) )
# }
# 
# png.chrtowh = function(x){
#   x = unlist(x)
#   ux = unique(x)
#   ref = data.frame( ux, 1:length(ux)-1 )
#   unlist( sapply( x , function(y) ref[ref[,1]%in%y,2] ) )
# }


png.findchr = function(x, pattern){
  res = x[ grepl(pattern, x) ]
  names(res) = NULL
  res
}


png.replace = function(x, list, b, nonamed=FALSE){
  res = x
  res[list]=b
  if(nonamed) names(res)=NULL
  res
}


png.cut = function(x, n, nonamed=FALSE){
  res = x
  x = as.numeric(na.omit(x))	
  res[!is.na(res)] = cut( x, seq( floor(min(x)), ceiling(max(x)*1.000000001), length.out=n+1), right=FALSE)
  if(nonamed) names(res) = NULL
  res
}


png.math.na = function(x,y,operator,nonamed=FALSE){
  x = as.numeric(x)
  y = as.numeric(y)
  res = data.frame(x,y)
  res = ( apply( res , 1, function(X) ifelse( is.na(X[1])|is.na(X[2]) , NA , eval( parse( text= paste0( X[1], operator, X[2] ) ) ) ) ) )
  if(nonamed) names(res)=NULL
  res
}


		
		
png.order <- function(x){
	library(gtools)
	mixedorder(x)
}
