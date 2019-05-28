Venn2 = function(x,y,duplicated=FALSE){
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

Venn3 = function(x,y,z){
  
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


nonamed = function(x){
  res = x
  names(res) = NULL
  res
}

natochr = function(x){
  x = as.vector(x)
  rep( x[!is.na(x)], diff( c( which( !is.na(x) ), length(x)+1 ) , 1 ) )
}

chrtowh = function(x){
  x = unlist(x)
  ux = unique(x)
  ref = data.frame( ux, 1:length(ux)-1 )
  unlist( sapply( x , function(y) ref[ref[,1]%in%y,2] ) )
}


findchr = function(x, pattern){
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
  res[!is.na(res)] = cut( x, seq( floor(min(x)), ceiling(max(x)), length.out=n+1), right=FALSE)
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



rcolor = function(n, gradient=FALSE){
  
  if( gradient ){
    colfunc <- colorRampPalette(c("white", "red"))
    res = colfunc(n)
  } else {
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    res = col_vector[1:n]
  }
  res
}


