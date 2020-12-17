# image segmentation demo
# author: rwalker (ryan@ryanwalker.us)
# license: MIT

library("jpeg")
library("png")
library("graphics")
library("ggplot2")
library("gridExtra")
#######################################################################################
# This script demonstrates a very simple image segmenter on color scheme
#######################################################################################

png.image.arrange <- function(LIST){
  # helper function to review the results of segmentation visually
  img <- NULL
  for( i in 1:length(LIST) ){
    img[[i]] = rasterGrob(LIST[[i]])
  }
  # img1 = rasterGrob(image.raw)
  # img2 = rasterGrob(image.segmented)
  # grid.arrange(arrangeGrob(img1,img2, nrow=1))
  grid.arrange(arrangeGrob(grobs = img, nrow=1))
}
png.image.print <- function(img){
  rows <- c(0, dim(img)[1]); cols <- c(0, dim(img)[2])
  plot(rows, cols, type = "n", xlab = "", ylab = "")
  rasterImage(img, rows[1], cols[1], rows[2], cols[2])
}
png.image.shrink <- function(input, reduce=10, rotate=FALSE){
  # dim: width, height, depth, channel
  
  # img <- imager::load.image("./data/tiger.png")
  
  input_new <- input
  if( rotate ){
    dim(input_new) <- dim(input)[c(2,1,3,4)]
  } else {
    dim(input_new) <- dim(input)
  }
  for( k in 1:3 ){
    if( rotate ){
      input_new[,,,k] <- t(input[,,,k])
    } else {
      input_new[,,,k] <- input[,,,k]
    }
  }
  
  input_new <- input_new %>% {resize(., round(width(.)/reduce),round(height(.)/reduce))}
  input_new[,,1,1:3]
}

#
# Kmeans based segmenter
# 
segment_image = function(img, n){
  # create a flat, segmented image data set using kmeans
  # Segment an RGB image into n groups based on color values using Kmeans
  df = data.frame(
    red = matrix(img[,,1], ncol=1),
    green = matrix(img[,,2], ncol=1),
    blue = matrix(img[,,3], ncol=1)
  )
  K = kmeans(df,n)
  df$label = K$cluster
  
  # compute rgb values and color codes based on Kmeans centers
  colors = data.frame(
    label = 1:nrow(K$centers), 
    R = K$centers[,"red"],
    G = K$centers[,"green"],
    B = K$centers[,"blue"],
    color=rgb(K$centers)
  )
  
  # merge color codes on to df but maintain the original order of df
  df$order = 1:nrow(df)
  df = merge(df, colors)
  df = df[order(df$order),]
  df$order = NULL
  
  return(df)
  
}




spectral_clustering <- function(X, # matrix of data points
                                nn = 10, # the k nearest neighbors to consider
                                n_eig = 2, # m number of eignenvectors to keep
                                normalized = FALSE) 
{
  
  mutual_knn_graph <- function(X, nn = 10)
  {
    D <- as.matrix( dist(X) ) # matrix of euclidean distances between data points in X
    
    # intialize the knn matrix
    knn_mat <- matrix(0,
                      nrow = nrow(X),
                      ncol = nrow(X))
    
    # find the 10 nearest neighbors for each point
    for (i in 1: nrow(X)) {
      neighbor_index <- order(D[i,])[2:(nn + 1)]
      knn_mat[i,][neighbor_index] <- 1 
    }
    
    # Now we note that i,j are neighbors iff K[i,j] = 1 or K[j,i] = 1 
    knn_mat <- knn_mat + t(knn_mat) # find mutual knn
    
    knn_mat[ knn_mat == 2 ] = 1
    
    return(knn_mat)
  }
  
  graph_laplacian <- function(W, normalized = TRUE)
  {
    stopifnot(nrow(W) == ncol(W)) 
    
    g = colSums(W) # degrees of vertices
    n = nrow(W)
    
    if(normalized)
    {
      D_half = diag(1 / sqrt(g) )
      return( diag(n) - D_half %*% W %*% D_half )
    }
    else
    {
      return( diag(g) - W )
    }
  }
  
  # gglasso( dist(X) )
  # fit.glasso <- glasso::glassopath( as.matrix(dist(X)) ) 
  
  W = mutual_knn_graph(X) # 1. matrix of similarities
  L = graph_laplacian(W, normalized=FALSE) # 2. compute graph laplacian
  ei = eigen(L, symmetric = TRUE) # 3. Compute the eigenvectors and values of L
  n = nrow(L)
  return(ei$vectors[,(n - n_eig):(n - 1)]) # return the eigenvectors of the n_eig smallest eigenvalues
  
}


spectral_segment_image = function(img, n){
  # create a flat, segmented image data set using kmeans
  # Segment an RGB image into n groups based on color values using Kmeans
  
  imgDim <- dim(img)
  
  df = data.frame(
    x = rep(1:imgDim[2], each = imgDim[1]),
    y = rep(imgDim[1]:1, imgDim[2]),
    red = matrix(img[,,1], ncol=1),
    green = matrix(img[,,2], ncol=1),
    blue = matrix(img[,,3], ncol=1)
  )
  

  # # similarity --------------------------------------------------------------
  # 
  # pixdiff=2
  # sigma2=.01 #var(imgvec[,3])
  # simmatrix=matrix(0,nrow(imgvec),nrow(imgvec))
  # for(r in 1:nrow(imgvec)) {
  #   cat(r,"out of",nrow(imgvec),"\n")
  #   simmatrix[r,]=ifelse(abs(imgvec[r,1]-imgvec[,1])<=pixdiff & abs(imgvec[r,2]-imgvec[,2])<=pixdiff,exp(-(imgvec[r,3]-imgvec[,3])^2/sigma2),0)
  # }
  
  fit.spectral <- spectral_clustering(df)  
  K = kmeans(fit.spectral,n)
  df$label = K$cluster
  
  # compute rgb values and color codes based on Kmeans centers
  colors <- NULL
  for( cc in unique(K$cluster) ){
    
    df.tmp <- data.frame(
                label = cc, 
                R = df[K$cluster == cc,"red"] %>% median,
                G = df[K$cluster == cc,"green"] %>% median,
                B = df[K$cluster == cc,"blue"] %>% median
              )
    df.tmp$color = rgb(df.tmp$R, df.tmp$G, df.tmp$B)
    
    colors = rbind.data.frame(colors, df.tmp)
  }
  
  
  # merge color codes on to df but maintain the original order of df
  df$order = 1:nrow(df)
  df = merge(df, colors)
  df = df[order(df$order),]
  df$order = NULL
  
  return(df)
  
}

#
# reconstitue the segmented images to RGB matrix
#
build_segmented_image = function(df, img){
  # reconstitue the segmented images to RGB array
  
  # get mean color channel values for each row of the df.
  R = matrix(df$R, nrow=dim(img)[1])
  G = matrix(df$G, nrow=dim(img)[1])
  B = matrix(df$B, nrow=dim(img)[1])
  
  # reconsitute the segmented image in the same shape as the input image
  img_segmented = array(dim=dim(img))
  img_segmented[,,1] = R
  img_segmented[,,2] = G
  img_segmented[,,3] = B
  
  return(img_segmented)
}

#
# 2D projection for visualizing the kmeans clustering
#
project2D_from_RGB = function(df){
  # Compute the projection of the RGB channels into 2D
  PCA = prcomp(df[,c("red","green","blue")], center=TRUE, scale=TRUE)
  pc2 = PCA$x[,1:2]
  df$x = pc2[,1]
  df$y = pc2[,2]
  return(df[,c("x","y","label","R","G","B", "color")])
}

#
# Create the projection plot of the clustered segments
#
plot_projection <- function(df, sample.size){
  # plot the projection of the segmented image data in 2D, using the
  # mean segment colors as the colors for the points in the projection
  index = sample(1:nrow(df), sample.size)
  return(ggplot(df[index,], aes(x=x, y=y, col=color)) + geom_point(size=2) + scale_color_identity())
}

#
# Inspect
#
inspect_segmentation <- function(image.raw, image.segmented, image.proj, samplesize){
  # helper function to review the results of segmentation visually
  img1 = rasterGrob(image.raw)
  img2 = rasterGrob(image.segmented)
  plt = plot_projection(image.proj, samplesize)
  grid.arrange(arrangeGrob(img1,img2, nrow=1),plt)
}



# tiger -------------------------------------------------------------------

# we can work with both JPEGs and PNGS.  For simplicty, we'll always write out to PNG though.
raw.img <- png.image.shrink( imager::load.image("./data/tiger.png"), reduce=1, rotate=TRUE)
reduced.img <- png.image.shrink( imager::load.image("./data/tiger.png"), reduce=10, rotate=TRUE)
png.image.print(raw.img)
png.image.print(reduced.img)

# segment -- tune the number of segments for each image
kmeans.seg.df = segment_image(reduced.img, 2)
spectral.seg.df = spectral_segment_image(reduced.img, 2)
png.image.arrange(
  list(raw.img, 
       reduced.img, 
       build_segmented_image(kmeans.seg.df, reduced.img)
       )
  )

# write the segmented images to disk
writePNG(img.segmented, "img_segmented.png")

# inspect the results
dev.new()
inspect_segmentation(img, img.segmented, img.proj, 100)



png.image.arrange(raw.img, img, build_segmented_image(seg.df, img))




# HotAirBalloon -----------------------------------------------------------

# we can work with both JPEGs and PNGS.  For simplicty, we'll always write out to PNG though.
raw.img <- png.image.shrink( imager::load.image("./data/HotAirBalloon.png"), reduce=1, rotate=TRUE)
reduced.img <- png.image.shrink( imager::load.image("./data/HotAirBalloon.png"), reduce=20, rotate=TRUE)
png.image.arrange(list(raw.img, reduced.img))

# segment -- tune the number of segments for each image
kmeans.seg.df = segment_image(reduced.img, 2)
spectral.seg.df = spectral_segment_image(reduced.img, 2)


png.image.arrange(
  list(raw.img, 
       reduced.img, 
       build_segmented_image(kmeans.seg.df, reduced.img)
  )
)






# ColorfulBird -----------------------------------------------------------

# we can work with both JPEGs and PNGS.  For simplicty, we'll always write out to PNG though.
raw.img <- png.image.shrink( imager::load.image("./data/ColorfulBird.jpg"), reduce=1, rotate=TRUE)
reduced.img <- png.image.shrink( imager::load.image("./data/ColorfulBird.jpg"), reduce=20, rotate=TRUE)
png.image.arrange(list(raw.img, reduced.img))

# segment -- tune the number of segments for each image
kmeans.seg.df = segment_image(reduced.img, 4)
spectral.seg.df = spectral_segment_image(reduced.img, 4)


png.image.arrange(
  list(raw.img, 
       reduced.img, 
       build_segmented_image(kmeans.seg.df, reduced.img)
  )
)


png.image.arrange(
  list(raw.img, 
       reduced.img, 
       build_segmented_image(spectral.seg.df, reduced.img)
  )
)
