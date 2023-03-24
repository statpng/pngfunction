# functions

png.mesh3d <- function(mesh, landmark=NULL, radius=0.1){
  library(SurfaceReconstruction)
  library(rgl)
  
  if( !class(mesh)[1] %in% "mesh3d" ){
    mesh <- SurfaceReconstruction::AFSreconstruction(mesh[,1:3] %>% as.matrix)
  }
  
  open3d()
  shade3d(mesh, color = "grey70")
  
  if(!is.null(landmark)){
    landmark %>% spheres3d(col="red", radius=radius)
    # text3d(coord.title,coord.title,coord.title,"mesh.rot", size=20)
  }
  
  
  
}








png.read_plts <- function(ListFiles, nmax=10, center=TRUE, type="triangle"){
  library(data.table)
  
  if(FALSE){
    path <- "/Users/png/Downloads/721.plt"
  }
  
  mesh <- NULL
  for( i in 1:length(ListFiles) ){
    x = ListFiles[i]
    mesh[[i]] <- png.read_plt(x, type=type)
    
    if(center){
      mesh[[i]]$vb[1:3,] <- mesh[[i]]$vb[1:3,] %>% t %>% scale(scale=FALSE) %>% t
    }
    
    
    if(i > nmax) break
  }
  
  mesh
}




png.read_plt <- function(path, type="triangle"){
  library(data.table)
  
  if(FALSE){
    path <- "/Users/png/Downloads/721.plt"
  }
  
  lines <- readLines(path)
  header <- lines[2]
  lines <- lines[-(1:2)]
  
  change_index = as.numeric( stringr::str_extract(header, "[0-9]+") )
  
  first_part_lines <- lines[1:change_index]
  second_part_lines <- lines[change_index:length(lines)]
  
  # 각 부분을 data.table로 변환
  first_part_data <- data.table::fread(text = paste(first_part_lines, collapse = "\n"), header = FALSE, sep = " ") %>% as_tibble()
  second_part_data <- data.table::fread(text = paste(second_part_lines, collapse = "\n"), header = FALSE, sep = " ") %>% as_tibble()
  
  data.table::setnames(first_part_data, c("x", "y", "z", "value"))
  data.table::setnames(second_part_data, c("x", "y", "z"))
  
  if(type == "triangle"){
    mesh <- rgl::tmesh3d(vertices = t(first_part_data[,1:3]), indices = t(second_part_data))
  } else {
    mesh <- rgl::mesh3d(vertices = t(first_part_data[,1:3]), quads = t(second_part_data))
  }
  
  mesh$value <- as.matrix(first_part_data[,4])
  
  mesh
  
}


png.mesh3d.chcol <- function(mesh, newcol=1:3){
  mesh$vb <- mesh$vb[c(newcol,4),]
  mesh$it <- mesh$it[newcol,]
  mesh$normals <- mesh$normals[c(newcol,4),]
  mesh
}


png.3d.rotmat <- function(angle=0){
  # Convert the rotation angle to radians
  rotation_angle_radians <- angle * (pi / 180)
  
  # Create the 3x3 orthogonal matrix for rotation around the Z axis
  Q <- matrix(c(cos(rotation_angle_radians), -sin(rotation_angle_radians), 0,
                sin(rotation_angle_radians), cos(rotation_angle_radians), 0, 0, 0, 1), nrow = 3, ncol = 3, byrow = TRUE)
  
  Q
}


png.mesh3d.rot <- function(mesh, angle=0){
  
  # Convert the rotation angle to radians
  rotation_angle_radians <- angle * (pi / 180)
  
  # Create the 3x3 orthogonal matrix for rotation around the Z axis
  Q <- matrix(c(cos(rotation_angle_radians), -sin(rotation_angle_radians), 0,
                sin(rotation_angle_radians), cos(rotation_angle_radians), 0, 0, 0, 1), nrow = 3, ncol = 3, byrow = TRUE)
  
  
  mesh$vb[1:3,] <- (t(mesh$vb[1:3,])%*%Q) %>% t
  # mesh$it <- mesh$it[newcol,]
  # mesh$normals <- mesh$normals[c(newcol,4),]
  mesh
}


png.plotly3d <- function(df, size=10){
  if(any(class(df) %in% "mesh3d")){
    df <- df$vb[1:3,] %>% t
  }
  
  png.plotly.scatter3d <- function(df, ...){
    library(plotly)  
    colnames(df) <- c("x", "y", "z")
    plot_ly(x=~x, y=~y, z=~z, type = 'scatter3d', data=df %>% as.data.frame(), ...)  
  }
  
  p <- df %>% as.data.frame %>% 
    png.plotly.scatter3d(alpha=1, marker = list(
      symbol = 3,
      size = size,
      line = list(
        color = "grey50",
        width = 2
      )))
  p
}




png.plotly.add_landmark <- function(p, landmark, size=10){
  for( i in 1:nrow(landmark) ){
    x=landmark[i,1]
    y=landmark[i,2]
    z=landmark[i,3]
    
    params <- list( symbol = 3, size = size, 
                    line = list(color="red"))
    p <- plotly::add_markers(p, x,y,z, marker=params)
  }
  p
}




png.landmark <- function(df, landmark){
  if(FALSE){
    landmark
    
  }
  
  data(nose)
  nose1 <- png.mesh3d.chcol(shortnose.mesh, newcol=c(3,1,2))
  nose2 <- longnose.mesh
  
  
  out <- NULL
  for( i in 1:nrow(landmark) ){
    sum_abs <- apply(df, 1, function(x) sum(abs(x-landmark[i,])))
    out[[i]] <- list( min=min(sum_abs), 
                      wh.min=which.min(sum_abs) )
  }
  
  landmark_new <- df[map_int(out, "wh.min"),]
  landmark_new
}




png.compare.icp2 <- function(mesh.org, mesh.target, landmark.org, landmark.target, radius=c(1,1,1)*0.1, seed=1234){
  # Example

  if(FALSE){
    data(nose)
    
    mesh.org <- shortnose.mesh
    mesh.target <- longnose.mesh %>% png.mesh3d.chcol(newcol=c(3,1,2))
    landmark.org <- shortnose.lm
    landmark.target <- longnose.lm[,c(3,1,2)]
    
    
    mesh.org <- mesh[[idx[1]]] %>% png.mesh3d.chcol(c(3,1,2))
    mesh.target <- mesh[[idx[2]]]
    landmark.org <- landmark[[idx[1]]][,c(3,1,2)]
    landmark.target <- landmark[[idx[2]]]
    
    
  }
  
  library(rgl)
  
  n.org <- nrow(landmark.org)
  n.target <- nrow(landmark.target)
  n.min <- min( nrow(landmark.org), nrow(landmark.target) )
  
  if( n.org != n.target ){
    set.seed(seed)
    landmark.org <- landmark.org[ sample(n.org, n.min), ]
    set.seed(seed)
    landmark.target <- landmark.target[ sample(n.target, n.min), ]
  }
  
  icp <- Morpho::icpmat(landmark.org, landmark.target, iterations=1000)
  
  rotmesh <- rotmesh.onto(mesh.org,
                          landmark.org,
                          icp, scale=TRUE, reflection=TRUE)
  
  
  open3d()
  mesh.org %>% shade3d(col="grey90")
  spheres3d(landmark.org, col="red", radius=radius[1])
  # text3d(coord.title,coord.title,coord.title,"mesh.org", size=20)
  open3d()
  mesh.target %>% shade3d(col="blue")
  spheres3d(landmark.target, col="red", radius=radius[2])
  # text3d(coord.title,coord.title,coord.title,"mesh.target", size=20)
  open3d()
  rotmesh$mesh %>% shade3d(col="green")
  rotmesh$yrot %>% spheres3d(col="red", radius=radius[3])
  # text3d(coord.title,coord.title,coord.title,"mesh.rot", size=20)
  
}






png.compare.icp <- function(mesh.org, mesh.target, landmark.org, landmark.target, radius=c(1,1,1), seed=1234){
  # Example
  if(FALSE){
    mesh.org = mesh[[idx[1]]]
    mesh.target = mesh[[idx[2]]]
    landmark.org = landmark[[idx[1]]]
    landmark.target = landmark[[idx[2]]]
    
    # mesh.org, mesh.target, landmark.org, landmark.target
    mesh.org=df3 %>% png.mesh
    mesh.target=df2 %>% png.mesh
    landmark.org=landmark3
    landmark.target=landmark2
    
    data(nose)
    nose1 <- png.mesh3d.chcol(shortnose.mesh, newcol=c(3,1,2))
    nose2 <- longnose.mesh
    
    
    
    require(rgl)
    data(boneData)
    ## rotate, translate and scale the mesh belonging to the first specimen
    ## onto the landmark configuration of the 10th specimen
    shortnose.mesh %>% png.mesh3d.chcol(c(3,1,2)) %>% shade3d(col=2, specular=1)
    longnose.mesh %>% png.mesh3d.chcol(c(3,1,2)) %>% shade3d(col=2, specular=1)
    
    rotmesh <- rotmesh.onto(shortnose.mesh,
                            shortnose.lm,
                            longnose.lm[,c(3,1,2)], scale=TRUE)
    
    
    ## render original mesh
    open3d()
    shade3d(shortnose.mesh, col="grey80", specular=1)
    spheres3d(shortnose.lm, col="red")
    
    
    ## render target mesh
    open3d()
    shade3d(longnose.mesh %>% png.mesh3d.chcol(c(3,1,2)), col="blue", specular=1)
    spheres3d(longnose.lm[,c(3,1,2)], col="red")
    
    
    ## render rotated mesh and landmarks
    open3d()
    shade3d(rotmesh$mesh, col="green", specular=1)
    spheres3d(rotmesh$yrot, col="red")
  }
  
  
  if(FALSE){
    data(nose)
    
    mesh.org <- shortnose.mesh
    mesh.target <- longnose.mesh %>% png.mesh3d.chcol(newcol=c(3,1,2))
    landmark.org <- shortnose.lm
    landmark.target <- longnose.lm[,c(3,1,2)]
  }
  
  library(rgl)
  
  n.org <- nrow(landmark.org)
  n.target <- nrow(landmark.target)
  n.min <- min( nrow(landmark.org), nrow(landmark.target) )
  
  if( n.org != n.target ){
    set.seed(seed)
    landmark.org <- landmark.org[ sample(n.org, n.min), ]
    set.seed(seed)
    landmark.target <- landmark.target[ sample(n.target, n.min), ]
  }
  
  
  rotmesh <- rotmesh.onto(mesh.org,
                          landmark.org,
                          landmark.target, scale=TRUE)
  
  
  open3d()
  mesh.org %>% shade3d(col="grey90")
  spheres3d(landmark.org, col="red", radius=radius[1])
  # text3d(coord.title,coord.title,coord.title,"mesh.org", size=20)
  open3d()
  mesh.target %>% shade3d(col="blue")
  spheres3d(landmark.target, col="red", radius=radius[2])
  # text3d(coord.title,coord.title,coord.title,"mesh.target", size=20)
  open3d()
  rotmesh$mesh %>% shade3d(col="green")
  rotmesh$yrot %>% spheres3d(col="red", radius=radius[3])
  # text3d(coord.title,coord.title,coord.title,"mesh.rot", size=20)
  
}






png.mesh3d.sample <- function(mesh, idx){
  if(FALSE){
    mesh <- mesh3
  }
  
  mesh$vb[,idx]
  
}

