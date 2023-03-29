# functions

png.data.nose <- function(){
  library(Morpho)
  data(nose)
  mesh <- shortnose.mesh
  landmark <- shortnose.lm
  
  list(mesh=mesh, landmark=landmark)
}

png.mesh3d <- function(df){
  # devtools::install_github("stla/SurfaceReconstruction")
  mesh <- SurfaceReconstruction::AFSreconstruction(df[,1:3] %>% as.matrix)
  mesh
}

png.mesh3d.plot <- function(mesh, landmark = NULL, radius = 0.1) {
  library(SurfaceReconstruction)
  library(rgl)
  
  if (!inherits(mesh, "mesh3d")) {
    mesh <- SurfaceReconstruction::AFSreconstruction(as.matrix(mesh[, 1:3]))
  }
  
  open3d()
  shade3d(mesh, color = "grey70")
  
  if (!is.null(landmark)) {
    spheres3d(landmark, col = "red", radius = radius)
  }
  
}





png.read_plts <- function(file_list, nmax = 10, center = TRUE, type = "triangle") {
  library(data.table)
  
  mesh_list <- list()
  
  for (i in 1:length(file_list)) {
    if (i > nmax) break
    
    file_path <- file_list[i]
    mesh <- png.read_plt(file_path, type = type)
    
    if (center) {
      mesh$vb[1:3,] <- t(scale(t(mesh$vb[1:3,]), scale = FALSE))
    }
    
    mesh_list[[i]] <- mesh
  }
  
  return(mesh_list)
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


png.change_mesh_columns <- function(mesh, newcol = 1:3) {
  mesh$vb <- mesh$vb[c(newcol, 4),]
  mesh$it <- mesh$it[newcol,]
  mesh$normals <- mesh$normals[c(newcol, 4),]
  return(mesh)
}



png.create_rotmat3d <- function(angle = 0) {
  rotation_angle_radians <- angle * (pi / 180)
  Q <- matrix(c(cos(rotation_angle_radians), -sin(rotation_angle_radians), 0,
                sin(rotation_angle_radians), cos(rotation_angle_radians), 0, 0, 0, 1), nrow = 3, ncol = 3, byrow = TRUE)
  return(Q)
}

png.mesh3d.rotate <- function(mesh, angle = 0) {
  Q <- png.create_rotmat3d(angle)
  mesh$vb[1:3,] <- t( t(mesh$vb[1:3,]) %*% Q )
  return(mesh)
}

png.plotly.scatter3d <- function(mesh, size = 2) {
  if (inherits(mesh, "mesh3d")) {
    df <- t(mesh$vb[1:3,])
  } else {
    df <- mesh
  }
  colnames(df) <- c("x", "y", "z")
  
  p <- plot_ly(x = ~x, y = ~y, z = ~z, type = 'scatter3d', data = as.data.frame(df), 
               marker = list(symbol = 3, size = size, color = "grey70",
                             line = list(color="grey10", width=1) ) )
  
  return(p)
}


png.plotly.add_landmark <- function(p, landmark, size = 4) {
  for (i in 1:nrow(landmark)) {
    x <- landmark[i, 1]
    y <- landmark[i, 2]
    z <- landmark[i, 3]
    
    params <- list(symbol = 3, size = size, color = "red")
    p <- plotly::add_markers(p, x, y, z, marker = params)
  }
  
  return(p)
}



png.mesh2df <- function(mesh, type="matrix"){
  if(type=="matrix"){
    out <- as.matrix( t(mesh$vb[1:3,]) )
  } else {
    out <- as.data.frame( t(mesh$vb[1:3,]) )
  }
  colnames(out) <- letters[24:26]
  
  return(out)
}

png.df2mesh <- function(df, type="matrix"){
  
  df <- as.matrix(df[,1:3])
  colnames(df) <- letters[24:26]
  
  mesh <- SurfaceReconstruction::AFSreconstruction( df )
  
  return(mesh)
}



png.align_landmark <- function(mesh, landmark) {
  # example: png.landmark(shortnose.mesh, shortnose.lm)
  
  df <- png.mesh2df(mesh)
  
  out <- NULL
  for (i in 1:nrow(landmark)) {
    sum_abs <- apply(df, 1, function(x) sum(abs(x - landmark[i,])))
    out[[i]] <- list(min = min(sum_abs), wh.min = which.min(sum_abs))
  }
  
  landmark_new <- df[map_int(out, "wh.min"),]
  return(landmark_new)
}




png.compare_icp <- function(mesh.org, mesh.target, landmark.org, landmark.target, use.icp=TRUE, radius=c(1,1,1), iterations=100, seed=1234){
  # Example
  
  library(rgl)
  
  n.org <- nrow(landmark.org)
  n.target <- nrow(landmark.target)
  n.min <- min( n.org, n.target )
  
  if (n.org != n.target) {
    set.seed(seed)
    landmark.org <- landmark.org[sample(n.org, n.min), ]
    set.seed(seed)
    landmark.target <- landmark.target[sample(n.target, n.min), ]
  }
  
  if(use.icp){
    icp <- Morpho::icpmat(landmark.org, landmark.target, iterations = iterations)
  }
  mesh.rot <- rotmesh.onto(mesh.org, landmark.org, icp, scale = TRUE, reflection = TRUE)
  
  open3d()
  mesh.org %>% shade3d(col = "grey90")
  spheres3d(landmark.org, col = "red", radius = radius[1])
  
  open3d()
  mesh.target %>% shade3d(col = "blue")
  spheres3d(landmark.target, col = "red", radius = radius[2])
  
  open3d()
  mesh.rot$mesh %>% shade3d(col = "green")
  mesh.rot$yrot %>% spheres3d(col = "red", radius = radius[3])
  
  return( list(mesh.rot=mesh.rot$mesh, landmark.rot=mesh.rot$yrot, Q=mesh.rot$trafo) )
}



