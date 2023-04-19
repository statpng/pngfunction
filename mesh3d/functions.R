# functions

# devtools::install_github("zarquon42b/Morpho")
`%_%` <- function(x, y) {
      paste0(x, y)
}


png.flat.interpolate <- function(nodes, col_names = c("Voltage", "DF", "LAT", "Smax"), method=c("KNN", "IDW")){
  # col_names = c("Voltage", "DF", "LAT", "Smax");  method="KNN"
  
  df_nodes <- df_flat_cleaned$nodes
  
  out <- as.list(1:length(col_names))
  names(out) <- col_names
  for( i in 1:length(col_names) ){
    col <- col_names[i]
    
    df_new <- df_nodes %>% select(x,y,col)
    wh.missing <- which(is.na(df_new[col]))
    
    if( method == "KNN" ){
      fit <- interpolate_data_KNN(df_new[-wh.missing,], df_new[-wh.missing,], col)[[1]]
    } else {
      fit <- interpolate_data_IDW(df_new[-wh.missing,], df_new[-wh.missing,], col, k=2)[[1]]
    }
    colnames(fit)[3] <- col
    if( i == 1 ){
      fit <- fit[,1:3]
    } else {
      fit <- fit[,3]
    }
    
    out[[i]] <- fit
  }
  
  
  out_df <- dplyr::bind_cols(out, .name_repair="minimal")
  
  out_df
}


png.flat.clean <- function(df){
  # df <- df_flat
  nodes <- df$nodes[,1:3]
  edges <- df$edges[,1:3]
  x <- nodes[,1]
  y <- nodes[,2]
  z <- 0 # nodes[,3]
  outliers <- which(x^2 + y^2 + z^2 > 0.5^2+1e-4)
  # nonoutliers <- (1:nrow(nodes))[-outliers]
  nodes_filtered <- df$nodes
  nodes_filtered[outliers,] <- NA
  edges_filtered <- edges %>% filter_all( all_vars(!. %in% outliers) )
  
  edges %>% dim %>% print
  edges_filtered %>% dim %>% print
  
  list(nodes=nodes_filtered, edges=edges_filtered)
}




png.plt2stl <- function(path){
  # path <- "/Volumes/T7/1.Mesh3d/temp1/voltmap.plt/"
  ListFiles <- list.files(path, "plt", full.names = TRUE)
  FileNames <- gsub(".plt", "", map_chr(ListFiles, ~strsplit(.x, "/")[[1]] %>% {.[length(.)]}))
  MESH <- png.read_plts( ListFiles, nmax=Inf, center=TRUE )
  MESH %>% map2( FileNames, ~ Rvcg::vcgStlWrite(.x, filename=paste0(gsub(".plt","",.y)) ) )
}



png_read.vtk <- function(filename, from_python = TRUE){

    if(!file.exists(filename)) stop("Cannot read: ",filename)
  con=file(filename,open='rb',encoding='ASCII')
  on.exit(close(con))
  magic=readLines(con,n=1)
  if(regexpr("# vtk DataFile Version [2345]",magic,ignore.case =T)<0)
    stop("Bad header line in file: ",filename)
  
  title=readLines(con,1)
  encoding=readLines(con,1)
  if(regexpr("ASCII",encoding,ignore.case=TRUE)<0)
    stop("Can only read ASCII encoded VTK pointsets")
  
  datasetLine=toupper(readLines(con,1))
  if(regexpr("^DATASET",datasetLine)<0)
    stop("Missing DATASET line")
  
  datasetType=sub("DATASET\\s+(\\w+)","\\1",datasetLine)
  
  validDatasetTypes<-c("STRUCTURED_POINTS", "STRUCTURED_GRID",
                       "UNSTRUCTURED_GRID", "POLYDATA", "RECTILINEAR_GRID", "FIELD")
  
  if(!datasetType%in%validDatasetTypes)
    stop(datasetType," is not a valid VTK dataset type")
  if(datasetType!="POLYDATA")
    stop("ReadVTKLandmarks can currently only read POLYDATA.",
         " See http://www.vtk.org/VTK/img/file-formats.pdf for details.")
  
  
  
  
  pointsLine=toupper(readLines(con,1))
  if(regexpr("POINTS",pointsLine)<0)
    stop("Missing POINTS definition line")
  ptinfo=unlist(strsplit(pointsLine,"\\s+",perl=TRUE))
  if(length(ptinfo)!=3)
    stop("Unable to extract points information from POINTS line",pointsLine)
  
  
  nummarkers=as.integer(ptinfo[2])
  if(is.na(nummarkers))
    stop("Unable to extract number of points from POINTS line:",pointsLine)
  datatype=ptinfo[3]
  if(!datatype%in%toupper(c("unsigned_char", "char", "unsigned_short", "short", "unsigned_int", "int",
                            "unsigned_long", "long", "float", "double")))
    stop("Unrecognised VTK datatype: ",datatype)
  
  points=scan(con,what=1.0,n=3*nummarkers,quiet=TRUE)
  
  # VTK seems to be hardcoded for 3D
  
  
  m=matrix(points,ncol=3,byrow=T)
  colnames(m)=c("X","Y","Z")
  attr(m,"file")=filename
  attr(m,"title")=title
  attr(m,"vtk_datatype")=datatype
  
  
  
  
  
  toupper(readLines(con,1))
  triangLine=toupper(readLines(con,1))
  if(length(triangLine)==0){
    warning("No data on polygons found")
    return(NULL)
  }
  if(regexpr("POLYGONS",triangLine)<0)
    stop("Missing POLYGONS definition line")
  lninfo=unlist(strsplit(triangLine,"\\s+",perl=TRUE))
  
  triangDataTypeLine=toupper(readLines(con,1))
  lnDataTypeinfo=unlist(strsplit(triangDataTypeLine,"\\s+",perl=TRUE))
  if(length(lninfo)!=3)
    stop("Unable to extract connection information from POLYGONS line:",triangLine)
  
  nummconns=as.integer(lninfo[2])
  if(is.na(nummconns))
    stop("Unable to extract number of connections from POLYGONS line:",triangLine)
  datatype=lnDataTypeinfo[2]
  
  triang=scan(con,what=1.0,n=nummconns,quiet=TRUE)
  triangle_df <- matrix(triang,ncol=3,byrow=T)
  attr(triangle_df,"file")=filename
  attr(triangle_df,"title")=title
  attr(triangle_df,"vtk_datatype")= "int"
  

  
  
  num_edges=as.integer(lninfo[3])
  toupper(readLines(con,1))
  edges = scan(con,what=1.0,n=num_edges,quiet=TRUE)
  
  nodes_df <- matrix(points, ncol=3, byrow=TRUE) %>% as.data.frame
  edges_df <- matrix(edges, ncol=3, byrow=TRUE) %>% as.data.frame
  colnames(nodes_df) <- colnames(edges_df) <- letters[24:26]
  
  if( from_python ){
    edges_df <- edges_df+1
  }
  
  list(nodes=nodes_df, edges=edges_df)
}




png.coord2mesh <- function(nodes, edges, value=NULL, type="triangle"){
  nodes <- nodes[,1:3]
  edges <- edges[,1:3]
      
  data.table::setnames(nodes, c("x", "y", "z"))
  data.table::setnames(edges, c("x", "y", "z"))
  
  if(type == "triangle"){
    mesh <- rgl::tmesh3d(vertices = t(nodes[,1:3]), indices = t(edges))
  } else {
    mesh <- rgl::mesh3d(vertices = t(nodes[,1:3]), quads = t(edges))
  }
  
  if(!is.null(value)){
    mesh$value <- as.matrix(value)
  }
  
  mesh
  
}


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

png.mesh3d.rotate <- function(mesh, angle = 0) {
  Q <- png.create_rotmat3d(angle)
  mesh$vb[1:3,] <- t( t(mesh$vb[1:3,]) %*% Q )
  return(mesh)
}




png.mesh3d.scatter2d <- function(mesh, theta=135, phi=60, alpha=0.1){
  library(plot3D)
  library(ggplot2)
  
  if (inherits(mesh, "mesh3d")) {
    df <- as.matrix(t(mesh$vb[1:3,]))
  }
  colnames(df) <- c("x", "y", "z")
  
  trans_3d_2d <- function(data, theta=135, phi=60) {
    pmat <- plot3D::perspbox(z=diag(2), plot=F, theta=theta, phi=phi)
    
    data <- as.data.frame(data)
    XY <- plot3D::trans3D(
      x = data$x,
      y = data$y,
      z = data$z,
      pmat = pmat) %>%
      data.frame()
    
    data$x <- XY$x
    data$y <- XY$y
    
    return(data[, c('x', 'y')])
  }

  trans_3d_2d(df[,1:3], theta=theta, phi=phi) %>% 
    ggplot(aes(x, y, color=x)) + geom_point(alpha=alpha)  + 
    coord_flip() + 
    theme(legend.position="none") + 
    scale_y_reverse()
  
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




                     
png.cuvia <- function(path, level=1, start=4){
  # start <- 4
  
  for( i in 1:level ){
    if( i > 1 ){
      start <- out$edge_end
    }
    names(start) <- NULL
    rows <- read.csv(path, skip=start, nrows=1, sep=" ", header = FALSE) %>% as.numeric
    
    node_start <- start+2
    node_end <- start+rows[1]+1
    edge_start <- start+rows[1]+2
    edge_end <- start+rows[1]+rows[2]+1
    out <- list(edge_start=edge_start, 
             edge_end=edge_end, 
             node_start=node_start, 
             node_end=node_end)
    
  }
  
  df_edge <- read.csv(path, skip=out$edge_start-1, nrows = with(out, edge_end-edge_start+1), header=F, sep=" " )
  df_node <- read.csv(path, skip=out$node_start-1, nrows = with(out, node_end-node_start+1), header=F, sep=" " )
  
  df_edge[,1:3] <- df_edge[,1:3]+1
  
  list(edge=as_tibble(df_edge), 
       node=as_tibble(df_node),
       row=out)
}

                     
                     
                     
png.cuvia2mesh <- function(fit.cuvia, type="triangle"){
  node <- fit.cuvia$node[,1:4]
  edge <- fit.cuvia$edge[,1:3]
  
  data.table::setnames(node, c("x", "y", "z", "value"))
  data.table::setnames(edge, c("x", "y", "z"))
  
  if(type == "triangle"){
    mesh <- rgl::tmesh3d(vertices = t(node[,1:3]), indices = t(edge))
  } else {
    mesh <- rgl::mesh3d(vertices = t(node[,1:3]), quads = t(edge))
  }
  
  mesh$value <- as.matrix(node[,4])
  
  mesh
  
}

                     
                     
                     
                     
png.cuvia2stl <- function(path, level=1){
  library(dplyr)
  path %>% 
    png.cuvia( level=level ) %>% 
    png.cuvia2mesh() %>% 
    Rvcg::vcgStlWrite(filename= paste0( "cuv2stl_lev", level, "_", strsplit(path,"/")[[1]] %>% {.[length(.)]} %>% gsub(".cuv", "", .) ) )
}

                     
                     
                     
                     
