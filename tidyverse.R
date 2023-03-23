library(ggplot2)

capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}
png.labeller <- labeller(
  vore = capitalize,
  # conservation = conservation_status,
  conservation2 = label_wrap_gen(10),
  .default = label_both
)
# usage: facet_grid( Group1~Group2, labeller = png.labeller )


# For library(patchwork)
png.label_plot <- function(label, angle) {
    # usage: png.label_plot("LA", 0)
    ggplot() + geom_text(aes(x = 0, y = 0, label = label), size = 6, fontface = 2, angle=angle) + theme_void()
  }
  
  

png.plot3d <- function(df, theta=135, phi=60, alpha=0.1){
  # df: x, y, z cooridnates
  # library(plot3D)
  # library(ggplot2)
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





png.mesh <- function(df){
  # devtools::install_github("stla/SurfaceReconstruction")
  mesh <- SurfaceReconstruction::AFSreconstruction(df[,1:3] %>% as.matrix)
  mesh
}
png.mesh3d <- function(df){
  library(SurfaceReconstruction)
  library(rgl)
  mesh <- SurfaceReconstruction::AFSreconstruction(df[,1:3] %>% as.matrix)
  open3d()
  wire3d(mesh, color = "grey90")
}
png.plotly3d <- function(df){
  library(plotly)  
  colnames(df) <- c("x", "y", "z")
  plot_ly(x=~x, y=~y, z=~z, type = 'scatter3d', data=df %>% as.data.frame())  
}


# devtools::install_github("zarquon42b/Morpho")
