
# Inverse Distance Weighting
interpolate_data_IDW <- function(person1_data, person2_data, value_column, grid_resolution = 100) {
  library(sp)
  library(gstat)
  
  radius <- 0.5
  # Create a common grid
  x_range <- seq(-radius, radius, length.out = grid_resolution)
  y_range <- seq(-radius, radius, length.out = grid_resolution)
  grid <- expand.grid(x = x_range, y = y_range)
  
  # Filter the grid to keep only points within the circle x^2 + y^2 <= 1
  grid <- grid[grid$x^2 + grid$y^2 <= radius^2, ]
  
  # Define the spatial data
  coordinates(person1_data) <- ~x + y
  coordinates(person2_data) <- ~x + y
  coordinates(grid) <- ~x + y
  
  # Perform IDW interpolation for each person's data
  idw1 <- idw(as.formula(paste(value_column, "~ 1")), person1_data, grid)
  idw2 <- idw(as.formula(paste(value_column, "~ 1")), person2_data, grid)
  
  # Convert the results to data frames
  idw1_df <- data.frame(coordinates(idw1), value = idw1@data$var1.pred)
  idw2_df <- data.frame(coordinates(idw2), value = idw2@data$var1.pred)
  
  return(list(person1 = idw1_df, person2 = idw2_df))
}


# K-nearest neighbours
interpolate_data_KNN <- function(person1_data, person2_data, value_column, k = 5, grid_resolution = 100) {
  library(sp)
  library(FNN)
  
  radius <- 0.5
  # Create a common grid
  x_range <- seq(-radius, radius, length.out = grid_resolution)
  y_range <- seq(-radius, radius, length.out = grid_resolution)
  grid <- expand.grid(x = x_range, y = y_range)
  
  # Filter the grid to keep only points within the circle x^2 + y^2 <= 1
  grid <- grid[grid$x^2 + grid$y^2 <= radius^2, ]
  
  # Define the spatial data
  coordinates(person1_data) <- ~x + y
  coordinates(person2_data) <- ~x + y
  coordinates(grid) <- ~x + y
  
  # Perform kNN interpolation for each person's data
  knn1 <- knn.reg(train = person1_data@coords[, c("x", "y")], test = grid %>% as.data.frame, y = person1_data[[value_column]], k = k)
  knn2 <- knn.reg(train = person2_data@coords[, c("x", "y")], test = grid %>% as.data.frame, y = person2_data[[value_column]], k = k)
  
  # Convert the results to data frames
  knn1_df <- data.frame(coordinates(grid), value = knn1$pred)
  knn2_df <- data.frame(coordinates(grid), value = knn2$pred)
  
  return(list(person1 = knn1_df, person2 = knn2_df))
}



png_disk_plot <- function(df, col_name = "LAT"){
  
  library(ggplot2)
  data.frame(x=df$x, y=df$y, value=df[col_name]) %>% 
    ggplot() +
    geom_point(aes(x,y,col=value), size=1) +
    # scale_color_gradient2(low="red", high="blue", mid="yellow", na.value="grey90",
    # midpoint=mean(df_value, na.rm=TRUE)) + 
    scale_color_gradientn(name=col_name,
                          colours = ( RColorBrewer::brewer.pal(5,c("Spectral","RdYlBu","RdYlGn")[1]) ), 
                          na.value="grey90") +
    theme_bw()
}
