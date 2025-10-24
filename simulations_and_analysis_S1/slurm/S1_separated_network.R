library(SpatialExperiment)
library(dbscan)
library(alphashape3d)
library(apcluster)
library(plotly)
library(dplyr)
library(reshape2)
library(gtools)
library(cowplot)
library(Hmisc)

### Functions -----


add_spe_metadata3D <- function(spe, metadata, plot_image = TRUE) {
  
  # Ignore the 'background' element in metadata
  metadata[['background']] <- NULL
  
  for (i in seq(length(metadata))) {
    metadata_cluster <- metadata[[i]]
    
    if (!(is.character(metadata_cluster$cluster_type) && length(metadata_cluster$cluster_type) == 1)) {
      stop(paste("cluster_type parameter found in the metadata cluster list", i,"is not a character."))
    }
    
    if (metadata_cluster$cluster_type == "regular") {
      spe <- simulate_clusters3D(spe, list(metadata_cluster), plot_image = F)
    }
    else if (metadata_cluster$cluster_type == "ring") {
      spe <- simulate_rings3D(spe, list(metadata_cluster), plot_image = F)
    }
    else if (metadata_cluster$cluster_type == "double ring") {
      spe <- simulate_double_rings3D(spe, list(metadata_cluster), plot_image = F)
    }
    else {
      stop("cluster_type parameter must be either 'regular', 'ring' or 'double ring'.")
    }
  }
  
  if (plot_image) {
    fig <- plot_cells3D(spe)
    methods::show(fig)
  }
  
  return(spe)
}
plot_cells3D <- function(spe,
                         plot_cell_types = NULL,
                         plot_colours = NULL,
                         feature_colname = "Cell.Type") {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (!is.null(plot_cell_types) && !is.character(plot_cell_types)) {
    stop("`plot_cell_types` is not a character vector or NULL.")
  } 
  if (!is.null(plot_colours) && !is.character(plot_colours)) {
    stop("`plot_colours` is not a character vector or NULL.")
  } 
  if (is.character(plot_colours)) {
    non_colours <- plot_colours[which(!(sapply(plot_colours, function(X) {
      tryCatch(is.matrix(col2rgb(X)), 
               error = function(e) FALSE)
    })))]
    if (length(non_colours) > 0) {
      stop(paste("The following plot_colours are not colours:\n   ",
                 paste(non_colours, collapse = ", ")))
    } 
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste(feature_colname, "is not a valid column in your spe object."))
  }
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(spe), "Cell.Type" = spe[[feature_colname]])
  
  ## If no cell types chosen, use all cell types found in data frame
  if (is.null(plot_cell_types)) {
    warning("plot_cell_types not specified, all cell types found in the spe object will be used.")
    plot_cell_types <- unique(df[["Cell.Type"]])
  }
  ## If no colours inputted, use viridis (D) palette
  if (is.null(plot_colours)) {
    warning("plot_colours not specified, viridis (D) palette will be used.")
    plot_colours <- viridis::viridis(n = length(plot_cell_types), option = "D")
  }
  ## User inputs mismatching cell types and colours
  if (length(plot_cell_types) != length(plot_colours)) {
    stop("Length of plot_cell_types is not equal to length of plot_colours")
  }
  
  ## If cell types have been chosen, check they are found in the spe object
  spe_cell_types <- unique(spe[[feature_colname]])
  unknown_cell_types <- setdiff(plot_cell_types, spe_cell_types)
  
  if (length(unknown_cell_types) == length(plot_cell_types)) {
    stop("None of the plot_cell_types are found in the spe object")
  }
  
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following plot_cell_types are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    plot_colours <- plot_colours[which(plot_cell_types %in% spe_cell_types)]
    plot_cell_types <- intersect(plot_cell_types, spe_cell_types)
  }
  
  ## Factor for feature column
  df[, "Cell.Type"] <- factor(df[, "Cell.Type"],
                              levels = plot_cell_types)
  
  ## Plot
  fig <- plot_ly(df,
                 type = "scatter3d",
                 mode = 'markers',
                 x = ~Cell.X.Position,
                 y = ~Cell.Y.Position,
                 z = ~Cell.Z.Position,
                 color = ~Cell.Type,
                 colors = plot_colours,
                 marker = list(size = 2))
  
  fig <- fig %>% layout(scene = list(xaxis = list(title = 'x', showgrid = T, showaxeslabels = F, showticklabels = T, gridwidth = 5, 
                                                  titlefont = list(size = 20), tickfont = list(size = 15)),
                                     yaxis = list(title = 'y', showgrid = T, showaxeslabels = F, showticklabels = T, gridwidth = 5,
                                                  titlefont = list(size = 20), tickfont = list(size = 15)),
                                     zaxis = list(title = 'z', showgrid = T, showaxeslabels = F, showticklabels = T, gridwidth = 5,
                                                  titlefont = list(size = 20), tickfont = list(size = 15))))
  
  return(fig)
}

simulate_clusters3D <- function(spe,
                                cluster_properties_list,
                                plot_image = TRUE,
                                plot_cell_types = NULL,
                                plot_colours = NULL) {
  
  # Check shape variable of cluster_properties
  shapes <- sapply(cluster_properties_list, function(x) {return(x[["shape"]])})
  n_invalid_shapes <- sum(!(shapes %in% c("sphere", "ellipsoid", "cylinder", "network")))
  if (n_invalid_shapes > 0) {
    stop("`cluster_properties_list` contains invalid shape parameters or no shape parameters.")
  }
  
  for (i in seq(length(cluster_properties_list))) { 
    
    shape <- shapes[[i]]
    
    ### Sphere shape
    if (shape == "sphere") {
      spe <- simulate_sphere_cluster(spe, cluster_properties_list[[i]])
    } 
    
    ### Ellipsoid shape
    if (shape == "ellipsoid") {
      spe <- simulate_ellipsoid_cluster(spe, cluster_properties_list[[i]])
    }
    
    ### Cylinder shape
    if (shape == "cylinder") {
      spe <- simulate_cylinder_cluster(spe, cluster_properties_list[[i]])
    }
    
    ### Network shape
    if (shape == "network") {
      spe <- simulate_network_cluster(spe, cluster_properties_list[[i]])
    }
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe, 
                        plot_cell_types,
                        plot_colours)
    methods::show(fig)
  }
  
  return(spe)
}
simulate_cylinder_cluster <- function(spe, cluster_properties) {
  
  # Check input parameters
  input_parameters <- cluster_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get cylinder properties
  cluster_cell_types <- cluster_properties$cluster_cell_types
  cluster_cell_proportions <- cluster_properties$cluster_cell_proportions
  radius <- cluster_properties$radius
  start_loc <- cluster_properties$start_loc
  end_loc <- cluster_properties$end_loc
  
  ## Check if start and end coordinates of the cylinder are the same
  if (identical(start_loc, end_loc)) warning("Start and end coordinates of the cylinder are the same.")
  
  ## Change cell types in the cylinder cluster
  spe_coords <- spatialCoords(spe)
  
  # Get directional vector
  v1 <- end_loc - start_loc
  
  # Get 'd values of planes' at start_loc and end_loc
  d1 <- sum(v1 * start_loc)
  d2 <- sum(v1 * end_loc)
  
  # Get vector between from each cell to start_loc
  v2 <- sweep(spe_coords, 2, end_loc, "-")
  
  spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                  rowSums(sweep(spe_coords, 2, v1, "*")) >= d1 & rowSums(sweep(spe_coords, 2, v1, "*")) <= d2) & # Cell must be between the planes
                                 (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) <= radius^2), # Cell must be close enough to the cylinder line
                               sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(cluster_properties$cluster_type)) cluster_properties <- append(list(cluster_type = "regular"), cluster_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- cluster_properties
  
  return(spe)
}
simulate_cylinder_dr <- function(spe, dr_properties) {
  
  # Check input parameters
  input_parameters <- dr_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get cylinder double ring properties
  cluster_cell_types <- dr_properties$cluster_cell_types
  cluster_cell_proportions <- dr_properties$cluster_cell_proportions
  radius <- dr_properties$radius
  start_loc <- dr_properties$start_loc
  end_loc <- dr_properties$end_loc
  inner_ring_cell_types <- dr_properties$inner_ring_cell_types
  inner_ring_cell_proportions <- dr_properties$inner_ring_cell_proportions
  inner_ring_width <- dr_properties$inner_ring_width
  outer_ring_cell_types <- dr_properties$outer_ring_cell_types
  outer_ring_cell_proportions <- dr_properties$outer_ring_cell_proportions
  outer_ring_width <- dr_properties$outer_ring_width
  
  ## Check if start and end coordinates of the cylinder are the same
  if (identical(start_loc, end_loc)) warning("Start and end coordinates of the cylinder are the same.")
  
  ## Change cell types in the cylinder cluster
  spe_coords <- spatialCoords(spe)
  
  # Get directional vector
  v1 <- end_loc - start_loc
  
  # Get 'd values of planes' at start_loc and end_loc
  d1 <- sum(v1 * start_loc)
  d2 <- sum(v1 * end_loc)
  
  # Get vector between from each cell to start_loc
  v2 <- sweep(spe_coords, 2, end_loc, "-")
  
  # Start with cells in outer ring
  spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                  rowSums(sweep(spe_coords, 2, v1, "*")) >= d1 & rowSums(sweep(spe_coords, 2, v1, "*")) <= d2) & # Cell must be between the planes
                                 (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) <= (radius + inner_ring_width + outer_ring_width)^2), # Cell must be close enough to the cylinder line
                               sample(outer_ring_cell_types, size = ncol(spe), replace = TRUE, prob = outer_ring_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Start with cells in inner ring
  spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                  rowSums(sweep(spe_coords, 2, v1, "*")) >= d1 & rowSums(sweep(spe_coords, 2, v1, "*")) <= d2) & # Cell must be between the planes
                                 (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) <= (radius + inner_ring_width)^2), # Cell must be close enough to the cylinder line
                               sample(inner_ring_cell_types, size = ncol(spe), replace = TRUE, prob = inner_ring_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Then do cells in the cluster 
  spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                  rowSums(sweep(spe_coords, 2, v1, "*")) >= d1 & rowSums(sweep(spe_coords, 2, v1, "*")) <= d2) & # Cell must be between the planes
                                 (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) <= radius^2), # Cell must be close enough to the cylinder line
                               sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(dr_properties$cluster_type)) dr_properties <- append(list(cluster_type = "double ring"), dr_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- dr_properties
  
  return(spe)
}
simulate_cylinder_ring <- function(spe, ring_properties) {
  
  # Check input parameters
  input_parameters <- ring_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get cylinder ring properties
  cluster_cell_types <- ring_properties$cluster_cell_types
  cluster_cell_proportions <- ring_properties$cluster_cell_proportions
  radius <- ring_properties$radius
  start_loc <- ring_properties$start_loc
  end_loc <- ring_properties$end_loc
  ring_cell_types <- ring_properties$ring_cell_types
  ring_cell_proportions <- ring_properties$ring_cell_proportions
  ring_width <- ring_properties$ring_width
  
  ## Check if start and end coordinates of the cylinder are the same
  if (identical(start_loc, end_loc)) warning("Start and end coordinates of the cylinder are the same.")
  
  ## Change cell types in the cylinder cluster
  spe_coords <- spatialCoords(spe)
  
  # Get directional vector
  v1 <- end_loc - start_loc
  
  # Get 'd values of planes' at start_loc and end_loc
  d1 <- sum(v1 * start_loc)
  d2 <- sum(v1 * end_loc)
  
  # Get vector between from each cell to start_loc
  v2 <- sweep(spe_coords, 2, end_loc, "-")
  
  # Start with cells in ring
  spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                  rowSums(sweep(spe_coords, 2, v1, "*")) >= d1 & rowSums(sweep(spe_coords, 2, v1, "*")) <= d2) & # Cell must be between the planes
                                 (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) <= (radius + ring_width)^2), # Cell must be close enough to the cylinder line
                               sample(ring_cell_types, size = ncol(spe), replace = TRUE, prob = ring_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Then do cells in the cluster 
  spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                  rowSums(sweep(spe_coords, 2, v1, "*")) >= d1 & rowSums(sweep(spe_coords, 2, v1, "*")) <= d2) & # Cell must be between the planes
                                 (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) <= radius^2), # Cell must be close enough to the cylinder line
                               sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                               spe[["Cell.Type"]])
  
  
  # Update current meta data
  if (is.null(ring_properties$cluster_type)) ring_properties <- append(list(cluster_type = "ring"), ring_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- ring_properties
  
  return(spe)
}
simulate_double_rings3D <- function(spe,
                                    dr_properties_list,
                                    plot_image = TRUE,
                                    plot_cell_types = NULL,
                                    plot_colours = NULL) {
  
  # Check shape variable of dr_properties_list
  shapes <- sapply(dr_properties_list, function(x) {return(x[["shape"]])})
  n_invalid_shapes <- sum(!(shapes %in% c("sphere", "ellipsoid", "cylinder", "network")))
  if (n_invalid_shapes > 0) {
    stop("`dr_properties_list` contains invalid shape parameters or no shape parameters.")
  }
  
  for (i in seq(length(dr_properties_list))) { 
    
    shape <- shapes[[i]]
    
    ### Sphere shape with double ring
    if (shape == "sphere") {
      spe <- simulate_sphere_dr(spe, dr_properties_list[[i]])
    } 
    
    ### Ellipsoid shape with double ring
    if (shape == "ellipsoid") {
      spe <- simulate_ellipsoid_dr(spe, dr_properties_list[[i]])
    }
    
    ### Cylinder shape with double ring
    if (shape == "cylinder") {
      spe <- simulate_cylinder_dr(spe, dr_properties_list[[i]])
    }
    
    ### Network shape with double ring
    if (shape == "network") {
      spe <- simulate_network_dr(spe, dr_properties_list[[i]])
    }
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe, 
                        plot_cell_types,
                        plot_colours)
    methods::show(fig)
  }
  
  return(spe)
}
simulate_ellipsoid_cluster <- function(spe, cluster_properties) {
  
  # Check input parameters
  input_parameters <- cluster_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get ellipsoid properties
  cluster_cell_types <- cluster_properties$cluster_cell_types
  cluster_cell_proportions <- cluster_properties$cluster_cell_proportions
  x_radius <- cluster_properties$radii[1]
  y_radius <- cluster_properties$radii[2]
  z_radius <- cluster_properties$radii[3]
  centre_loc <- cluster_properties$centre_loc
  theta <- cluster_properties$axes_rotation[1] * (pi/180) # rotation in x-axis
  alpha <- cluster_properties$axes_rotation[2] * (pi/180) # rotation in y-axis
  beta  <- cluster_properties$axes_rotation[3] * (pi/180) # rotation in z-axis
  
  # Get rotation matrices for rotation in the y-z plane (T2), x-z plane (T3) and x-y plane (T4)
  T1 <- matrix(data = c(1, 0, 0,
                        0, cos(theta), -sin(theta),
                        0, sin(theta), cos(theta)), nrow = 3, ncol = 3, byrow = TRUE)
  T2 <- matrix(data = c(cos(alpha), 0, -sin(alpha),
                        0, 1, 0,
                        sin(alpha), 0, cos(alpha)), nrow = 3, ncol = 3, byrow = TRUE)
  T3 <- matrix(data = c(cos(beta), -sin(beta), 0,
                        sin(beta), cos(beta), 0,
                        0, 0, 1), nrow = 3, ncol = 3, byrow = TRUE)
  
  # Get translation matrix from ellipsoid centre (same as centre...)
  T4 <- centre_loc
  
  ## Change cell types in the ellipsoid cluster
  # Get spatial coords from spe (rows are x, y, z, columns are each cell)
  spe_coords <- t(spatialCoords(spe))
  
  # Apply transformations to spe_coords'
  spe_coords <- solve(T1) %*% solve(T2) %*% solve(T3) %*% (spe_coords - T4)
  x <- spe_coords[1, ]
  y <- spe_coords[2, ]
  z <- spe_coords[3, ]
  
  spe[["Cell.Type"]] <- ifelse((x / x_radius)^2 +
                                 (y / y_radius)^2 +
                                 (z / z_radius)^2 <= 1,
                               sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(cluster_properties$cluster_type)) cluster_properties <- append(list(cluster_type = "regular"), cluster_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- cluster_properties
  
  return(spe)
}

simulate_ellipsoid_dr <- function(spe, dr_properties) {
  
  # Check input parameters
  input_parameters <- dr_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get ellipsoid double ring properties
  cluster_cell_types <- dr_properties$cluster_cell_types
  cluster_cell_proportions <- dr_properties$cluster_cell_proportions
  x_radius <- dr_properties$radii[1]
  y_radius <- dr_properties$radii[2]
  z_radius <- dr_properties$radii[3]
  centre_loc <- dr_properties$centre_loc
  inner_ring_cell_types <- dr_properties$inner_ring_cell_types
  inner_ring_cell_proportions <- dr_properties$inner_ring_cell_proportions
  inner_ring_width <- dr_properties$inner_ring_width
  outer_ring_cell_types <- dr_properties$outer_ring_cell_types
  outer_ring_cell_proportions <- dr_properties$outer_ring_cell_proportions
  outer_ring_width <- dr_properties$outer_ring_width
  theta <- dr_properties$axes_rotation[1] * (pi/180) # rotation in x-axis
  alpha <- dr_properties$axes_rotation[2] * (pi/180) # rotation in y-axis
  beta  <- dr_properties$axes_rotation[3] * (pi/180) # rotation in z-axis
  
  # Get rotation matrices for rotation in the y-z plane (T2), x-z plane (T3) and x-y plane (T4)
  T1 <- matrix(data = c(1, 0, 0,
                        0, cos(theta), -sin(theta),
                        0, sin(theta), cos(theta)), nrow = 3, ncol = 3, byrow = TRUE)
  T2 <- matrix(data = c(cos(alpha), 0, -sin(alpha),
                        0, 1, 0,
                        sin(alpha), 0, cos(alpha)), nrow = 3, ncol = 3, byrow = TRUE)
  T3 <- matrix(data = c(cos(beta), -sin(beta), 0,
                        sin(beta), cos(beta), 0,
                        0, 0, 1), nrow = 3, ncol = 3, byrow = TRUE)
  
  # Get translation matrix from ellipsoid centre (same as centre...)
  T4 <- centre_loc
  
  ## Change cell types in the ellipsoid cluster
  # Get spatial coords from spe (rows are x, y, z, columns are each cell)
  spe_coords <- t(spatialCoords(spe))
  
  # Apply transformations to spe_coords'
  spe_coords <- solve(T1) %*% solve(T2) %*% solve(T3) %*% (spe_coords - T4)
  x <- spe_coords[1, ]
  y <- spe_coords[2, ]
  z <- spe_coords[3, ]
  
  
  # Start with cells in outer ring  
  spe[["Cell.Type"]] <- ifelse((x / (x_radius + inner_ring_width + outer_ring_width))^2 +
                                 (y / (y_radius + inner_ring_width + outer_ring_width))^2 +
                                 (z / (z_radius + inner_ring_width + outer_ring_width))^2 <= 1,
                               sample(outer_ring_cell_types, size = ncol(spe), replace = TRUE, prob = outer_ring_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Then do cells in inner ring  
  spe[["Cell.Type"]] <- ifelse((x / (x_radius + inner_ring_width))^2 +
                                 (y / (y_radius + inner_ring_width))^2 +
                                 (z / (z_radius + inner_ring_width))^2 <= 1,
                               sample(inner_ring_cell_types, size = ncol(spe), replace = TRUE, prob = inner_ring_cell_proportions),
                               spe[["Cell.Type"]])
  
  
  # Then do cells in the cluster  
  spe[["Cell.Type"]] <- ifelse((x / x_radius)^2 +
                                 (y / y_radius)^2 +
                                 (z / z_radius)^2 <= 1,
                               sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(dr_properties$cluster_type)) dr_properties <- append(list(cluster_type = "double ring"), dr_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- dr_properties
  
  return(spe)
}

simulate_ellipsoid_ring <- function(spe, ring_properties) {
  
  # Check input parameters
  input_parameters <- ring_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get ellipsoid ring properties
  cluster_cell_types <- ring_properties$cluster_cell_types
  cluster_cell_proportions <- ring_properties$cluster_cell_proportions
  x_radius <- ring_properties$radii[1]
  y_radius <- ring_properties$radii[2]
  z_radius <- ring_properties$radii[3]
  centre_loc <- ring_properties$centre_loc
  ring_cell_types <- ring_properties$ring_cell_types
  ring_cell_proportions <- ring_properties$ring_cell_proportions
  ring_width <- ring_properties$ring_width
  theta <- ring_properties$axes_rotation[1] * (pi/180) # rotation in x-axis
  alpha <- ring_properties$axes_rotation[2] * (pi/180) # rotation in y-axis
  beta  <- ring_properties$axes_rotation[3] * (pi/180) # rotation in z-axis
  
  # Get rotation matrices for rotation in the y-z plane (T2), x-z plane (T3) and x-y plane (T4)
  T1 <- matrix(data = c(1, 0, 0,
                        0, cos(theta), -sin(theta),
                        0, sin(theta), cos(theta)), nrow = 3, ncol = 3, byrow = TRUE)
  T2 <- matrix(data = c(cos(alpha), 0, -sin(alpha),
                        0, 1, 0,
                        sin(alpha), 0, cos(alpha)), nrow = 3, ncol = 3, byrow = TRUE)
  T3 <- matrix(data = c(cos(beta), -sin(beta), 0,
                        sin(beta), cos(beta), 0,
                        0, 0, 1), nrow = 3, ncol = 3, byrow = TRUE)
  
  # Get translation matrix from ellipsoid centre (same as centre...)
  T4 <- centre_loc
  
  ## Change cell types in the ellipsoid cluster
  # Get spatial coords from spe (rows are x, y, z, columns are each cell)
  spe_coords <- t(spatialCoords(spe))
  
  # Apply transformations to spe_coords'
  spe_coords <- solve(T1) %*% solve(T2) %*% solve(T3) %*% (spe_coords - T4)
  x <- spe_coords[1, ]
  y <- spe_coords[2, ]
  z <- spe_coords[3, ]
  
  # Start with cells in ring  
  spe[["Cell.Type"]] <- ifelse((x / (x_radius + ring_width))^2 +
                                 (y / (y_radius + ring_width))^2 +
                                 (z / (z_radius + ring_width))^2 <= 1,
                               sample(ring_cell_types, size = ncol(spe), replace = TRUE, prob = ring_cell_proportions),
                               spe[["Cell.Type"]])
  
  
  # Then do cells in the cluster  
  spe[["Cell.Type"]] <- ifelse((x / x_radius)^2 +
                                 (y / y_radius)^2 +
                                 (z / z_radius)^2 <= 1,
                               sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                               spe[["Cell.Type"]])
  
  
  # Update current meta data
  if (is.null(ring_properties$cluster_type)) ring_properties <- append(list(cluster_type = "ring"), ring_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- ring_properties
  
  return(spe)
}

simulate_mixing3D <- function(spe,
                              cell_types,
                              cell_proportions,
                              plot_image = TRUE,
                              plot_cell_types = NULL,
                              plot_colours = NULL) {
  
  # Check input parameters
  input_parameters <- list("spe" = spe,
                           "cell_types" = cell_types,
                           "cell_proportions" = cell_proportions,
                           "plot_image" = plot_image)
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Apply mixing
  spe[["Cell.Type"]] <- sample(cell_types, size = ncol(spe), replace = TRUE, prob = cell_proportions)
  
  spe@metadata[["simulation"]][["background"]][["cell_types"]] <- cell_types
  spe@metadata[["simulation"]][["background"]][["cell_proportions"]] <- cell_proportions
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe,
                        plot_cell_types,
                        plot_colours)
    methods::show(fig)
  }
  
  return(spe)
}
simulate_network_cluster <- function(spe, cluster_properties) {  
  
  # Check input parameters
  input_parameters <- cluster_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get network properties
  cluster_cell_types <- cluster_properties$cluster_cell_types
  cluster_cell_proportions <- cluster_properties$cluster_cell_proportions
  n_edges <- cluster_properties$n_edges
  width <- cluster_properties$width
  centre_loc <- cluster_properties$centre_loc
  radius <- cluster_properties$radius
  
  # Number of vertices is always one more than the number of edges for the MST will we make
  n_vertices <- n_edges + 1 
  
  ## Generate n_vertices random points with coords inside a sphere with given radius and centre loc. 
  # Starting with 1000 points inside a cube should be a good enough buffer, unless the user wants more than 1000 edges...
  # Lets stop them from inputting more than 99
  max_edges <- 99
  if (n_edges > max_edges) stop("Only networks with less than 100 edges can be simulated")
  random_coords <- data.frame(x = runif(1000, centre_loc[1] - radius, centre_loc[1] + radius),
                              y = runif(1000, centre_loc[2] - radius, centre_loc[2] + radius),
                              z = runif(1000, centre_loc[3] - radius, centre_loc[3] + radius))
  
  # Then subset points which are inside the sphere
  random_coords <- random_coords[(random_coords$x - centre_loc[1])^2 +
                                   (random_coords$y - centre_loc[2])^2 +
                                   (random_coords$z- centre_loc[3])^2 <= radius^2, ]
  
  ## Subset further and pick 'n_vertices' coords to represent the vertices
  random_coords <- sample_n(random_coords, n_vertices)
  
  ## Get adjacency matrix from points (pairwise distance between points)
  # Assume all points have an edge between each other
  # Assume weight of each edge is equal to the distance between points
  adj_mat <- -1 * apcluster::negDistMat(random_coords)
  
  ## Use prim's algorithm to get edges (i.e. the cells connected by each edge)
  tree_edges <- prims_algorithm(adj_mat)
  
  ### Determine width of cylinders so that cylinders further away are thinner
  tree_edges <- get_tree_depth(tree_edges)
  
  ## Get cluster properties using edge data
  network_cluster_properties <- list()
  max_depth <- max(tree_edges[["depth"]])
  
  for (i in seq(n_edges)) {
    start_loc <- as.numeric(random_coords[tree_edges[i, "vertex1"], ])
    end_loc <- as.numeric(random_coords[tree_edges[i, "vertex2"], ])
    curr_width <- (1 - 0.10 * (max_depth - tree_edges[i, "depth"])) * width # 10% decrease with each depth
    
    # Very unlikely case when width is negative, just ignore these cylinders
    if (curr_width < 0) curr_width <- 0
    
    network_cluster_properties[[i]] <- list(shape = "cylinder",
                                            cluster_cell_types = cluster_cell_types,
                                            cluster_cell_proportions = cluster_cell_proportions,
                                            radius = curr_width,
                                            start_loc = start_loc,
                                            end_loc = end_loc)
  }
  
  network_spe <- simulate_clusters3D(spe,
                                     cluster_properties = network_cluster_properties,
                                     plot_image = F)
  
  # Update current meta data
  metadata <- spe@metadata
  if (is.null(cluster_properties$cluster_type)) cluster_properties <- append(list(cluster_type = "regular"), cluster_properties)
  cluster_properties[["cylinders"]] <- network_cluster_properties # Include metadata of cylinders used to make up network
  metadata[["simulation"]][[paste("cluster", length(metadata[["simulation"]]), sep = "_")]] <- cluster_properties
  
  network_spe@metadata <- metadata
  
  return(network_spe)
}
simulate_network_dr <- function(spe, dr_properties) {  
  
  # Check input parameters
  input_parameters <- dr_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get network double ring properties
  cluster_cell_types <- dr_properties$cluster_cell_types
  cluster_cell_proportions <- dr_properties$cluster_cell_proportions
  n_edges <- dr_properties$n_edges
  width <- dr_properties$width
  centre_loc <- dr_properties$centre_loc
  radius <- dr_properties$radius
  inner_ring_cell_types <- dr_properties$inner_ring_cell_types
  inner_ring_cell_proportions <- dr_properties$inner_ring_cell_proportions
  inner_ring_width <- dr_properties$inner_ring_width
  outer_ring_cell_types <- dr_properties$outer_ring_cell_types
  outer_ring_cell_proportions <- dr_properties$outer_ring_cell_proportions
  outer_ring_width <- dr_properties$outer_ring_width
  
  # Number of vertices is always one more than the number of edges for the MST will we make
  n_vertices <- n_edges + 1 
  
  ## Generate n_vertices random points with coords inside a sphere with given radius and centre loc. 
  # Starting with 1000 points inside a cube should be a good enough buffer, unless the user wants more than 1000 edges...
  # Let's stop them from inputting more than 99
  max_edges <- 99
  if (n_edges > max_edges) stop("Only networks with less than 100 edges can be simulated.")
  random_coords <- data.frame(x = runif(1000, centre_loc[1] - radius, centre_loc[1] + radius),
                              y = runif(1000, centre_loc[2] - radius, centre_loc[2] + radius),
                              z = runif(1000, centre_loc[3] - radius, centre_loc[3] + radius))
  
  # Then subset points which are inside the sphere
  random_coords <- random_coords[(random_coords$x - centre_loc[1])^2 +
                                   (random_coords$y - centre_loc[2])^2 +
                                   (random_coords$z- centre_loc[3])^2 <= radius^2, ]
  
  ## Subset further and pick 'n_vertices' coords to represent the vertices
  random_coords <- sample_n(random_coords, n_vertices)
  
  ## Get adjacency matrix from points (pairwise distance between points)
  # Assume all points have an edge between each other
  # Assume weight of each edge is equal to the distance between points
  adj_mat <- -1 * apcluster::negDistMat(random_coords)
  
  ## Use prim's algorithm to get edges (i.e. the cells connected by each edge)
  tree_edges <- prims_algorithm(adj_mat)
  
  ### Determine width of cylinders so that cylinders further away are thinner
  tree_edges <- get_tree_depth(tree_edges)
  
  ## Get cluster properties using edge data
  network_dr_properties <- list()
  max_depth <- max(tree_edges[["depth"]])
  
  for (i in seq(n_edges)) {
    start_loc <- as.numeric(random_coords[tree_edges[i, "vertex1"], ])
    end_loc <- as.numeric(random_coords[tree_edges[i, "vertex2"], ])
    curr_width <- (1 - 0.10 * (max_depth - tree_edges[i, "depth"])) * width # 10% decrease with each depth
    
    # Very unlikely case when width is negative, just ignore these cylinders
    if (width < 0) {
      width <- 0
    }
    
    network_dr_properties[[i]] <- list(shape = "cylinder",
                                       cluster_cell_types = cluster_cell_types,
                                       cluster_cell_proportions = cluster_cell_proportions,
                                       radius = curr_width,
                                       start_loc = start_loc,
                                       end_loc = end_loc,
                                       inner_ring_cell_types = inner_ring_cell_types,
                                       inner_ring_cell_proportions = inner_ring_cell_proportions,
                                       inner_ring_width = inner_ring_width,
                                       outer_ring_cell_types = outer_ring_cell_types,
                                       outer_ring_cell_proportions = outer_ring_cell_proportions,
                                       outer_ring_width = outer_ring_width)
  }
  
  network_spe <- simulate_double_rings3D(spe,
                                         dr_properties = network_dr_properties,
                                         plot_image = F)
  
  # Update current meta data
  metadata <- spe@metadata
  if (is.null(dr_properties$cluster_type)) dr_properties <- append(list(cluster_type = "double ring"), dr_properties)
  dr_properties[["cylinders"]] <- network_dr_properties # Include metadata of cylinders used to make up network
  metadata[["simulation"]][[paste("cluster", length(metadata[["simulation"]]), sep = "_")]] <- dr_properties
  
  network_spe@metadata <- metadata
  
  return(network_spe)
}
simulate_network_ring <- function(spe, ring_properties) {  
  
  # Check input parameters
  input_parameters <- ring_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get network ring properties
  cluster_cell_types <- ring_properties$cluster_cell_types
  cluster_cell_proportions <- ring_properties$cluster_cell_proportions
  n_edges <- ring_properties$n_edges
  width <- ring_properties$width
  centre_loc <- ring_properties$centre_loc
  radius <- ring_properties$radius
  ring_cell_types <- ring_properties$ring_cell_types
  ring_cell_proportions <- ring_properties$ring_cell_proportions
  ring_width <- ring_properties$ring_width
  
  # Number of vertices is always one more than the number of edges for the MST will we make
  n_vertices <- n_edges + 1 
  
  ## Generate n_vertices random points with coords inside a sphere with given radius and centre loc. 
  # Starting with 1000 points inside a cube should be a good enough buffer, unless the user wants more than 1000 edges...
  # Lets stop them from inputting more than 99
  max_edges <- 99
  if (n_edges > max_edges) stop("Only networks with less than 100 edges can be simulated")
  random_coords <- data.frame(x = runif(1000, centre_loc[1] - radius, centre_loc[1] + radius),
                              y = runif(1000, centre_loc[2] - radius, centre_loc[2] + radius),
                              z = runif(1000, centre_loc[3] - radius, centre_loc[3] + radius))
  
  # Then subset points which are inside the sphere
  random_coords <- random_coords[(random_coords$x - centre_loc[1])^2 +
                                   (random_coords$y - centre_loc[2])^2 +
                                   (random_coords$z- centre_loc[3])^2 <= radius^2, ]
  
  ## Subset further and pick 'n_vertices' coords to represent the vertices
  random_coords <- sample_n(random_coords, n_vertices)
  
  ## Get adjacency matrix from points (pairwise distance between points)
  # Assume all points have an edge between each other
  # Assume weight of each edge is equal to the distance between points
  adj_mat <- -1 * apcluster::negDistMat(random_coords)
  
  ## Use prim's algorithm to get edges (i.e. the cells connected by each edge)
  tree_edges <- prims_algorithm(adj_mat)
  
  ### Determine width of cylinders so that cylinders further away are thinner
  tree_edges <- get_tree_depth(tree_edges)
  
  ## Get cluster properties using edge data
  network_ring_properties <- list()
  max_depth <- max(tree_edges[["depth"]])
  
  for (i in seq(n_edges)) {
    start_loc <- as.numeric(random_coords[tree_edges[i, "vertex1"], ])
    end_loc <- as.numeric(random_coords[tree_edges[i, "vertex2"], ])
    curr_width <- (1 - 0.10 * (max_depth - tree_edges[i, "depth"])) * width # 10% decrease with each depth
    
    # Very unlikely case when width is negative, just ignore these cylinders
    if (width < 0) {
      width <- 0
    }
    
    network_ring_properties[[i]] <- list(shape = "cylinder",
                                         cluster_cell_types = cluster_cell_types,
                                         cluster_cell_proportions = cluster_cell_proportions,
                                         radius = curr_width,
                                         start_loc = start_loc,
                                         end_loc = end_loc,
                                         ring_cell_types = ring_cell_types,
                                         ring_cell_proportions = ring_cell_proportions,
                                         ring_width = ring_width)
  }
  
  network_spe <- simulate_rings3D(spe,
                                  ring_properties = network_ring_properties,
                                  plot_image = F)
  
  # Update current meta data
  metadata <- spe@metadata
  if (is.null(ring_properties$cluster_type)) ring_properties <- append(list(cluster_type = "ring"), ring_properties)
  ring_properties[["cylinders"]] <- network_ring_properties # Include metadata of cylinders used to make up network
  metadata[["simulation"]][[paste("cluster", length(metadata[["simulation"]]), sep = "_")]] <- ring_properties
  
  network_spe@metadata <- metadata
  
  return(network_spe)
}
simulate_ordered_background_cells3D <- function(n_cells, 
                                                length, 
                                                width, 
                                                height,
                                                jitter_proportion = 0.25,
                                                background_cell_type = "Others", 
                                                plot_image = TRUE) {
  
  # Check input parameters
  input_parameters <- list("n_cells" = n_cells,
                           "length" = length,
                           "width" = width,
                           "height" = height,
                           "jitter_proportion" = jitter_proportion,
                           "background_cell_type" = background_cell_type,
                           "plot_image" = plot_image)
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Obtain distance between each point using MAGIC formula
  d_cells <- ((sqrt(2) * length * width * height)/n_cells)^(1/3)
  
  # Get distance between rows, columns and layers using 'd_cells'
  d_rows <- d_cells
  d_cols <- (sqrt(3) / 2) * d_cells
  d_lays <- (sqrt(6) / 3) * d_cells
  
  # Get number of rows, columns and layers
  n_rows <- round(length / d_rows)
  n_cols <- round(width / d_cols)
  n_lays <- round(height / d_lays)
  
  # Step 0. Assume points are on a 3D rectangular grid
  rows <- rep(seq(n_rows), n_cols * n_lays) * d_rows
  cols <- rep(rep(seq(n_cols), each = n_rows), n_lays) * d_cols
  lays <- rep(seq(n_lays), each = n_rows * n_cols) * d_lays
  
  # Step 1. For every odd sheet, every even row shifts by d_cells/2 right
  if (n_cols %% 2 == 0) {
    shift <- rep(c(rep(0, n_rows), rep(d_cells/2, n_rows)), n_cols/2)
  } 
  else {
    shift <- c(rep(c(rep(0, n_rows), rep(d_cells/2, n_rows)), n_cols/2), rep(0, n_rows))
  }
  rows <- rows + c(shift, rep(0, n_rows * n_cols)) # Shift each even row by d_cells/2 right
  
  # Step 2. For every even sheet, odd rows shift d_cells/2 right, all rows shift d_cells/(2*sqrt(3)) up
  if (n_cols %% 2 == 0) {
    shift <- rep(c(rep(d_cells/2, n_rows), rep(0, n_rows)), n_cols/2)
  } 
  else {
    shift <- c(rep(c(rep(d_cells/2, n_rows), rep(0, n_rows)), n_cols/2), rep(d_cells/2, n_rows))
  }
  rows <- rows + c(rep(0, n_rows * n_cols), shift) # Shift each odd row by d_cells/2 right
  cols <- cols + rep(c(0, d_cells/(2 * sqrt(3))), each = n_rows * n_cols) # Shift all rows by d_cells/(2*sqrt(3)) up
  
  # Get total number of cells (should be roughly equal to n_cells)
  n_total <- n_rows * n_cols * n_lays
  
  # Add randomness to the location of the cells
  jitter <- jitter_proportion * d_cells # Jitter is proportional to distance between points in hexagonal grid
  jitter_row <- runif(n_total, -jitter, jitter)
  jitter_col <- runif(n_total, -jitter, jitter)
  jitter_lay <- runif(n_total, -jitter, jitter)
  
  rows <- rows + jitter_row
  cols <- cols + jitter_col
  lays <- lays + jitter_lay
  
  # Put data into data frame
  df <- data.frame("Cell.X.Position" = rows,
                   "Cell.Y.Position" = cols,
                   "Cell.Z.Position" = lays,
                   "Cell.Type" = background_cell_type)
  df$Cell.ID <- paste("Cell", seq(nrow(df)), sep = "_")
  
  # Get metadata
  background_metadata <- list("background_type" = "normal",
                              "n_cells" = n_cells,
                              "length" = length,
                              "width" = width,
                              "height" = height,
                              "jitter_proportion" = jitter_proportion,
                              "cell_types" = background_cell_type,
                              "cell_proportions" = 1)
  simulation_metadata <- list(background = background_metadata)
  
  ## Convert data frame to spe object
  spe <- SpatialExperiment(
    assay = matrix(data = NA, nrow = nrow(df), ncol = nrow(df)),
    colData = df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = list(simulation = simulation_metadata))
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe,
                        background_cell_type,
                        "lightgray")
    methods::show(fig)
  }
  
  return(spe)
}
simulate_random_background_cells3D <- function(n_cells, 
                                               length, 
                                               width, 
                                               height, 
                                               minimum_distance_between_cells,
                                               background_cell_type = "Others", 
                                               plot_image = TRUE) {
  
  # Check input parameters
  input_parameters <- list("n_cells" = n_cells,
                           "length" = length,
                           "width" = width,
                           "height" = height,
                           "minimum_distance_between_cells" = minimum_distance_between_cells,
                           "background_cell_type" = background_cell_type,
                           "plot_image" = plot_image)
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Need to over-sample as cells which are too close will be removed later
  n_cells_inflated <- n_cells * 2
  
  # Use poisson distribution to sample points
  pois_df <- poisson_distribution3D(n_cells = n_cells_inflated, 
                                    length = length, 
                                    width = width, 
                                    height = height)
  
  # Add integer rownames to data frame - each cell is labelled by an integer
  rownames(pois_df) <- seq(nrow(pois_df)) 
  
  ### Check if all other cells are to close to the current cell 
  # Use frNN function: for each point, get all points within min_d of it
  pois_df_distances <- dbscan::frNN(pois_df, 
                                    eps = minimum_distance_between_cells,
                                    query = NULL, 
                                    sort = FALSE)
  
  # For each cell, get all other cells which were within 'minimum_distance_between_cells'
  pois_df_distances_ids <- pois_df_distances$id
  
  # Filter out zero length cells
  pois_df_distances_ids <- Filter(function(x) length(x) != 0, pois_df_distances_ids)
  
  # Get integer labels for the remaining cells
  pois_df_distances_ids_names <- as.integer(names(pois_df_distances_ids))
  
  # Determine which cells should be chosen from pois_df
  cells_chosen <- rep(T, nrow(pois_df))
  for (i in seq_len(length(pois_df_distances_ids))) {
    cells_too_close <- pois_df_distances_ids[[i]]
    
    if (cells_chosen[pois_df_distances_ids_names[i]]) cells_chosen[cells_too_close] <- F
  }
  
  pois_df <- pois_df[cells_chosen, ]
  
  # If number of cells remaining is still higher than n_cells, randomly subset n_cells cells
  if (nrow(pois_df) > n_cells) pois_df <- dplyr::sample_n(pois_df, n_cells)
  
  # Add Cell.Type and Cell.ID
  pois_df$Cell.Type <- background_cell_type
  pois_df$Cell.ID <- paste("Cell", seq(nrow(pois_df)), sep = "_")
  
  # Get metadata
  background_metadata <- list("background_type" = "random",
                              "n_cells" = n_cells,
                              "length" = length,
                              "width" = width,
                              "height" = height,
                              "minimum_distance_between_cells" = minimum_distance_between_cells,
                              "cell_types" = background_cell_type,
                              "cell_proportions" = 1)
  simulation_metadata <- list(background = background_metadata)
  
  ## Convert data frame to spe object
  spe <- SpatialExperiment(
    assay = matrix(data = NA, nrow = nrow(pois_df), ncol = nrow(pois_df)),
    colData = pois_df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = list(simulation = simulation_metadata))
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe,
                        background_cell_type,
                        "lightgray")
    methods::show(fig)
  }
  
  return(spe)
}
simulate_rings3D <- function(spe,
                             ring_properties_list,
                             plot_image = TRUE,
                             plot_cell_types = NULL,
                             plot_colours = NULL) {
  
  # Check shape variable of ring_properties_list
  shapes <- sapply(ring_properties_list, function(x) {return(x[["shape"]])})
  n_invalid_shapes <- sum(!(shapes %in% c("sphere", "ellipsoid", "cylinder", "network")))
  if (n_invalid_shapes > 0) {
    stop("`ring_properties_list` contains invalid shape parameters or no shape parameters.")
  }
  
  for (i in seq(length(ring_properties_list))) { 
    
    shape <- shapes[[i]]
    
    ### Sphere shape with ring
    if (shape == "sphere") {
      spe <- simulate_sphere_ring(spe, ring_properties_list[[i]])
    } 
    
    ### Ellipsoid shape with ring
    else if (shape == "ellipsoid") {
      spe <- simulate_ellipsoid_ring(spe, ring_properties_list[[i]])
    }
    
    ### Cylinder shape with ring
    else if (shape == "cylinder") {
      spe <- simulate_cylinder_ring(spe, ring_properties_list[[i]])
    }
    
    ### Network shape with ring
    else if (shape == "network") {
      spe <- simulate_network_ring(spe, ring_properties_list[[i]])
    }
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe, 
                        plot_cell_types,
                        plot_colours)
    methods::show(fig)
  }
  
  return(spe)
}
simulate_spe_metadata3D <- function(spe_metadata, 
                                    plot_image = TRUE, 
                                    plot_cell_types = NULL,
                                    plot_colours = NULL) {
  
  # First element should contain background metadata
  bg_metadata <- spe_metadata[[1]]
  if (!(is.character(bg_metadata$background_type) && length(bg_metadata$background_type) == 1)) {
    stop("background_type parameter found in the metadata background list is not a character.")
  }
  
  if (bg_metadata$background_type == "random") {
    spe <- simulate_random_background_cells3D(bg_metadata$n_cells,
                                              bg_metadata$length,
                                              bg_metadata$width,
                                              bg_metadata$height,
                                              bg_metadata$minimum_distance_between_cells,
                                              plot_image = F)    
  }
  else if (bg_metadata$background_type == "ordered") {
    spe <- simulate_ordered_background_cells3D(bg_metadata$n_cells,
                                               bg_metadata$length,
                                               bg_metadata$width,
                                               bg_metadata$height,
                                               bg_metadata$jitter_proportion,
                                               plot_image = F) 
  }
  else {
    stop("background_type parameter found in the first list must be 'random' or 'ordered'.")
  }
  # Apply background mixing
  spe <- simulate_mixing3D(spe,
                           bg_metadata$cell_types,
                           bg_metadata$cell_proportions,
                           plot_image = F)
  
  ### If there is only background metadata, we are done
  if (length(spe_metadata) == 1) {
    
    # Plot
    if (plot_image) {
      fig <- plot_cells3D(spe,
                          plot_cell_types,
                          plot_colours)
      methods::show(fig)
    }
    
    return(spe)
  }
  
  ### All other elements should help to simulate clusters 
  for (i in 2:length(spe_metadata)) {
    cluster_metadata <- spe_metadata[[i]]
    
    if (!(is.character(cluster_metadata$cluster_type) && length(cluster_metadata$cluster_type) == 1)) {
      stop(paste("cluster_type parameter found in the metadata cluster list", i,"is not a character."))
    }
    
    if (cluster_metadata$cluster_type == "regular") {
      spe <- simulate_clusters3D(spe, list(cluster_metadata), plot_image = F)
    }
    else if (cluster_metadata$cluster_type == "ring") {
      spe <- simulate_rings3D(spe, list(cluster_metadata), plot_image = F)      
    }
    else if (cluster_metadata$cluster_type == "double ring") {
      spe <- simulate_double_rings3D(spe, list(cluster_metadata), plot_image = F)
    }
    else {
      stop("cluster_type parameter must be either 'regular', 'ring' or 'double ring'.")
    }
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe,
                        plot_cell_types,
                        plot_colours)
    methods::show(fig)
  }
  
  return(spe)
}

simulate_sphere_cluster <- function(spe, cluster_properties) {
  
  # Check input parameters
  input_parameters <- cluster_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get sphere properties
  cluster_cell_types <- cluster_properties$cluster_cell_types
  cluster_cell_proportions <- cluster_properties$cluster_cell_proportions
  radius <- cluster_properties$radius
  centre_loc <- cluster_properties$centre_loc
  
  # Change cell types in the sphere cluster
  spe_coords <- data.frame(spatialCoords(spe))
  
  spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                 (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                 (spe_coords$Cell.Z.Position - centre_loc[3])^2 <= radius^2,
                               sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(cluster_properties$cluster_type)) cluster_properties <- append(list(cluster_type = "regular"), cluster_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- cluster_properties
  
  return(spe)
}
simulate_sphere_dr <- function(spe, dr_properties) {
  
  # Check input parameters
  input_parameters <- dr_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get sphere double ring properties
  cluster_cell_types <- dr_properties$cluster_cell_types
  cluster_cell_proportions <- dr_properties$cluster_cell_proportions
  radius <- dr_properties$radius
  centre_loc <- dr_properties$centre_loc
  inner_ring_cell_types <- dr_properties$inner_ring_cell_types
  inner_ring_cell_proportions <- dr_properties$inner_ring_cell_proportions
  inner_ring_width <- dr_properties$inner_ring_width
  outer_ring_cell_types <- dr_properties$outer_ring_cell_types
  outer_ring_cell_proportions <- dr_properties$outer_ring_cell_proportions
  outer_ring_width <- dr_properties$outer_ring_width
  
  ## Change cell types in the sphere ringed cluster
  spe_coords <- data.frame(spatialCoords(spe))
  
  # Start with cells in outer ring  
  spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                 (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                 (spe_coords$Cell.Z.Position - centre_loc[3])^2 <= (radius + inner_ring_width + outer_ring_width)^2,
                               sample(outer_ring_cell_types, size = ncol(spe), replace = TRUE, prob = outer_ring_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Then do cells in inner ring  
  spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                 (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                 (spe_coords$Cell.Z.Position - centre_loc[3])^2 <= (radius + inner_ring_width)^2,
                               sample(inner_ring_cell_types, size = ncol(spe), replace = TRUE, prob = inner_ring_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Then do cells in the cluster 
  spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                 (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                 (spe_coords$Cell.Z.Position - centre_loc[3])^2 <= radius^2,
                               sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(dr_properties$cluster_type)) dr_properties <- append(list(cluster_type = "double ring"), dr_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- dr_properties
  
  return(spe)
}
simulate_sphere_ring <- function(spe, ring_properties) {
  
  # Check input parameters
  input_parameters <- ring_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get sphere ring properties
  cluster_cell_types <- ring_properties$cluster_cell_types
  cluster_cell_proportions <- ring_properties$cluster_cell_proportions
  radius <- ring_properties$radius
  centre_loc <- ring_properties$centre_loc
  ring_cell_types <- ring_properties$ring_cell_types
  ring_cell_proportions <- ring_properties$ring_cell_proportions
  ring_width <- ring_properties$ring_width
  
  ## Change cell types in the sphere ringed cluster
  spe_coords <- data.frame(spatialCoords(spe))
  
  # Start with cells in ring  
  spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                 (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                 (spe_coords$Cell.Z.Position - centre_loc[3])^2 <= (radius + ring_width)^2,
                               sample(ring_cell_types, size = ncol(spe), replace = TRUE, prob = ring_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Then do cells in the cluster 
  spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                 (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                 (spe_coords$Cell.Z.Position - centre_loc[3])^2 <= radius^2,
                               sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(ring_properties$cluster_type)) ring_properties <- append(list(cluster_type = "ring"), ring_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- ring_properties
  
  return(spe)
}
spaSim3D_background_integrator <- function() {
  
  ### Message strings
  message_background <- paste("Hello spaSim-3D user, how do you want your background cells to look like?\n
          1. Random pattern\n
          2. Ordered pattern\n\n",
                              "In a random pattern, cells are placed randomly...\n",
                              "In a ordered pattern, cells follow a regularly spaced in a hexagonal grid\n",
                              "To choose, please enter 1 or 2.\n", sep = "")
  
  message_background_random <- paste("We will need a few parameters before we can obtain the simulation\n",
                                     "    Window size - length, width and height (e.g. 100 x 100 x 100)\n",
                                     "    Number of cells (e.g. 10000 cells)\n",
                                     "    Minimum distance between cells (e.g. minimum distance of 2)\n",
                                     "If you want to change your inputs, you'll be able to at the end.\n", sep = "")
  
  message_background_ordered <- paste("We will need a few parameters before we can obtain the simulation\n",
                                      "    Window size - length, width and height (e.g. 100 x 100 x 100)\n",
                                      "    Number of cells (e.g. 10000 cells)\n",
                                      "    Amount of jitter (choose to give a bit or a lot of randomness)\n",
                                      "If you want to change your inputs, you'll be able to at the end.\n", sep = "")
  
  message_mixing <- paste("Would you like to MIX the background cells with chosen cell types randomly?\n")
  
  # Ask if user wants a 'random' or 'ordered' patterned background
  message(message_background)
  user_input_background <- get_integer_input_from_options(c(1, 2))
  
  ### Simulate random pattern
  if (user_input_background == 1) {
    
    # Get required parameters for a random background from user
    message(message_background_random)
    parameter_values <- list("length" = get_positive_numeric_input("length"),
                             "width" = get_positive_numeric_input("width"),
                             "height" = get_positive_numeric_input("height"),
                             "number of cells" = get_positive_numeric_input("number of cells"),
                             "minimum distance between cells" = get_non_negative_numeric_input("minimum distance between cells"))
    display_parameters(parameter_values)
    
    # Generate random background simulation using these parameters
    message("Generating simulation...")
    simulated_spe <- simulate_random_background_cells3D(parameter_values[["number of cells"]],
                                                        parameter_values[["length"]],
                                                        parameter_values[["width"]],
                                                        parameter_values[["height"]],
                                                        parameter_values[["minimum distance between cells"]])
    
    # Allow user the option to change their input parameters
    message("Would you like to change your input parameters?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      user_input_parameter_choice <- get_integer_input_from_options(seq(length(parameter_values)))
      
      if (user_input_parameter_choice == 1) parameter_values[["length"]] <- get_positive_numeric_input("length")
      if (user_input_parameter_choice == 2) parameter_values[["width"]] <- get_positive_numeric_input("width")
      if (user_input_parameter_choice == 3) parameter_values[["height"]] <- get_positive_numeric_input("height")
      if (user_input_parameter_choice == 4) parameter_values[["number of cells"]] <- get_positive_numeric_input("number of cells")
      if (user_input_parameter_choice == 5) parameter_values[["minimum distance between cells"]] <- get_non_negative_numeric_input("minimum distance between cells")
      
      # Generate random background simulation using updated parameters
      display_parameters(parameter_values)
      message("Generating simulation...")
      simulated_spe <- simulate_random_background_cells3D(parameter_values[["number of cells"]],
                                                          parameter_values[["length"]],
                                                          parameter_values[["width"]],
                                                          parameter_values[["height"]],
                                                          parameter_values[["minimum distance between cells"]])
      
      message("Would you like to change your inputs?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
  }
  ### Simulate ordered pattern
  else if (user_input_background == 2) {
    
    # Get required parameters for a ordered background from user
    message(message_background_ordered)
    parameter_values <- list("length" = get_positive_numeric_input("length"),
                             "width" = get_positive_numeric_input("width"),
                             "height" = get_positive_numeric_input("height"),
                             "number of cells" = get_positive_numeric_input("number of cells"),
                             "amount of jitter" = get_numeric_between_input("amount of jitter", 0, 1))
    display_parameters(parameter_values)
    
    # Generate ordered background simulation using these parameters
    message("Generating simulation...")
    simulated_spe <- simulate_ordered_background_cells3D(parameter_values[["number of cells"]],
                                                         parameter_values[["length"]],
                                                         parameter_values[["width"]],
                                                         parameter_values[["height"]],
                                                         parameter_values[["amount of jitter"]])
    
    # Allow user the option to change their input parameters
    message("Would you like to change your inputs?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      user_input_parameter_choice <- get_integer_input_from_options(seq(length(parameter_values))) # 5 different parameters
      
      if (user_input_parameter_choice == 1) parameter_values[["length"]] <- get_positive_numeric_input("length")
      if (user_input_parameter_choice == 2) parameter_values[["width"]] <- get_positive_numeric_input("width")
      if (user_input_parameter_choice == 3) parameter_values[["height"]] <- get_positive_numeric_input("height")
      if (user_input_parameter_choice == 4) parameter_values[["number of cells"]] <- get_positive_numeric_input("number of cells")
      if (user_input_parameter_choice == 5) parameter_values[["amount of jitter"]] <- get_numeric_between_input("amount of jitter", 0, 1)
      
      # Generate ordered background simulation using updated parameters
      display_parameters(parameter_values)
      message("Generating simulation...")
      simulated_spe <- simulate_ordered_background_cells3D(parameter_values[["number of cells"]],
                                                           parameter_values[["length"]],
                                                           parameter_values[["width"]],
                                                           parameter_values[["height"]],
                                                           parameter_values[["amount of jitter"]])
      message("Would you like to change your inputs?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
  }
  
  ### Simulate mixing
  message(message_mixing)
  choose_cell_types_y_or_n <- get_y_or_n_input()
  if (choose_cell_types_y_or_n == "y") {
    simulated_spe <- get_cell_types_and_proportions_for_mixing(simulated_spe) 
  }
  message("All done!")
  
  return(simulated_spe)
}
spaSim3D_cluster_integrator <- function(simulated_spe = NULL) {
  
  ### Message strings
  message_no_simulated_spe <- paste("Hello spaSim-3D user. Please input your simulated spe object into this function.\n",
                                    "If you don't have any, you can use the spaSim3D_background_integrator function")
  
  message_shape_choice <- paste("Hello spaSim-3D user, hopefully you can see a plot of your current spe object. What type of shape do you want your cluster to be?\n
          1. Sphere\n
          2. Ellipsoid\n
          3. Cylinder\n
          4. Network\n\n",
                                "To choose, please enter 1, 2, 3 or 4.\n", sep = "")
  
  message_sphere_cluster <- paste("We will need a few parameters to generate a sphere cluster\n",
                                  "    Radius\n",
                                  "    Coordinates of sphere centre: x, y and z\n",
                                  "If you want to change your inputs, you'll be able to at the end.\n", sep = "")
  
  message_ellipsoid_cluster <- paste("We will need a few parameters to generate a ellipsoid cluster\n",
                                     "    Radii: x, y and z\n",
                                     "    Coordinates of ellipsoid centre: x, y and z\n",
                                     "    Angle of rotation in the x-axis, y-axis and z-axis\n",
                                     "If you want to change your inputs, you'll be able to at the end.\n", sep = "")
  
  message_cylinder_cluster <- paste("We will need a few parameters to generate a cylinder cluster\n",
                                    "    Radius\n",
                                    "    Coordinates of the cylinder start point: x, y and z\n",
                                    "    Coordinates of the cylinder end point: x, y and z\n",
                                    "If you want to change your inputs, you'll be able to at the end.\n", sep = "")
  
  message_network_cluster <- paste("We will need a few parameters to generate a network cluster\n",
                                   "    Number of branches\n",
                                   "    Width of each branch: x, y and z\n",
                                   "    Radius spanned by the whole network: x, y and z\n",
                                   "    Coordinates of network centre: x, y and z\n",
                                   "If you want to change your inputs, you'll be able to at the end.\n", sep = "")
  
  
  message_cluster_choice <- paste("You can customise your cluster further if you'd like:\n
          1. Add a ring\n
          2. Add a double ring\n
          3. Continue\n\n",
                                  "To choose, please enter 1, 2 or 3.\n", sep = "")
  
  
  ## Start with checking if the user has inputted spe object
  if (class(simulated_spe) != "SpatialExperiment") {
    stop(message_no_simulated_spe)
  }
  
  ## Plot the user's data so they can see what they already have
  fig <- plot_cells3D(simulated_spe)
  methods::show(fig)
  
  ## Get user's choice for shape type (sphere, ellipsoid, cylinder or network)
  message(message_shape_choice)
  user_input_shape <- get_integer_input_from_options(1:4)
  
  ### Sphere
  if (user_input_shape == 1) {
    # Get required parameters for a sphere cluster from user
    message(message_sphere_cluster)
    parameter_values <- list("radius" = get_positive_numeric_input("radius"),
                             "centre x coordinate" = get_non_negative_numeric_input("centre x coordinate"),
                             "centre y coordinate" = get_non_negative_numeric_input("centre y coordinate"),
                             "centre z coordinate" = get_non_negative_numeric_input("centre z coordinate"))
    display_parameters(parameter_values)
    
    # Generate sphere cluster simulation using these parameters
    cluster_properties <- list(list(shape = "sphere",
                                    cluster_cell_types = "Cluster",
                                    cluster_cell_proportions = 1,
                                    radius = parameter_values[["radius"]],
                                    centre_loc = c(parameter_values[["centre x coordinate"]],
                                                   parameter_values[["centre y coordinate"]],
                                                   parameter_values[["centre z coordinate"]])))
    message("Generating simulation...")
    simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                             cluster_properties,
                                             plot_image = TRUE,
                                             plot_cell_types = NULL,
                                             plot_colours = NULL)
    
    # Allow user the option to change their input parameters
    message("Would you like to change your input parameters?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      user_input_parameter_choice <- get_integer_input_from_options(seq(length(parameter_values)))
      
      if (user_input_parameter_choice == 1) parameter_values[["radius"]] <- get_positive_numeric_input("radius")
      if (user_input_parameter_choice == 2) parameter_values[["centre x coordinate"]] <- get_non_negative_numeric_input("centre x coordinate")
      if (user_input_parameter_choice == 3) parameter_values[["centre y coordinate"]] <- get_non_negative_numeric_input("centre y coordinate")
      if (user_input_parameter_choice == 4) parameter_values[["centre z coordinate"]] <- get_non_negative_numeric_input("centre z coordinate")
      
      display_parameters(parameter_values)
      
      # Generate sphere cluster simulation using updated parameters
      cluster_properties <- list(list(shape = "sphere",
                                      cluster_cell_types = "Cluster",
                                      cluster_cell_proportions = 1,
                                      radius = parameter_values[["radius"]],
                                      centre_loc = c(parameter_values[["centre x coordinate"]],
                                                     parameter_values[["centre y coordinate"]],
                                                     parameter_values[["centre z coordinate"]])))
      message("Generating simulation...")
      simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                               cluster_properties,
                                               plot_image = TRUE,
                                               plot_cell_types = NULL,
                                               plot_colours = NULL)
      
      message("Would you like to change your inputs?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
  }
  ### Ellipsoid
  else if (user_input_shape == 2) {
    # Get required parameters for an ellipsoid cluster from user
    message(message_ellipsoid_cluster)
    parameter_values <- list("x radius" = get_positive_numeric_input("x radius"),
                             "y radius" = get_positive_numeric_input("y radius"),
                             "z radius" = get_positive_numeric_input("z radius"),
                             "centre x coordinate" = get_non_negative_numeric_input("centre x coordinate"),
                             "centre y coordinate" = get_non_negative_numeric_input("centre y coordinate"),
                             "centre z coordinate" = get_non_negative_numeric_input("centre z coordinate"),
                             "x-axis rotation angle" = get_non_negative_numeric_input("x-axis rotation angle"),
                             "y-axis rotation angle" = get_non_negative_numeric_input("y-axis rotation angle"),
                             "z-axis rotation angle" = get_non_negative_numeric_input("z-axis rotation angle"))
    display_parameters(parameter_values)
    
    # Generate ellipsoid cluster simulation using these parameters
    cluster_properties <- list(list(shape = "ellipsoid",
                                    cluster_cell_types = "Cluster",
                                    cluster_cell_proportions = 1,
                                    radii = c(parameter_values[["x radius"]], 
                                              parameter_values[["y radius"]], 
                                              parameter_values[["z radius"]]),
                                    centre_loc = c(parameter_values[["centre x coordinate"]],
                                                   parameter_values[["centre y coordinate"]],
                                                   parameter_values[["centre z coordinate"]]),
                                    axes_rotation = c(parameter_values[["x-axis rotation angle"]], 
                                                      parameter_values[["y-axis rotation angle"]], 
                                                      parameter_values[["z-axis rotation angle"]])))
    message("Generating simulation...")
    simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                             cluster_properties,
                                             plot_image = TRUE,
                                             plot_cell_types = NULL,
                                             plot_colours = NULL)
    
    # Allow user the option to change their input parameters
    message("Would you like to change your input parameters?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      user_input_parameter_choice <- get_integer_input_from_options(seq(length(parameter_values)))
      
      if (user_input_parameter_choice == 1) parameter_values[["x radius"]] <- get_positive_numeric_input("x radius")
      if (user_input_parameter_choice == 2) parameter_values[["y radius"]] <- get_positive_numeric_input("y radius")
      if (user_input_parameter_choice == 3) parameter_values[["z radius"]] <- get_positive_numeric_input("z radius")
      if (user_input_parameter_choice == 4) parameter_values[["centre x coordinate"]] <- get_non_negative_numeric_input("centre x coordinate")
      if (user_input_parameter_choice == 5) parameter_values[["centre y coordinate"]] <- get_non_negative_numeric_input("centre y coordinate")
      if (user_input_parameter_choice == 6) parameter_values[["centre z coordinate"]] <- get_non_negative_numeric_input("centre z coordinate")
      if (user_input_parameter_choice == 7) parameter_values[["x-axis rotation angle"]] <- get_non_negative_numeric_input("x-axis rotation angle")
      if (user_input_parameter_choice == 8) parameter_values[["y-axis rotation angle"]] <- get_non_negative_numeric_input("y-axis rotation angle")
      if (user_input_parameter_choice == 9) parameter_values[["z-axis rotation angle"]] <- get_non_negative_numeric_input("z-axis rotation angle")
      
      display_parameters(parameter_values)
      
      # Generate ellipsoid cluster simulation using updated parameters
      cluster_properties <- list(list(shape = "ellipsoid",
                                      cluster_cell_types = "Cluster",
                                      cluster_cell_proportions = 1,
                                      radii = c(parameter_values[["x radius"]], 
                                                parameter_values[["y radius"]], 
                                                parameter_values[["z radius"]]),
                                      centre_loc = c(parameter_values[["centre x coordinate"]],
                                                     parameter_values[["centre y coordinate"]],
                                                     parameter_values[["centre z coordinate"]]),
                                      axes_rotation = c(parameter_values[["x-axis rotation angle"]], 
                                                        parameter_values[["y-axis rotation angle"]], 
                                                        parameter_values[["z-axis rotation angle"]])))
      
      message("Generating simulation...")
      simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                               cluster_properties,
                                               plot_image = TRUE,
                                               plot_cell_types = NULL,
                                               plot_colours = NULL)
      
      message("Would you like to change your inputs?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
  }
  ### Cylinder
  else if (user_input_shape == 3) {
    # Get required parameters for a cylinder cluster from user
    message(message_cylinder_cluster)
    parameter_values <- list("radius" = get_positive_numeric_input("radius"),
                             "start x coordinate" = get_non_negative_numeric_input("start x coordinate"),
                             "start y coordinate" = get_non_negative_numeric_input("start y coordinate"),
                             "start z coordinate" = get_non_negative_numeric_input("start z coordinate"),
                             "end x coordinate" = get_non_negative_numeric_input("end x coordinate"),
                             "end y coordinate" = get_non_negative_numeric_input("end y coordinate"),
                             "end z coordinate" = get_non_negative_numeric_input("end z coordinate"))
    display_parameters(parameter_values)
    
    # Generate cylinder cluster simulation using these parameters
    cluster_properties <- list(list(shape = "cylinder",
                                    cluster_cell_types = "Cluster",
                                    cluster_cell_proportions = 1,
                                    radius = parameter_values[["radius"]],
                                    start_loc = c(parameter_values[["start x coordinate"]],
                                                  parameter_values[["start y coordinate"]],
                                                  parameter_values[["start z coordinate"]]),
                                    end_loc = c(parameter_values[["end x coordinate"]],
                                                parameter_values[["end y coordinate"]],
                                                parameter_values[["end z coordinate"]])))
    message("Generating simulation...")
    simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                             cluster_properties,
                                             plot_image = TRUE,
                                             plot_cell_types = NULL,
                                             plot_colours = NULL)
    
    # Allow user the option to change their input parameters
    message("Would you like to change your input parameters?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      user_input_parameter_choice <- get_integer_input_from_options(seq(length(parameter_values)))
      
      if (user_input_parameter_choice == 1) parameter_values[["radius"]] <- get_positive_numeric_input("radius")
      if (user_input_parameter_choice == 2) parameter_values[["start x coordinate"]] <- get_non_negative_numeric_input("start x coordinate")
      if (user_input_parameter_choice == 3) parameter_values[["start y coordinate"]] <- get_non_negative_numeric_input("start y coordinate")
      if (user_input_parameter_choice == 4) parameter_values[["start z coordinate"]] <- get_non_negative_numeric_input("start z coordinate")
      if (user_input_parameter_choice == 5) parameter_values[["end x coordinate"]] <- get_non_negative_numeric_input("end x coordinate")
      if (user_input_parameter_choice == 6) parameter_values[["end y coordinate"]] <- get_non_negative_numeric_input("end y coordinate")
      if (user_input_parameter_choice == 7) parameter_values[["end z coordinate"]] <- get_non_negative_numeric_input("end z coordinate")
      
      display_parameters(parameter_values)
      
      # Generate cylinder cluster simulation using updated parameters
      cluster_properties <- list(list(shape = "cylinder",
                                      cluster_cell_types = "Cluster",
                                      cluster_cell_proportions = 1,
                                      radius = parameter_values[["radius"]],
                                      start_loc = c(parameter_values[["start x coordinate"]],
                                                    parameter_values[["start y coordinate"]],
                                                    parameter_values[["start z coordinate"]]),
                                      end_loc = c(parameter_values[["end x coordinate"]],
                                                  parameter_values[["end y coordinate"]],
                                                  parameter_values[["end z coordinate"]])))
      message("Generating simulation...")
      simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                               cluster_properties,
                                               plot_image = TRUE,
                                               plot_cell_types = NULL,
                                               plot_colours = NULL)
      
      message("Would you like to change your inputs?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
  }
  ### Network
  else if (user_input_shape == 4) {
    # Get required parameters for a network cluster from user
    message(message_network_cluster)
    parameter_values <- list("number of branches" = get_integer_greater_than_or_equal_input("number of branches", 2),
                             "width of branch" = get_positive_numeric_input("width of branch"),
                             "radius spanned by network" = get_positive_numeric_input("radius spanned by network"),
                             "centre x coordinate" = get_non_negative_numeric_input("centre x coordinate"),
                             "centre y coordinate" = get_non_negative_numeric_input("centre y coordinate"),
                             "centre z coordinate" = get_non_negative_numeric_input("centre z coordinate"))
    display_parameters(parameter_values)
    
    # Generate network cluster simulation using these parameters
    cluster_properties <- list(list(shape = "network",
                                    cluster_cell_types = "Cluster",
                                    cluster_cell_proportions = 1,
                                    n_edges = parameter_values[["number of branches"]],
                                    width = parameter_values[["width of branch"]],
                                    radius = parameter_values[["radius spanned by network"]],
                                    centre_loc = c(parameter_values[["centre x coordinate"]],
                                                   parameter_values[["centre y coordinate"]],
                                                   parameter_values[["centre z coordinate"]])))
    message("Generating simulation...")
    simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                             cluster_properties,
                                             plot_image = TRUE,
                                             plot_cell_types = NULL,
                                             plot_colours = NULL)
    
    # Allow user the option to change their input parameters
    message("Would you like to change your input parameters?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      user_input_parameter_choice <- get_integer_input_from_options(seq(length(parameter_values)))
      
      if (user_input_parameter_choice == 1) parameter_values[["number of branches"]] <- get_integer_greater_than_or_equal_input("number of branches", 2)
      if (user_input_parameter_choice == 2) parameter_values[["width of branch"]] <- get_positive_numeric_input("width of branch")
      if (user_input_parameter_choice == 3) parameter_values[["radius spanned by network"]] <- get_positive_numeric_input("radius spanned by network")
      if (user_input_parameter_choice == 4) parameter_values[["centre x coordinate"]] <- get_non_negative_numeric_input("centre x coordinate")
      if (user_input_parameter_choice == 5) parameter_values[["centre y coordinate"]] <- get_non_negative_numeric_input("centre y coordinate")
      if (user_input_parameter_choice == 6) parameter_values[["centre z coordinate"]] <- get_non_negative_numeric_input("centre z coordinate")
      
      display_parameters(parameter_values)
      
      # Generate sphere cluster simulation using updated parameters
      cluster_properties <- list(list(shape = "network",
                                      cluster_cell_types = "Cluster",
                                      cluster_cell_proportions = 1,
                                      n_edges = parameter_values[["number of branches"]],
                                      width = parameter_values[["width of branch"]],
                                      radius = parameter_values[["radius spanned by network"]],
                                      centre_loc = c(parameter_values[["centre x coordinate"]],
                                                     parameter_values[["centre y coordinate"]],
                                                     parameter_values[["centre z coordinate"]])))
      
      message("Generating simulation...")
      simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                               cluster_properties,
                                               plot_image = TRUE,
                                               plot_cell_types = NULL,
                                               plot_colours = NULL)
      
      message("Would you like to change your inputs?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
  }
  
  # Allow user to change the cell composition of the cluster
  message("Let's change the cell composition of this cluster")
  simulated_spe_new_and_properties <- get_cell_types_and_proportions_for_clusters(simulated_spe_new,
                                                                                  simulate_clusters3D,
                                                                                  cluster_properties,
                                                                                  "cluster_cell_types",
                                                                                  "cluster_cell_proportions",
                                                                                  "Cluster")
  simulated_spe_new <- simulated_spe_new_and_properties[["data"]]
  cluster_properties <- simulated_spe_new_and_properties[["properties"]]
  
  
  
  ## Get user's choice for cluster type (ringed, double ringed or continue)
  message(message_cluster_choice)
  user_input_cluster <- get_integer_input_from_options(1:3)
  
  ### Ring
  if (user_input_cluster == 1) {
    # Get width of ring from user
    message("For a single ring, we needs its width.\n")
    
    # Generate cluster with ring simulation using this width
    cluster_properties[[1]][["ring_width"]] <- get_positive_numeric_input("ring width")
    cluster_properties[[1]][["ring_cell_types"]] <- c("Ring")
    cluster_properties[[1]][["ring_cell_proportions"]] <- 1
    
    message("Generating simulation...")
    simulated_spe_new <- simulate_rings3D(simulated_spe,
                                          cluster_properties,
                                          plot_image = TRUE,
                                          plot_cell_types = NULL,
                                          plot_colours = NULL)
    
    # Allow user the option to change the ring width
    message("Would you like to change the ring width?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      cluster_properties[[1]][["ring_width"]] <- get_positive_numeric_input("ring width")
      
      message("Generating simulation...")
      simulated_spe_new <- simulate_rings3D(simulated_spe,
                                            cluster_properties,
                                            plot_image = TRUE,
                                            plot_cell_types = NULL,
                                            plot_colours = NULL)
      
      message("Would you like to the ring width?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
    
    # Allow user to change the cell composition of the ring
    message("Let's change the cell composition of the ring")
    
    simulated_spe_new_and_properties <- get_cell_types_and_proportions_for_clusters(simulated_spe_new,
                                                                                    simulate_rings3D,
                                                                                    cluster_properties,
                                                                                    "ring_cell_types",
                                                                                    "ring_cell_proportions",
                                                                                    "Ring")
    
    simulated_spe_new <- simulated_spe_new_and_properties[["data"]]
  }
  ### Double ring
  else if (user_input_cluster == 2) {
    # Get width of inner and outer ring from user
    message("For a double ring, we needs the width of the inner and outer ring.\n")
    
    # Generate cluster with double ring simulation using both widths
    parameter_values <- list("inner ring width" = get_positive_numeric_input("inner ring width"),
                             "outer ring width" = get_positive_numeric_input("outer ring width"))
    display_parameters(parameter_values)
    
    cluster_properties[[1]][["inner_ring_width"]] <- parameter_values[["inner ring width"]]
    cluster_properties[[1]][["outer_ring_width"]] <- parameter_values[["outer ring width"]]
    cluster_properties[[1]][["inner_ring_cell_types"]] <- c("Inner ring")
    cluster_properties[[1]][["inner_ring_cell_proportions"]] <- 1
    cluster_properties[[1]][["outer_ring_cell_types"]] <- c("Outer ring")
    cluster_properties[[1]][["outer_ring_cell_proportions"]] <- 1
    
    message("Generating simulation...")
    simulated_spe_new <- simulate_double_rings3D(simulated_spe,
                                                 cluster_properties,
                                                 plot_image = TRUE,
                                                 plot_cell_types = NULL,
                                                 plot_colours = NULL)
    
    # Allow user the option to change the widths of the inner or outer ring
    message("Would you like to change the widths of the inner or outer ring?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      if (user_input_parameter_choice == 1) parameter_values[["inner ring width"]] <- get_positive_numeric_input("inner ring width")
      if (user_input_parameter_choice == 2) parameter_values[["outer ring width"]] <- get_positive_numeric_input("outer ring width")
      
      cluster_properties[[1]][["inner_ring_width"]] <- parameter_values[["inner ring width"]]
      cluster_properties[[1]][["outer_ring_width"]] <- parameter_values[["outer ring width"]]
      
      display_parameters(parameter_values)
      
      message("Generating simulation...")
      simulated_spe_new <- simulate_double_rings3D(simulated_spe,
                                                   cluster_properties,
                                                   plot_image = TRUE,
                                                   plot_cell_types = NULL,
                                                   plot_colours = NULL)
      
      message("Would you like to change the widths of the inner or outer ring?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
    
    # Allow user to change the cell composition of the inner ring
    message("Let's change the cell composition of the inner ring")
    simulated_spe_new_and_properties <- get_cell_types_and_proportions_for_clusters(simulated_spe_new,
                                                                                    simulate_double_rings3D,
                                                                                    cluster_properties,
                                                                                    "inner_ring_cell_types",
                                                                                    "inner_ring_cell_proportions",
                                                                                    "Inner ring")
    
    simulated_spe_new <- simulated_spe_new_and_properties[["data"]]
    cluster_properties <- simulated_spe_new_and_properties[["properties"]]
    
    # Allow user to change the cell composition of the outer ring
    message("Let's change the cell composition of the outer ring")
    simulated_spe_new_and_properties <- get_cell_types_and_proportions_for_clusters(simulated_spe_new,
                                                                                    simulate_double_rings3D,
                                                                                    cluster_properties,
                                                                                    "outer_ring_cell_types",
                                                                                    "outer_ring_cell_proportions",
                                                                                    "Outer ring")
    
    simulated_spe_new <- simulated_spe_new_and_properties[["data"]]
  }
  ### Continue
  else if (user_input_cluster == 3) {
    
  }
  
  message("All done!")
  return(simulated_spe_new) 
}

spe_metadata_background_template <- function(background_type, original_spe_metadata = NULL) {
  
  if (background_type == "random") {
    background_metadata <- list(background = list(background_type = "random",
                                                  n_cells = 20000,
                                                  length = 600,
                                                  width = 600,
                                                  height = 300,
                                                  minimum_distance_between_cells = 10,
                                                  cell_types = c("Tumour", "Others"),
                                                  cell_proportions = c(0.05, 0.95)))
  }
  else if (background_type == "ordered") {
    background_metadata <- list(background = list(background_type = "ordered",
                                                  n_cells = 20000,
                                                  length = 600,
                                                  width = 300,
                                                  height = 300,
                                                  jitter_proportion = 0.25,
                                                  cell_types = c("Immune", "Others"),
                                                  cell_proportions = c(0.05, 0.95)))
  }
  else {
    stop("background_type parameter must be 'random' or 'ordered'.")
  }
  
  
  # If original_spe_metadata input is not null, replace its background metadata with new background metadata
  if (!is.null(original_spe_metadata) && !is.null(original_spe_metadata[["background"]])) {
    original_spe_metadata[["background"]] <- background_metadata    
    return(original_spe_metadata)
  }
  else if (!is.null(original_spe_metadata) && is.null(original_spe_metadata[["background"]])) {
    original_spe_metadata <- c(background_metadata, original_spe_metadata)
    return(original_spe_metadata)
  }
  
  # Else, just return the background_metadata
  return(background_metadata)
}
spe_metadata_cluster_template <- function(cluster_type, shape, original_spe_metadata = NULL) {
  
  ### Get template for different shapes
  if (shape == "sphere") {
    cluster_metadata <- list(shape = "sphere",
                             cluster_cell_types = c("Tumour", "Immune", "Others"),
                             cluster_cell_proportions = c(0.8, 0.15, 0.05),
                             radius = 100,
                             centre_loc = c(200, 150, 200))
  }
  else if (shape == "ellipsoid") {
    cluster_metadata <- list(shape = "ellipsoid",
                             cluster_cell_types = c("Tumour", "Immune", "Others"),
                             cluster_cell_proportions = c(0.8, 0.15, 0.05),
                             radii = c(75, 100, 125),
                             centre_loc = c(450, 300, 100),
                             axes_rotation = c(0, 45, 0))
  }
  else if (shape == "cylinder") {
    cluster_metadata <- list(shape = "cylinder",
                             cluster_cell_types = c("Endothelial", "Others"),
                             cluster_cell_proportions = c(0.95, 0.05),
                             radius = 40,
                             start_loc = c(400, 0, 0),
                             end_loc   = c(600, 400, 200)) 
  }
  else if (shape == "network") {
    cluster_metadata <- list(shape = "network",
                             cluster_cell_types = c("Immune", "Others"),
                             cluster_cell_proportions = c(0.95, 0.05),
                             n_edges = 20,
                             width = 30,
                             centre_loc = c(200, 400, 150),
                             radius = 200)
  }
  else {
    stop("shape parameter must be 'sphere', 'ellipsoid', 'cylinder' or 'network'")
  }
  
  ### Add extra metadata for different cluster types
  if (cluster_type == "regular") {
    cluster_metadata <- append(list(cluster_type = "regular"), cluster_metadata)    
  }
  else if (cluster_type == "ring") {
    cluster_metadata <- append(list(cluster_type = "ring"), cluster_metadata)
    cluster_metadata$ring_cell_types <- c("Immune1", "Others")
    cluster_metadata$ring_cell_proportions <- c(0.85, 0.15)
    cluster_metadata$ring_width <- 12
  }
  else if (cluster_type == "double ring") {
    cluster_metadata <- append(list(cluster_type = "double ring"), cluster_metadata)
    cluster_metadata$inner_ring_cell_types <- c("Immune1", "Others")
    cluster_metadata$inner_ring_cell_proportions <- c(0.85, 0.15)
    cluster_metadata$inner_ring_width <- 10
    cluster_metadata$outer_ring_cell_types <- c("Immune2", "Others")
    cluster_metadata$outer_ring_cell_proportions <- c(0.85, 0.15)
    cluster_metadata$outer_ring_width <- 10
  }
  else {
    stop("cluster_type parameter must be 'regular', 'ring' or 'double ring'")
  }
  
  # If original_spe_metadata input is not null, add new cluster_metadata to it
  if (!is.null(original_spe_metadata) && !is.null(original_spe_metadata[["background"]])) {
    original_spe_metadata[[paste("cluster", length(original_spe_metadata), sep="_")]] <- cluster_metadata    
    return(original_spe_metadata)
  }
  else if (!is.null(original_spe_metadata) && is.null(original_spe_metadata[["background"]])) {
    original_spe_metadata[[paste("cluster", length(original_spe_metadata) + 1, sep="_")]] <- cluster_metadata
    return(original_spe_metadata)
  }
  
  # Else, just return the new cluster_metadata
  return(list("cluster_1" = cluster_metadata))
}

poisson_distribution3D <- function(n_cells, length, width, height)  {
  
  # Choose lambda
  lambda <- 5
  
  # Set number of rows, columns and layers
  nRows <- nCols <- nLays <- round((n_cells/lambda)^(1/3))
  
  # Get number of cubes in grid
  nCubes <- nRows * nCols * nLays
  
  # Get pois vector
  pois <- rpois(nCubes, lambda)
  
  # Get points for each prism region
  x <- c()
  y <- c()
  z <- c()
  
  for (row in seq(nRows)) {
    
    for (col in seq(nCols)) {
      
      for (lay in seq(nLays)) {
        current_cube_index <- nRows^2 * (row - 1) + nCols * (col - 1) + lay
        
        x <- append(x, runif(pois[current_cube_index], row - 1, row))
        y <- append(y, runif(pois[current_cube_index], col - 1, col))
        z <- append(z, runif(pois[current_cube_index], lay - 1, lay))
      }
    }
  }
  x <- x * length / nRows
  y <- y * width / nCols
  z <- z * height / nLays
  
  df <- data.frame("Cell.X.Position" = x, 
                   "Cell.Y.Position" = y, 
                   "Cell.Z.Position" = z)
  
  return(df)
}

## Prim's algorithm function
# Input is the adjacency matrix of the graph (i.e. output from -1 * apcluster::negDistMat(df of coords))
prims_algorithm <- function(graph) {
  
  # Number of vertices is number of points
  num_vertices <- nrow(graph)
  
  # Start with no vertices selected except first
  selected <- rep(FALSE, num_vertices)
  selected[1] <- TRUE
  
  # Create tree_edge matrix. Currently zero, each row represents the two vertices the edge joins
  tree_edges <- matrix(0, 
                       nrow = num_vertices - 1,
                       ncol = 2)
  
  # Iterate until we select enough edges (one less than the number of vertices for a MST)
  num_edges <- 0
  while (num_edges < num_vertices - 1) {
    # Set initial temp values for weight and vertex
    min_weight <- Inf
    min_vertex <- -1
    
    # Iterate through each currently selected vertex
    for (i in seq(num_vertices)) {
      
      # Found a currently selected vertex
      if (selected[i] == TRUE) {
        
        # Iterate through each unselected vertex and find the nearest one
        for (j in seq(num_vertices)) {
          if (!selected[j] && graph[i, j] < min_weight) {
            min_weight <- graph[i, j]
            min_vertex <- j
            curr_vertex <- i
          }
        }
      }
    }
    
    # Current edge connects the min_vertex and curr_vertex
    tree_edges[num_edges + 1, ] <- c(min_vertex, curr_vertex)
    selected[min_vertex] <- TRUE
    num_edges <- num_edges + 1
  }
  return(tree_edges)
}

get_tree_depth <- function(tree_edges) {
  
  tree_edges <- data.frame(tree_edges)
  colnames(tree_edges) <- c("vertex1", "vertex2")
  
  # Set the initial depth of each tree_edge to be NA.
  tree_edges$depth <- NA
  
  # Get vertices on the 'outskirts' of MST (leaf_vertices which have a depth of 1)
  tree_vertices <- c(tree_edges[ , 1], tree_edges[ , 2])
  leaf_vertices <- as.numeric(names(table(tree_vertices))[table(tree_vertices) == 1])
  
  # Start with leaf_vertices
  curr_vertices <- leaf_vertices
  curr_depth <- 1
  
  while (NA %in% tree_edges$depth) {
    
    # New vertices will be those adjacent to the current vertices
    new_vertices <- c()
    
    # Check each current vertex
    for (vertex in curr_vertices) {
      # Start with vertex1
      curr_edges <- which(tree_edges$vertex1 == vertex)
      tree_edges[curr_edges, "depth"][is.na(tree_edges[curr_edges, "depth"])] <- curr_depth
      new_vertices <- c(new_vertices, tree_edges[curr_edges, "vertex2"])
      
      # Then vertex2
      curr_edges <- which(tree_edges$vertex2 == vertex)
      tree_edges[curr_edges, "depth"][is.na(tree_edges[curr_edges, "depth"])] <- curr_depth
      new_vertices <- c(new_vertices, tree_edges[curr_edges, "vertex1"])
      
      # Only keep unique vertices
      new_vertices <- unique(new_vertices)
    }
    
    curr_depth <- curr_depth + 1
    curr_vertices <- new_vertices
  }
  
  return(tree_edges)
}

check_input_parameters <- function(input_parameters) {
  
  input_parameter_names <- names(input_parameters)
  
  check_value <- 0
  
  for (input_parameter_name in input_parameter_names) {
    input_parameter <- input_parameters[[input_parameter_name]]
    
    # spe
    if (input_parameter_name == "spe" && class(input_parameter) != "SpatialExperiment") {
      check_value <- 1
      break
    }
    # Positive integer
    if (input_parameter_name %in% c("n_cells", "n_edges") && !(is.integer(input_parameter) && length(input_parameter) == 1 || (is.numeric(input_parameter) && length(input_parameter) == 1 && input_parameter > 0 && input_parameter%%1 == 0))) {
      check_value <- 2
      break
    }  
    # Positive numeric
    if (input_parameter_name %in% c("length", "width", "height", "radius", "x_radius", "y_radius", "z_radius", "ring_width", "inner_ring_width", "outer_ring_width") && !(is.numeric(input_parameter) && length(input_parameter) == 1 && input_parameter > 0)) {
      check_value <- 3
      break
    }
    # Non-negative numeric
    if (input_parameter_name %in% c("minimum_distance_between_cells") && !(is.numeric(input_parameter) && length(input_parameter) == 1 && input_parameter >= 0)) {
      check_value <- 4
      break
    }
    # Numeric between 0 and 1
    if (input_parameter_name %in% c("jitter_proportion") && !(is.numeric(input_parameter) && length(input_parameter) == 1 && input_parameter >= 0 && input_parameter <= 1)) {
      check_value <- 5
      break
    }
    # Character
    if (input_parameter_name %in% c("background_cell_type") && !(is.character(input_parameter)) && length(input_parameter) == 1) {
      check_value <- 6
      break
    }
    # Logical
    if (input_parameter_name %in% c("plot_image") && !(is.logical(input_parameter)) && length(input_parameter) == 1) {
      check_value <- 7
      break
    }
    # Character vector
    if (input_parameter_name %in% c("cell_types", "cluster_cell_types", "ring_cell_types", "inner_ring_cell_types", "outer_ring_cell_types") && 
        !(is.character(input_parameter))) {
      check_value <- 8
      break
    }
    # Numeric vector
    if (input_parameter_name %in% c("cell_proportions", "cluster_cell_proportions", "ring_cell_proportions", "inner_ring_cell_proportions", "outer_ring_cell_proportions") && 
        !(is.numeric(input_parameter))) {
      check_value <- 9
      break
    }
    # Numeric vector contains values between 0 and 1
    if (input_parameter_name %in% c("cell_proportions", "cluster_cell_proportions", "ring_cell_proportions", "inner_ring_cell_proportions", "outer_ring_cell_proportions") && 
        sum(input_parameter < 0 | input_parameter > 1) != 0) {
      check_value <- 10
      break
    }
    # Numeric vector contains values that sum to 1
    if (input_parameter_name %in% c("cell_proportions", "cluster_cell_proportions", "ring_cell_proportions", "inner_ring_cell_proportions", "outer_ring_cell_proportions") && 
        !is_equal_with_tolerance(sum(input_parameter), 1)) {
      check_value <- 11
      break
    }
    # Numeric vector of length 3
    if (input_parameter_name %in% c("centre_loc", "start_loc", "end_loc") && !(is.numeric(input_parameter) && length(input_parameter) == 3)) {
      check_value <- 12
      break
    }
    # Numeric
    if (input_parameter_name %in% c("y_z_rotation", "x_z_rotation", "x_y_rotation") && !(is.numeric(input_parameter) && length(input_parameter) == 1)) {
      check_value <- 13
      break
    }
  }
  
  # Two vectors match in length
  if (check_value == 0) {
    pairs <- data.frame(name1 = c("cell_types", "cluster_cell_types", "ring_cell_types", "inner_ring_cell_types", "outer_ring_cell_types"),
                        name2 = c("cell_proportions", "cluster_cell_proportions", "ring_cell_proportions", "inner_ring_cell_proportions", "outer_ring_cell_proportions"))  
    
    for (i in seq(nrow(pairs))) {
      name1 <- pairs[["name1"]][i]
      name2 <- pairs[["name2"]][i]
      if (name1 %in% input_parameter_names && name2 %in% input_parameter_names && length(input_parameters[[name1]]) != length(input_parameters[[name2]])) {
        check_value <- 14
        input_parameter_name <- c(name1, name2)
        break
      }
    }
  }
  
  # If check_value equals 0, all inputs are valid.
  if (check_value == 0) {
    return(TRUE)
  }
  # At least one input is not valid, return the first invalid input. 
  else {
    return(list(input_parameter_name = input_parameter_name, check_value = check_value)) 
  }
}

input_parameter_error_message <- function(input_parameter_check_value) {
  
  input_parameter_name <- input_parameter_check_value[[1]]
  check_value <- input_parameter_check_value[[2]]
  
  error_message <- switch(check_value,
                          "1" = paste(input_parameter_name, "is not a SpatialExperiment object."),
                          "2" = paste(input_parameter_name, "is not a positive integer."),
                          "3" = paste(input_parameter_name, "is not a positive numeric."),
                          "4" = paste(input_parameter_name, "is not a non-negative numeric."),
                          "5" = paste(input_parameter_name, "is not a numeric between 0 and 1."),
                          "6" = paste(input_parameter_name, "is not a character."),
                          "7" = paste(input_parameter_name, "is not a logical (TRUE or FALSE)."),
                          "8" = paste(input_parameter_name, "is not a character vector."),
                          "9" = paste(input_parameter_name, "is not a numeric vector."),
                          "10" = paste(input_parameter_name, "cannot be negative or greater than 1."),
                          "11" = paste(input_parameter_name, "does not sum to 1."),
                          "12" = paste(input_parameter_name, "is not a numeric vector of length 3."),
                          "13" = paste(input_parameter_name, "is not a numeric."),
                          "14" = paste(input_parameter_name[1], "and", input_parameter_name[2], "do not match in length."))
  
  return(error_message)
}

is_equal_with_tolerance <- function(x, y, tolerance = 1e-6) {
  abs(x - y) <= tolerance
}

get_integer_greater_than_or_equal_input <- function(parameter, lower) {
  
  prompt <- paste("Enter an integer in value greater than or equal to ", lower, " for the ", parameter, ": ", sep = "")
  
  valid_input <- FALSE
  while (!valid_input) {
    user_input <- readline(prompt = prompt)
    # Try converting to numeric
    integer_value <- tryCatch({as.numeric(user_input)}, error = function(e) NA)
    
    # Non-numeric input
    if (is.na(integer_value)) {
      message("Invalid input. Please enter a numeric integer value.")
    }
    # Numeric but not integer
    else if (integer_value%%1 != 0) {
      message("Non-integer input. Please enter an integer value.")
    }
    # Integer but below the lower bound
    else if (integer_value < lower) {
      message("Out of bounds input. Please a number greater than or equal to ", lower)
    }
    else {
      valid_input <- TRUE
      message("Valid input received!")
    }
  }
  
  return(integer_value)
}

get_integer_input_from_options <- function(integer_options) {
  
  first_integers <- integer_options[1:(length(integer_options) - 1)]
  last_integer <- integer_options[length(integer_options)]
  integers_string <- paste(paste(first_integers, collapse = ", "), "or", last_integer)
  
  prompt <- paste("Enter either ", integers_string, ": ", sep = "")
  invalid_input_message <- paste("Invalid input. Please enter only", integers_string)
  
  valid_input <- FALSE
  while (!valid_input) {
    user_input <- readline(prompt = prompt)
    # Try converting to integer
    int_value <- tryCatch({as.integer(user_input)}, error = function(e) NA)
    
    # Check if conversion was successful and value is in integer_options
    if (!is.na(int_value) && int_value %in% integer_options) {
      valid_input <- TRUE
      message("Valid input received!")
    } 
    else {
      message(invalid_input_message)
    }
  }
  
  return(int_value)
}

get_non_negative_numeric_input <- function(parameter) {
  
  prompt <- paste("Enter a non-negative numeric value for the ", parameter, ": ", sep = "")
  
  valid_input <- FALSE
  while (!valid_input) {
    user_input <- readline(prompt = prompt)
    # Try converting to numeric
    non_negative_value <- tryCatch({as.numeric(user_input)}, error = function(e) NA)
    
    # Non-numeric input
    if (is.na(non_negative_value)) {
      message("Invalid input. Please enter a numeric value.")
    }
    # Negative input
    else if (non_negative_value < 0) {
      message("Negative input. Please enter a non-negative number") 
    }
    # Should be correct input
    else {
      valid_input <- TRUE
      message("Valid input received!") 
    }
  }
  
  return(non_negative_value)
}

get_numeric_between_input <- function(parameter, lower, upper) {
  
  prompt <- paste("Enter a numeric value between ", lower, " and ", upper, " for the ", parameter, ": ", sep = "")
  
  valid_input <- FALSE
  while (!valid_input) {
    user_input <- readline(prompt = prompt)
    # Try converting to numeric
    numeric_value <- tryCatch({as.numeric(user_input)}, error = function(e) NA)
    
    # Non-numeric input
    if (is.na(numeric_value)) {
      message("Invalid input. Please enter a numeric value.")
    }
    # Out of bounds input
    else if (numeric_value < lower || numeric_value > upper) {
      message("Out of bounds input. Please a number between ", lower, " and ", upper, ".", sep = "")
    }
    # Should be correct
    else {
      valid_input <- TRUE
      message("Valid input received!")
    }
  }
  
  return(numeric_value)
}

get_positive_numeric_input <- function(parameter) {
  
  prompt <- paste("Enter a positive numeric value for the ", parameter, ": ", sep = "")
  
  valid_input <- FALSE
  while (!valid_input) {
    user_input <- readline(prompt = prompt)
    # Try converting to numeric
    positive_numeric_value <- tryCatch({as.numeric(user_input)}, error = function(e) NA)
    
    # Non-numeric input
    if (is.na(positive_numeric_value)) {
      message("Invalid input. Please enter a numeric value.")
    }
    # Non-positive input
    else if (positive_numeric_value <= 0) {
      message("Non-positive input. Please enter a positive number") 
    }
    # Should be correct input
    else {
      valid_input <- TRUE
      message("Valid input received!") 
    }
  }
  
  return(positive_numeric_value)
}

get_y_or_n_input <- function() {
  
  valid_input <- FALSE
  while (!valid_input) {
    user_input <- readline(prompt = "Enter either y or n: ")
    
    if (user_input %in% c("y", "n")) {
      valid_input <- TRUE
      message("Valid input received!")
    }
    else {
      message("Invalid input. Please enter either y or n.")
    }
  }
  
  return(user_input)
}

display_parameters <- function(parameter_values) {
  
  message("Your current inputs are:\n")
  
  display_message <- ""
  
  for (i in seq(length(parameter_values))) {
    display_message <- paste(display_message, "    ", i, ". ", names(parameter_values)[i], ": ", parameter_values[[i]], '\n', sep = "")
  }
  message(display_message)
}

get_cell_types_and_proportions_for_mixing <- function(simulated_spe) {
  
  message_get_cell_types <- "Keep entering the name of cell types you would like (e.g. Tumour, Immune, etc.).\n    enter 'stop' to move on."
  
  ## Get cell types from user
  cell_types <- c()
  user_input <- ""
  message(message_get_cell_types)
  while (user_input != "stop") {
    
    user_input <- readline(prompt = "Enter a cell type, or enter 'stop': ")
    
    ## Ignore if user enters a blank string
    if (user_input == "") {
      
    }
    ## Add inputted cell type to cell_types vector
    else if (user_input != "stop") {
      cell_types <- c(cell_types, user_input)
      message(paste("Cell type added:", user_input))
    }
    ## User wants to stop but hasn't entered any cell types
    else if (user_input == "stop" && length(cell_types) == 0) {
      message("You have not entered any cell types. Try again\n")
      user_input <- ""
    }
    ## User wants to stop
    else {
      message(paste("Your cell types chosen are:", paste(cell_types, collapse = ", ")))
      
      ## Allow user to re-choose cell types
      message("Would like to re-choose these cell types?\n")
      user_input_y_or_n <- get_y_or_n_input()
      if (user_input_y_or_n == "y") {
        cell_types <- c()
        message(message_get_cell_types)
        user_input <- ""
      }
    }
  }
  
  ## Get cell proportions from user
  cell_proportions <- c()
  max_proportion <- 1
  i <- 1
  message("For each cell type, choose their proportion in the simulation. They must add to 1.\n")
  while (i <= length(cell_types)) {
    
    ## For the last cell type, we can figure out what the cell proportion must be
    if (i == length(cell_types)) {
      cell_proportions <- c(cell_proportions, max_proportion)
      message("Cell proportion for ", cell_types[i], " must be ", round(max_proportion, 5))
    }
    ## Add inputted cell proportion to cell_proportions vector
    else {
      cell_proportion <- get_numeric_between_input(paste("cell proportion of", cell_types[i], "cells"), 0, max_proportion)
      cell_proportions <- c(cell_proportions, cell_proportion)
      max_proportion <- 1 - sum(cell_proportions)
      message("Cell proportion for ", cell_types[i], " is ", cell_proportion)
    }
    i <- i + 1
    
    if (i > length(cell_types)) {
      ## Generate simulation
      message("Generating simulation...")
      simulated_spe <- simulate_mixing3D(simulated_spe,
                                         cell_types,
                                         cell_proportions,
                                         plot_image = F)
      
      fig <- plot_cells3D(simulated_spe)
      print(fig)
      
      if (length(cell_types) == 1) break # If there is only one cell type, proportion is always 1
      
      ## Allow user to re-choose cell proportions  
      message("Would like to re-choose these cell proportions?\n")
      user_input_y_or_n <- get_y_or_n_input()
      if (user_input_y_or_n == "y") {
        cell_proportions <- c()
        max_proportion <- 1
        i <- 1
        message("For each cell type, choose their proportion in the simulation. They must add to 1.\n")
      }
    }
  }
  
  return(simulated_spe)
}

get_cell_types_and_proportions_for_clusters <- function(simulated_spe, simulate_function, properties, cell_type_option, cell_proportion_option, temp_cell_type) {
  
  message_get_cell_types <- "Keep entering the name of cell types you would like (e.g. Tumour, Immune, etc.).\n    enter 'stop' to move on."
  
  ## Display the cell types currently found in simulated_spe to the user
  current_cell_types <- setdiff(unique(simulated_spe[["Cell.Type"]]), temp_cell_type)
  message("Your data currently has the following cell types:\n", paste(current_cell_types, collapse = ", "), "\n")
  
  ## Get cell types from user
  cell_types <- c()
  user_input <- ""
  message(message_get_cell_types)
  while (user_input != "stop") {
    
    user_input <- readline(prompt = "Enter a cell type, or enter 'stop': ")
    
    ## Ignore if user enters a blank string
    if (user_input == "") {
      
    }
    ## Add inputted cell type to cell_types vector
    else if (user_input != "stop") {
      cell_types <- c(cell_types, user_input)
      message(paste("Cell type added:", user_input))
    }
    ## User wants to stop but hasn't entered any cell types
    else if (user_input == "stop" && length(cell_types) == 0) {
      message("You have not entered any cell types. Try again\n")
      user_input <- ""
    }
    ## User wants to stop
    else {
      message(paste("Your cell types chosen are:", paste(cell_types, collapse = ", ")))
      
      ## Allow user to re-choose cell types
      message("Would like to re-choose these cell types?\n")
      user_input_y_or_n <- get_y_or_n_input()
      if (user_input_y_or_n == "y") {
        message("Your data currently has the following cell types:\n", paste(current_cell_types, collapse = ", "), "\n")
        cell_types <- c()
        message(message_get_cell_types)
        user_input <- ""
      }
    }
  }
  properties[[1]][[cell_type_option]] <- cell_types
  
  ## Get cell proportions from user
  cell_proportions <- c()
  max_proportion <- 1
  i <- 1
  message("For each cell type, choose their proportion in the simulation. They must add to 1.\n")
  while (i <= length(cell_types)) {
    
    ## For the last cell type, we can figure out what the cell proportion must be
    if (i == length(cell_types)) {
      cell_proportions <- c(cell_proportions, max_proportion)
      message("Cell proportion for ", cell_types[i], " must be ", round(max_proportion, 5))
    }
    ## Add inputted cell proportion to cell_proportions vector
    else {
      cell_proportion <- get_numeric_between_input(paste("cell proportion of", cell_types[i], "cells"), 0, max_proportion)
      cell_proportions <- c(cell_proportions, cell_proportion)
      max_proportion <- 1 - sum(cell_proportions)
      message("Cell proportion for ", cell_types[i], " is ", cell_proportion)
    }
    i <- i + 1
    
    if (i > length(cell_types)) {
      properties[[1]][[cell_proportion_option]] <- cell_proportions
      
      ## Convert spe object to data frame
      df <- data.frame(spatialCoords(simulated_spe), "Cell.Type" = simulated_spe[["Cell.Type"]])
      
      ## Just change the cell type of the temp_cell_type, no need to actually re-simulate
      df[["Cell.Type"]] <- ifelse(df[["Cell.Type"]] == temp_cell_type, 
                                  sample(cell_types, size = length(df[["Cell.Type"]]), replace = TRUE, prob = cell_proportions), 
                                  df[["Cell.Type"]])
      
      # Add Cell.ID column to data frame
      df$Cell.ID <- paste("Cell", seq(nrow(df)), sep = "_")
      
      # Update current meta data
      metadata <- simulated_spe@metadata
      metadata[["simulation"]][[length(metadata[["simulation"]])]][[cell_type_option]] <- cell_types
      metadata[["simulation"]][[length(metadata[["simulation"]])]][[cell_proportion_option]] <- cell_proportions
      
      # Convert data frame to spe object
      simulated_spe_new <- SpatialExperiment(
        assay = matrix(data = NA, nrow = nrow(df), ncol = nrow(df)),
        colData = df,
        spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
        metadata = metadata)
      
      ## Generate simulation
      message("Generating simulation...")
      fig <- plot_cells3D(simulated_spe_new)
      print(fig)
      
      if (length(cell_types) == 1) break # If there is only one cell type, proportion is always 1
      
      ## Allow user to re-choose cell proportions  
      message("Would like to re-choose these cell proportions?\n")
      user_input_y_or_n <- get_y_or_n_input()
      if (user_input_y_or_n == "y") {
        cell_proportions <- c()
        max_proportion <- 1
        i <- 1
        message("For each cell type, choose their proportion in the simulation. They must add to 1.\n")
      }
    }
  }
  
  return(list(data = simulated_spe_new, properties = properties))
}


library(alphashape3d)

alpha_hull_clustering3D <- function(spe, 
                                    cell_types_of_interest, 
                                    alpha, 
                                    minimum_cells_in_cluster,
                                    feature_colname = "Cell.Type", 
                                    plot_image = T) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  ## Check cell types of interst are found in the spe object
  unknown_cell_types <- setdiff(cell_types_of_interest, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  if (!(is.numeric(alpha) && length(alpha) == 1 && alpha > 0)) {
    stop("`alpha` is not a positive numeric.")
  }
  if (!(is.integer(minimum_cells_in_cluster) && length(minimum_cells_in_cluster) == 1 || (is.numeric(minimum_cells_in_cluster) && length(minimum_cells_in_cluster) == 1 && minimum_cells_in_cluster > 0 && minimum_cells_in_cluster%%1 == 0))) {
    stop("`minimum_cells_in_cluster` is not a positive integer.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  ## Subset for the chosen cell_types_of_interest
  spe_subset <- spe[ , spe[[feature_colname]] %in% cell_types_of_interest]
  spe_subset_coords <- spatialCoords(spe_subset)
  
  ## Get the alpha hull
  alpha_hull <- ashape3d(as.matrix(spe_subset_coords), alpha = alpha)
  
  if (sum(alpha_hull$triang[, 9]) == 0) stop("alpha value is too small? No alpha hulls identified")
  
  ## Determine which alpha hull cluster each cell_type_of_interest belongs to
  alpha_hull_clusters <- components_ashape3d(alpha_hull)
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(spe), colData(spe))
  
  df_cell_types_of_interest <- df[df[[feature_colname]] %in% cell_types_of_interest, ]
  df_other_cell_types <- df[!(df[[feature_colname]] %in% cell_types_of_interest), ]
  
  df_cell_types_of_interest$alpha_hull_cluster <- alpha_hull_clusters
  df_other_cell_types$alpha_hull_cluster <- 0
  
  ## Ignore cell_types_of_interest which belong to an alpha hull cluster with less than minimum_cells_in_cluster
  alpha_hull_clusters_table <- table(alpha_hull_clusters)
  maximium_alpha_hull_cluster <- Position(function(x) x < minimum_cells_in_cluster, alpha_hull_clusters_table)
  maximium_alpha_hull_cluster <- as.numeric(names(alpha_hull_clusters_table[maximium_alpha_hull_cluster]))
  
  if (!is.na(maximium_alpha_hull_cluster) && maximium_alpha_hull_cluster != -1) {
    spe_subset_coords <- spe_subset_coords[alpha_hull_clusters >= 1 & alpha_hull_clusters < maximium_alpha_hull_cluster, ]
    
    df_cell_types_of_interest$alpha_hull_cluster <- ifelse(alpha_hull_clusters >= 1 & alpha_hull_clusters < maximium_alpha_hull_cluster, 
                                                           alpha_hull_clusters, 0)
    
    ## Get the alpha hull again...
    alpha_hull <- ashape3d(as.matrix(spe_subset_coords), alpha = alpha)
  }
  
  ## Convert data frame to spe object
  df <- rbind(df_cell_types_of_interest, df_other_cell_types)
  
  spe <- SpatialExperiment(
    assay = matrix(data = NA, nrow = nrow(df), ncol = nrow(df)),
    colData = df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = spe@metadata)
  
  ## Get the information of the vertices and faces of the alpha hull (what 3 vertices make up each face triangle?)
  vertices <- alpha_hull$x
  faces <- alpha_hull$triang[alpha_hull$triang[, 9] == 2, c("tr1", "tr2", "tr3")]
  spe@metadata$alpha_hull <- list(vertices = vertices, faces = faces, ashape3d_object = alpha_hull)
  
  ## Plot
  if (plot_image) {
    fig <- plot_alpha_hull_clusters3D(spe, feature_colname = feature_colname)
    methods::show(fig)
  }
  
  return(spe)
}
calculate_all_gradient_cc_metrics3D <- function(spe, 
                                                reference_cell_type, 
                                                target_cell_types, 
                                                radii, 
                                                feature_colname = "Cell.Type", 
                                                plot_image = T) {
  
  # Define constants
  cross_K_df_colnames <- c("reference",
                           "expected",
                           target_cell_types)
  mixing_score_df_colnames <- c("ref_cell_type", 
                                "tar_cell_type", 
                                "n_ref_cells",
                                "n_tar_cells", 
                                "n_ref_tar_interactions",
                                "n_ref_ref_interactions", 
                                "mixing_score", 
                                "normalised_mixing_score")
  cross_G_df_colnames <- c("observed_cross_G",
                           "expected_cross_G")
  co_occurrence_df_colnames <- c("reference",
                                 target_cell_types)
  
  ## Define result
  result <- list("mixing_score" = list(),
                 "cells_in_neighbourhood" = data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types))),
                 "cells_in_neighbourhood_proportion" = data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types))),
                 "entropy" = data.frame(matrix(nrow = length(radii), ncol = 1)),
                 "cross_K" = data.frame(matrix(nrow = length(radii), ncol = length(cross_K_df_colnames))),
                 "cross_L" = data.frame(matrix(nrow = length(radii), ncol = length(cross_K_df_colnames))),
                 "cross_G" = list(),
                 "co_occurrence" = data.frame(matrix(nrow = length(radii), ncol = length(co_occurrence_df_colnames))))
  colnames(result[["cells_in_neighbourhood"]]) <- target_cell_types
  colnames(result[["cells_in_neighbourhood_proportion"]]) <- target_cell_types
  colnames(result[["entropy"]]) <- "entropy"
  colnames(result[["cross_K"]]) <- cross_K_df_colnames
  colnames(result[["cross_L"]]) <- cross_K_df_colnames
  colnames(result[["co_occurrence"]]) <- co_occurrence_df_colnames
  
  # Define individual data frames for mixing_score and cross_G
  for (target_cell_type in target_cell_types) {
    if (reference_cell_type != target_cell_type) {
      result[["mixing_score"]][[target_cell_type]] <- data.frame(matrix(nrow = length(radii), ncol = length(mixing_score_df_colnames)))
      colnames(result[["mixing_score"]][[target_cell_type]]) <- mixing_score_df_colnames
    }
    result[["cross_G"]][[target_cell_type]] <- data.frame(matrix(nrow = length(radii), ncol = length(cross_G_df_colnames)))
    colnames(result[["cross_G"]][[target_cell_type]]) <- cross_G_df_colnames
  }
  
  # Get gradient results for each metric
  for (i in seq(length(radii))) {
    df <- calculate_all_single_radius_cc_metrics3D(spe,
                                                   reference_cell_type,
                                                   target_cell_types,
                                                   radii[i],
                                                   feature_colname)
    
    if (is.null(df)) return(NULL)
    
    df[["cells_in_neighbourhood"]]$ref_cell_id <- NULL
    
    result[["cells_in_neighbourhood"]][i, ] <- apply(df[["cells_in_neighbourhood"]], 2, mean)
    result[["cells_in_neighbourhood_proportion"]][i, ] <- apply(df[["cells_in_neighbourhood_proportion"]][ , paste(target_cell_types, "_prop", sep = "")], 2, mean, na.rm = T)
    result[["entropy"]][i, "entropy"] <- mean(df[["entropy"]]$entropy, na.rm = T)
    result[["cross_K"]][i, ] <- df[["cross_K"]]
    result[["cross_L"]][i, ] <- df[["cross_L"]]
    result[["co_occurrence"]][i, ] <- df[["co_occurrence"]]
    
    for (target_cell_type in names(df[["mixing_score"]])) {
      result[["mixing_score"]][[target_cell_type]][i, ] <- df[["mixing_score"]][[target_cell_type]]
    }
    for (target_cell_type in names(df[["cross_G"]])) {
      result[["cross_G"]][[target_cell_type]][i, ] <- df[["cross_G"]][[target_cell_type]]
    }
  }
  
  # Add radius column to each data frame
  result[["cells_in_neighbourhood"]]$radius <- radii
  result[["cells_in_neighbourhood_proportion"]]$radius <- radii
  result[["entropy"]]$radius <- radii
  result[["cross_K"]]$radius <- radii
  result[["cross_L"]]$radius <- radii
  result[["co_occurrence"]]$radius <- radii
  for (target_cell_type in names(df[["mixing_score"]])) {
    result[["mixing_score"]][[target_cell_type]]$radius <- radii
  }
  for (target_cell_type in names(df[["cross_G"]])) {
    result[["cross_G"]][[target_cell_type]]$radius <- radii
  }
  
  
  ## Plot
  if (plot_image) {
    fig_ACIN <- plot_cells_in_neighbourhood_gradient3D(result[["cells_in_neighbourhood"]], reference_cell_type)
    methods::show(fig_ACIN)
    
    fig_ACINP <- plot_cells_in_neighbourhood_proportions_gradient3D(result[["cells_in_neighbourhood_proportion"]], reference_cell_type)
    methods::show(fig_ACINP)
    
    expected_entropy <- calculate_entropy_background3D(spe, target_cell_types, feature_colname)
    fig_AE <- plot_entropy_gradient3D(result[["entropy"]], expected_entropy, reference_cell_type, target_cell_types)
    methods::show(fig_AE)
    
    for (target_cell_type in names(result[["mixing_score"]])) {
      fig_NMS <- plot_mixing_scores_gradient3D(result[["mixing_score"]][[target_cell_type]], "NMS")
      fig_MS <- plot_mixing_scores_gradient3D(result[["mixing_score"]][[target_cell_type]], "MS")
      fig_NMS_MS <- plot_grid(fig_NMS, fig_MS, nrow = 2)
      methods::show(fig_NMS_MS)
    }
    fig_CK <- plot_cross_K_gradient3D(result[["cross_K"]])
    fig_CKR <- plot_cross_K_gradient_ratio3D(result[["cross_K"]])
    fig_CK_CKR <- plot_grid(fig_CK, fig_CKR, nrow = 2)
    methods::show(fig_CK_CKR)
    
    fig_CL <- plot_cross_L_gradient3D(result[["cross_L"]])
    fig_CLR <- plot_cross_L_gradient_ratio3D(result[["cross_L"]])
    fig_CL_CLR <- plot_grid(fig_CL, fig_CLR, nrow = 2)
    methods::show(fig_CL_CLR)
    
    for (target_cell_type in names(result[["cross_G"]])) {
      fig_CG <- plot_cross_G_gradient3D(result[["cross_G"]][[target_cell_type]], reference_cell_type, target_cell_type)
      methods::show(fig_CG)
    }
    
    fig_co_occ <- plot_co_occurrence_gradient3D(result[["co_occurrence"]])
    methods::show(fig_co_occ)
  }
  
  return(result)
}


### Calculate all single radius cell-colocalisation metrics
# If a function only requires one target cell type, iterate through each cell type in target_cell_types, else use all target_cell_types

calculate_all_single_radius_cc_metrics3D <- function(spe, 
                                                     reference_cell_type, 
                                                     target_cell_types, 
                                                     radius, 
                                                     feature_colname = "Cell.Type") {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (!(is.character(reference_cell_type) && length(reference_cell_type) == 1)) {
    stop("`reference_cell_type` is not a character.")
  }
  if (!is.character(target_cell_types)) {
    stop("`target_cell_types` is not a character vector.")
  }
  if (!(is.numeric(radius) && length(radius) == 1 && radius > 0)) {
    stop(paste(radius, " is not a positive numeric."))
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  
  ## For reference_cell_type, check it is found in the spe object
  if (!(reference_cell_type %in% spe[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spe object"))
    return(NULL)
  }
  ## For target_cell_types, check they are found in the spe object
  unknown_cell_types <- setdiff(target_cell_types, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
  }
  
  # Define result
  result <- list("cells_in_neighbourhood" = list(),
                 "cells_in_neighbourhood_proportion" = list(),
                 "entropy" = list(),
                 "mixing_score" = list(),
                 "cross_K" = list(),
                 "cross_L" = list(),
                 "cross_G" = list(),
                 "co_occurrence" = list())
  
  # Define other constants
  mixing_score_df_colnames <- c("ref_cell_type", 
                                "tar_cell_type", 
                                "n_ref_cells",
                                "n_tar_cells", 
                                "n_ref_tar_interactions",
                                "n_ref_ref_interactions", 
                                "mixing_score", 
                                "normalised_mixing_score")
  
  cross_K_df_colnames <- c("reference",
                           "expected",
                           target_cell_types)
  cross_G_df_colnames <- c("observed_cross_G",
                           "expected_cross_G")
  co_occurrence_df_colnames <- c("reference",
                                 target_cell_types)
  
  # Get rough dimensions of window for cross_K
  spe_coords <- data.frame(spatialCoords(spe))
  length <- round(max(spe_coords$Cell.X.Position) - min(spe_coords$Cell.X.Position))
  width  <- round(max(spe_coords$Cell.Y.Position) - min(spe_coords$Cell.Y.Position))
  height <- round(max(spe_coords$Cell.Z.Position) - min(spe_coords$Cell.Z.Position))
  ## Get volume of the window the cells are in
  volume <- length * width * height
  
  
  
  # All single radius cc metrics stem from calculate_entropy3D function
  entropy_df <- calculate_entropy3D(spe, 
                                    reference_cell_type, 
                                    target_cell_types, 
                                    radius, 
                                    feature_colname)  
  
  ## Cells in neighbourhood ----------
  result[["cells_in_neighbourhood"]] <- entropy_df[ , c("ref_cell_id", target_cell_types)]
  
  ## Cells in neighbourhood proportion ----------
  result[["cells_in_neighbourhood_proportion"]] <- entropy_df[ , c("ref_cell_id", target_cell_types, paste(target_cell_types, "_prop", sep = ""))]
  
  ## Entropy --------------
  result[["entropy"]] <- entropy_df
  
  ## Mixing score -----------------
  for (target_cell_type in target_cell_types) {
    mixing_score_df <- data.frame(matrix(nrow = 1, ncol = length(mixing_score_df_colnames)))
    colnames(mixing_score_df) <- mixing_score_df_colnames
    mixing_score_df$ref_cell_type <- reference_cell_type
    
    # No need to fill in mixing_score_df if the reference and target cell is the same
    if (reference_cell_type != target_cell_type) {
      mixing_score_df$tar_cell_type <- target_cell_type
      mixing_score_df$n_ref_cells <- sum(spe[[feature_colname]] == reference_cell_type)
      mixing_score_df$n_tar_cells <- sum(spe[[feature_colname]] == target_cell_type)
      mixing_score_df$n_ref_tar_interactions <- sum(entropy_df[[target_cell_type]])
      mixing_score_df$n_ref_ref_interactions <- sum(entropy_df[[reference_cell_type]])
      mixing_score_df$mixing_score <- mixing_score_df$n_ref_tar_interactions / (0.5 * mixing_score_df$n_ref_ref_interactions)
      mixing_score_df$normalised_mixing_score <- 0.5 * mixing_score_df$mixing_score * mixing_score_df$n_ref_cells / mixing_score_df$n_tar_cell
      if (is.infinite(mixing_score_df$mixing_score)) mixing_score_df$mixing_score <- NA
      if (is.infinite(mixing_score_df$normalised_mixing_score)) mixing_score_df$normalised_mixing_score <- NA
      result[["mixing_score"]][[target_cell_type]] <- mixing_score_df
    }
  }
  
  ## Cross_K ---------------------
  cross_K_df <- data.frame(matrix(nrow = 1, ncol = length(cross_K_df_colnames)))
  colnames(cross_K_df) <- cross_K_df_colnames
  cross_K_df$reference <- reference_cell_type
  cross_K_df$expected <- (4/3) * pi * radius^3
  
  for (target_cell_type in target_cell_types) {
    cross_K_df[[target_cell_type]] <- (((volume * sum(entropy_df[[target_cell_type]])) / sum(spe[[feature_colname]] == reference_cell_type)) / sum(spe[[feature_colname]] == target_cell_type)) 
  }
  result[["cross_K"]] <- cross_K_df
  
  ## Cross_L ---------------------
  cross_L_df <- cross_K_df
  cross_L_df[ , c("expected", target_cell_types)] <- (cross_L_df[ , c("expected", target_cell_types)] / (4 * pi / 3)) ^ (1/3)
  result[["cross_L"]] <- cross_L_df
  
  ## Cross_G ---------------------
  for (target_cell_type in target_cell_types) {
    cross_G_df <- data.frame(matrix(nrow = 1, ncol = length(cross_G_df_colnames)))
    colnames(cross_G_df) <- cross_G_df_colnames
    
    reference_target_interactions <- entropy_df[[target_cell_type]]
    n_target_cells <- sum(spe[[feature_colname]] == target_cell_type)
    target_cell_type_intensity <- n_target_cells / volume
    observed_cross_G <- sum(reference_target_interactions != 0) / length(reference_target_interactions)
    expected_cross_G <- 1 - exp(-1 * target_cell_type_intensity * (4 / 3) * pi * radius^3)
    
    cross_G_df$observed_cross_G <- observed_cross_G
    cross_G_df$expected_cross_G <- expected_cross_G
    result[["cross_G"]][[target_cell_type]] <- cross_G_df
  }
  
  
  ## Co_occurrence ---------------
  all_cell_types <- unique(spe[[feature_colname]])
  cells_in_neighbourhood_proportions_df <- calculate_cells_in_neighbourhood_proportions3D(spe,
                                                                                          reference_cell_type,
                                                                                          all_cell_types,
                                                                                          radius,
                                                                                          feature_colname)
  
  co_occurrence_df <- data.frame(matrix(nrow = 1, ncol = length(co_occurrence_df_colnames)))
  colnames(co_occurrence_df) <- co_occurrence_df_colnames
  co_occurrence_df$reference <- reference_cell_type
  
  n_cells_in_spe <- length(spe[[feature_colname]])
  n_cells_in_reference_cell_type_radius <- sum(cells_in_neighbourhood_proportions_df$total)
  
  for (target_cell_type in target_cell_types) {
    n_target_cells_in_reference_cell_type_radius <- sum(cells_in_neighbourhood_proportions_df[[target_cell_type]])
    target_cell_type_proportion_in_reference_cell_type_radius <- n_target_cells_in_reference_cell_type_radius / n_cells_in_reference_cell_type_radius
    n_target_cells_in_spe <- sum(spe[[feature_colname]] == target_cell_type)
    target_cell_type_proportion_in_spe <- n_target_cells_in_spe / n_cells_in_spe
    target_cell_type_co_occurrence <- target_cell_type_proportion_in_reference_cell_type_radius / target_cell_type_proportion_in_spe
    
    co_occurrence_df[[target_cell_type]] <- target_cell_type_co_occurrence
  }
  result[["co_occurrence"]] <- co_occurrence_df
  
  return(result)
}

calculate_border_of_clusters3D <- function(spe, 
                                           radius,
                                           cluster_colname, 
                                           feature_colname = "Cell.Type", 
                                           plot_image = T) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (!(is.numeric(radius) && length(radius) == 1 && radius > 0)) {
    stop("`radius` is not a positive numeric.")
  }
  if (!is.character(cluster_colname)) {
    stop("`cluster_colname` is not a character. This should be 'alpha_hull_cluster', 'dbscan_cluster', or 'grid_based_cluster', depending on the chosen method.")
  }
  if (is.null(spe[[cluster_colname]])) {
    stop(paste("No column called", cluster_colname, "found in spe object."))
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  ## Get spatial coords of spe
  spe_coords <- data.frame(spatialCoords(spe))
  
  ## Get coords of non-cluster cells
  non_cluster_coords <- spe_coords[spe[[cluster_colname]] == 0, ]
  
  # New column for spe object: 'cluster_border'. Default is 'outside'
  spe$cluster_border <- "outside"
  
  # Label cells part of a cluster (e.g. 'cluster1')
  spe$cluster_border[spe[[cluster_colname]] != 0] <- paste("inside_C", spe[[cluster_colname]][spe[[cluster_colname]] != 0], sep = "")
  
  ## Iterate for each cluster
  n_clusters <- max(spe[[cluster_colname]])
  
  for (i in seq_len(n_clusters)) {
    
    ## Subset for cells in the current cluster of interest
    cluster_coords <- spe_coords[spe[[cluster_colname]] == i, ]
    
    # For each cell in the current cluster, check how many other cells in the cluster are in its radius
    cluster_to_cluster_interactions <- dbscan::frNN(cluster_coords, radius)
    
    # Determine the median minimum number of cluster cells found in the radius of cluster cell. Use this as the threshold for non-cluster cells.
    non_cluster_threshold <- quantile(unlist(lapply(cluster_to_cluster_interactions$dist, length)), 0.5)
    
    # For each non-cluster cell, check how many cluster cells are in its radius.
    non_cluster_to_cluster_interactions <- dbscan::frNN(cluster_coords, radius, non_cluster_coords)
    
    # If number of cluster cells found in the radius of non-cluster cells is greater than threshold, non-cluster cell has probably infiltrated cluster too
    n_cluster_cells_in_non_cluster_cell_radius <- unlist(lapply(non_cluster_to_cluster_interactions$id, length))
    
    spe$cluster_border[as.numeric(names(non_cluster_to_cluster_interactions$id)[n_cluster_cells_in_non_cluster_cell_radius > non_cluster_threshold])] <- paste("infiltrated_C", i, sep = "")
    
    # If number of cluster cells found in the radius of non-cluster cells is less than threshold, but greater than 0, non-cluster cell is probably on the border
    spe$cluster_border[as.numeric(names(non_cluster_to_cluster_interactions$id)[n_cluster_cells_in_non_cluster_cell_radius > 0 & n_cluster_cells_in_non_cluster_cell_radius < non_cluster_threshold])] <- paste("border_C", i, sep = "")
  }
  
  ## Plot
  if (plot_image) {
    fig <- plot_cells3D(spe, feature_colname = "cluster_border")
    methods::show(fig)
  }
  
  return(spe)
}
calculate_cell_proportion_grid_metrics3D <- function(spe, 
                                                     n_splits,
                                                     reference_cell_types,
                                                     target_cell_types,
                                                     feature_colname = "Cell.Type",
                                                     plot_image = TRUE) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  ## Check reference_cell_types are found in the spe object
  unknown_cell_types <- setdiff(reference_cell_types, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in reference_cell_types are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  ## Check target_cell_types are found in the spe object
  unknown_cell_types <- setdiff(target_cell_types, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  # Check if there is intersection between reference_cell_types and target_cell_types
  if (length(intersect(reference_cell_types, target_cell_types)) > 0) {
    stop("Cannot have same cells in both reference_cell_types and target_cell_types")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Add grid metrics to spe
  spe <- get_spe_grid_metrics3D(spe, n_splits, feature_colname)
  
  # Get grid_prism_cell_matrix from spe
  grid_prism_cell_matrix <- spe@metadata$grid_metrics$grid_prism_cell_matrix
  
  ## Define data frame which contains all results
  n_grid_prisms <- n_splits^3
  result <- data.frame(row.names = seq(n_grid_prisms))
  
  # Fill in the result data frame
  if (length(reference_cell_types) == 1) {
    result$reference <- grid_prism_cell_matrix[[reference_cell_types]]
  }
  else {
    result$reference <- rowSums(grid_prism_cell_matrix[ , reference_cell_types])
  }
  if (length(target_cell_types) == 1) {
    result$target <- grid_prism_cell_matrix[[target_cell_types]]
  }
  else {
    result$target <- rowSums(grid_prism_cell_matrix[ , target_cell_types])
  }
  result$total <- result$reference + result$target
  result$proportion <- result$target / result$total
  
  # Add grid_prism_coordinates info to result
  result <- cbind(result, spe@metadata$grid_metrics$grid_prism_coordinates)
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_metrics_continuous3D(result, "proportion")
    methods::show(fig)
  }
  
  return(result)
}
calculate_cell_proportions_of_clusters3D <- function(spe, cluster_colname, feature_colname = "Cell.Type", plot_image = T) {
  
  # Get number of clusters
  n_clusters <- max(spe[[cluster_colname]])
  
  ## Get different cell types found in the clusters (alphabetical for consistency)
  cell_types <- unique(spe[[feature_colname]][spe[[cluster_colname]] != 0])
  cell_types <- cell_types[order(cell_types)]
  
  ## For each cluster, determine the size and cell proportion of each cluster
  result <- data.frame(matrix(nrow = n_clusters, ncol = 2 + length(cell_types)))
  colnames(result) <- c("cluster_number", "n_cells", cell_types)
  result$cluster_number <- as.character(seq(n_clusters))
  
  for (i in seq(n_clusters)) {
    cells_in_cluster <- spe[[feature_colname]][spe[[cluster_colname]] == i]
    result[i, "n_cells"] <- length(cells_in_cluster)
    
    for (cell_type in cell_types) {
      result[i, cell_type] <- sum(cells_in_cluster == cell_type) / result[i, "n_cells"]
    }
  }
  
  ## Plot
  if (plot_image) {
    plot_result <- reshape2::melt(result, id.vars = c("cluster_number", "n_cells"))
    fig <- ggplot(plot_result, aes(cluster_number, value, fill = variable)) +
      geom_bar(stat = "identity") +
      labs(title = "Cell proportions of each cluster", x = "", y = "Cell proportion") +
      scale_x_discrete(labels = paste("cluster_", result$cluster_number, ", n = ", result$n_cells, sep = "")) +
      guides(fill = guide_legend(title="Cell type")) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
    
    methods::show(fig)
  }
  
  return(result)
}

#' @title Calculate cell proportions in 3D spatial data.
#'
#' @description This function calculates the proportions of different cell types in a 3D SpatialExperiment object. 
#'    It can optionally plot a bar chart of the cell proportions.
#'
#' @param spe A SpatialExperiment object containing 3D spatial information for the cells.
#' @param cell_types_of_interest A character vector specifying the cell types of interest.
#'    If NULL, all cell types in the `feature_colname` column will be considered.
#' @param feature_colname A string specifying the name of the column in the `colData` slot of the SpatialExperiment
#'    object that contains the cell type information.
#' @param plot_image A logical indicating whether to plot violin plots of the minimum distances 
#'    between cell type pairs. Defaults to TRUE.
#'
#' @return A data frame containing the cell types, their frequencies, proportions, and percentages.
#'
#' @examples
#' cell_proportions <- calculate_cell_proportions3D(
#'     spe = SPIAT3D::simulated_spe, 
#'     cell_types_of_interest = NULL, 
#'     feature_colname = "Cell.Type", 
#'     plot_image = TRUE
#' )
#' 
#' @export


calculate_cell_proportions3D <- function(spe,
                                         cell_types_of_interest = NULL, 
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (ncol(spe) == 0) {
    stop("No cells found for calculating cell proportions.")
  }
  if (!(is.null(cell_types_of_interest) || is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Creates frequency/bar plot of all cell types in the entire image
  cell_proportions <- data.frame(table(spe[[feature_colname]]))
  names(cell_proportions) <- c("cell_type", 'frequency')
  
  # Only include cell types the user has chosen
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spe object
    unknown_cell_types <- setdiff(cell_types_of_interest, cell_proportions$cell_type)
    if (length(unknown_cell_types) != 0) {
      stop(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                 paste(unknown_cell_types, collapse = ", ")))
    }
    
    # Subset for cell types chosen by user
    cell_proportions <- cell_proportions[(cell_proportions$cell_type %in% cell_types_of_interest), ]
    
  }
  
  # Get frequency total for all cells
  cell_type_frequency_total <- sum(cell_proportions$frequency)
  
  # Get proportions and percentages
  cell_proportions$proportion <- cell_proportions$frequency / cell_type_frequency_total
  cell_proportions$percentage <- cell_proportions$proportion * 100
  
  # Order the cell types by proportion (highest cell proportion is first)
  cell_proportions <- cell_proportions[rev(order(cell_proportions$proportion)), ]
  rownames(cell_proportions) <- seq(nrow(cell_proportions))
  
  # Plot
  if (plot_image) {
    
    labels <- paste(round(cell_proportions$percentage, 1), "%", sep = "")
    
    fig <- ggplot(cell_proportions, aes(x = factor(cell_type, cell_type), y = percentage, fill = cell_type)) +
      geom_bar(stat='identity') + 
      theme_bw() +
      labs(title="Cell proportions", x = "Cell type", y = "Percentage") +
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "none") +
      geom_text(aes(label = labels), vjust = 0)
    
    methods::show(fig)
  }
  
  return(cell_proportions)
}
calculate_cells_in_neighbourhood_gradient3D <- function(spe, 
                                                        reference_cell_type, 
                                                        target_cell_types, 
                                                        radii, 
                                                        feature_colname = "Cell.Type",
                                                        plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types)))
  colnames(result) <- target_cell_types
  
  for (i in seq(length(radii))) {
    cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood3D(spe,
                                                                    reference_cell_type,
                                                                    target_cell_types,
                                                                    radii[i],
                                                                    feature_colname,
                                                                    FALSE,
                                                                    FALSE)
    
    if (is.null(cells_in_neighbourhood_df)) return(NULL)
    
    cells_in_neighbourhood_df$ref_cell_id <- NULL
    result[i, ] <- apply(cells_in_neighbourhood_df, 2, mean)
  }
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig <- plot_cells_in_neighbourhood_gradient3D(result, reference_cell_type)
    methods::show(fig)
  }
  
  return(result)
}

calculate_cells_in_neighbourhood_proportions_gradient3D <- function(spe, 
                                                                    reference_cell_type, 
                                                                    target_cell_types, 
                                                                    radii, 
                                                                    feature_colname = "Cell.Type",
                                                                    plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types)))
  colnames(result) <- target_cell_types
  
  for (i in seq(length(radii))) {
    cell_proportions_neighbourhood_proportions_df <- calculate_cells_in_neighbourhood_proportions3D(spe,
                                                                                                    reference_cell_type,
                                                                                                    target_cell_types,
                                                                                                    radii[i],
                                                                                                    feature_colname)
    
    if (is.null(cell_proportions_neighbourhood_proportions_df)) return(NULL)
    
    result[i, ] <- apply(cell_proportions_neighbourhood_proportions_df[ , paste(target_cell_types, "_prop", sep = "")], 2, mean, na.rm = T)
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  # Plot
  if (plot_image) {
    fig <- plot_cells_in_neighbourhood_proportions_gradient3D(result, reference_cell_type)
    methods::show(fig)
  }
  
  return(result)
}

calculate_cells_in_neighbourhood_proportions3D <- function(spe, 
                                                           reference_cell_type, 
                                                           target_cell_types, 
                                                           radius, 
                                                           feature_colname = "Cell.Type") {
  
  ## Get cells in neighbourhood df
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood3D(spe,
                                                                  reference_cell_type,
                                                                  target_cell_types,
                                                                  radius,
                                                                  feature_colname,
                                                                  FALSE,
                                                                  FALSE)
  
  if (is.null(cells_in_neighbourhood_df)) return(NULL)
  
  ## Get total number of target cells for each row (first column is the reference cell id column, so we exclude it)
  cells_in_neighbourhood_df$total <- apply(cells_in_neighbourhood_df[ , c(-1)], 1, sum)
  
  cells_in_neighbourhood_df[ , paste(target_cell_types, "_prop", sep = "")] <- cells_in_neighbourhood_df[ , target_cell_types] / cells_in_neighbourhood_df$total
  
  return(cells_in_neighbourhood_df)
}
calculate_cells_in_neighbourhood3D <- function(spe, 
                                               reference_cell_type, 
                                               target_cell_types, 
                                               radius, 
                                               feature_colname = "Cell.Type",
                                               show_summary = TRUE,
                                               plot_image = TRUE) {
  
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (!(is.character(reference_cell_type) && length(reference_cell_type) == 1)) {
    stop("`reference_cell_type` is not a character.")
  }
  if (!is.character(target_cell_types)) {
    stop("`target_cell_types` is not a character vector.")
  }
  if (!(is.numeric(radius) && length(radius) == 1 && radius > 0)) {
    stop("`radius` is not a positive numeric.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(show_summary)) {
    stop("`show_summary` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  ## For reference_cell_type, check it is found in the spe object
  if (!(reference_cell_type %in% spe[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spe object"))
    return(NULL)
  }
  ## For target_cell_types, check they are found in the spe object
  unknown_cell_types <- setdiff(target_cell_types, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
  }
  
  if (is.null(spe[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spe")
    spe$Cell.ID <- paste("Cell", seq(ncol(spe)), sep = "_")
  }  
  
  # Get spe coords
  spe_coords <- data.frame(spatialCoords(spe))
  
  # Get reference_cell_type coords
  reference_cell_type_coords <- spe_coords[spe[[feature_colname]] == reference_cell_type, ]
  
  result <- data.frame(matrix(nrow = nrow(reference_cell_type_coords), ncol = 0))
  
  for (target_cell_type in target_cell_types) {
    
    if (sum(spe[[feature_colname]] == target_cell_type) == 0) {
      result[[target_cell_type]] <- NA
      next
    }
    
    ## Get target_cell_type coords
    target_cell_type_coords <- spe_coords[spe[[feature_colname]] == target_cell_type, ]
    
    ## Determine number of target cells specified distance for each reference cell
    ref_tar_result <- dbscan::frNN(target_cell_type_coords, 
                                   eps = radius,
                                   query = reference_cell_type_coords, 
                                   sort = FALSE)
    
    n_targets <- rapply(ref_tar_result$id, length)
    
    
    # Don't want to include the reference cell as one of the target cells
    if (reference_cell_type == target_cell_type) n_targets <- n_targets - 1
    
    ## Add to data frame
    result[[target_cell_type]] <- n_targets
  }
  
  result <- data.frame(ref_cell_id = spe$Cell.ID[spe[[feature_colname]] == reference_cell_type], result)
  
  ## Print summary
  if (show_summary) {
    print(summarise_cells_in_neighbourhood3D(result))    
  }
  
  ## Plot
  if (plot_image) {
    fig <- plot_cells_in_neighbourhood_violin3D(result, reference_cell_type)
    methods::show(fig)
  }
  
  return(result)
}
### Assume that clusters have uniform density and that the centre of each cluster is defined by its centre of mass
### Centre of mass can be estimated by taking the average of the x, y, and z coordinates of cells in the cluster

calculate_center_of_clusters3D <- function(spe, cluster_colname) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (!is.character(cluster_colname)) {
    stop("`cluster_colname` is not a character. This should be 'alpha_hull_cluster', 'dbscan_cluster', or 'grid_based_cluster', depending on the chosen method.")
  }
  if (is.null(spe[[cluster_colname]])) {
    stop(paste("No column called", cluster_colname, "found in spe object."))
  }
  
  # Get number of clusters
  n_clusters <- max(spe[[cluster_colname]])
  
  # Get spe coords
  spe_coords <- spatialCoords(spe)
  
  ## For each cluster, determine the number of cells in each cluster of each cluster
  result <- data.frame(matrix(nrow = n_clusters, ncol = 4))
  colnames(result) <- c("cluster_number", "Centre.X.Position", "Centre.Y.Position", "Centre.Z.Position")
  
  result$cluster_number <- as.character(seq(n_clusters))
  for (i in seq(n_clusters)) {
    spe_cluster_coords <- spe_coords[spe[[cluster_colname]] == i, ]
    result[i, c("Centre.X.Position", "Centre.Y.Position", "Centre.Z.Position")] <- 
      apply(spe_cluster_coords, 2, mean)
  }
  
  return(result)
}
calculate_co_occurrence_gradient3D <- function(spe, 
                                               reference_cell_type, 
                                               target_cell_types, 
                                               radii, 
                                               feature_colname = "Cell.Type",
                                               plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types) + 1))
  colnames(result) <- c("reference", target_cell_types)
  
  for (i in seq(length(radii))) {
    co_occurrence_df <- calculate_co_occurrence3D(spe,
                                                  reference_cell_type,
                                                  target_cell_types,
                                                  radii[i],
                                                  feature_colname)
    
    result[i, ] <- co_occurrence_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig <- plot_co_occurrence_gradient3D(result)
    methods::show(fig)
  }
  
  return(result)
}
calculate_co_occurrence3D <- function(spe, 
                                      reference_cell_type, 
                                      target_cell_types, 
                                      radius, 
                                      feature_colname = "Cell.Type") {
  
  # Get all cell types in spe
  all_cell_types <- unique(spe[[feature_colname]])
  
  cells_in_neighbourhood_proportions_df <- calculate_cells_in_neighbourhood_proportions3D(spe,
                                                                                          reference_cell_type,
                                                                                          all_cell_types,
                                                                                          radius,
                                                                                          feature_colname)
  
  result <- data.frame(reference = reference_cell_type)
  
  # Get total number of cells in spe
  n_cells_in_spe <- length(spe[[feature_colname]])
  
  # Get total number of cells in radius around reference cell type
  n_cells_in_reference_cell_type_radius <- sum(cells_in_neighbourhood_proportions_df$total)
  
  for (target_cell_type in target_cell_types) {
    
    # Get total number of target cells in radius around reference cell type
    n_target_cells_in_reference_cell_type_radius <- sum(cells_in_neighbourhood_proportions_df[[target_cell_type]])
    
    # Get proportion of target cells in radius around reference cell type
    target_cell_type_proportion_in_reference_cell_type_radius <- n_target_cells_in_reference_cell_type_radius / n_cells_in_reference_cell_type_radius
    
    # Get proportion of target cell type in spe
    n_target_cells_in_spe <- sum(spe[[feature_colname]] == target_cell_type)
    target_cell_type_proportion_in_spe <- n_target_cells_in_spe / n_cells_in_spe
    
    # Get co-occurence value for taget cell type
    target_cell_type_co_occurrence <- target_cell_type_proportion_in_reference_cell_type_radius / target_cell_type_proportion_in_spe
    
    # Add to result data frame
    result[[target_cell_type]] <- target_cell_type_co_occurrence
  }
  
  return(result)
}
calculate_cross_G_gradient3D <- function(spe, 
                                         reference_cell_type, 
                                         target_cell_type, 
                                         radii, 
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = 2))
  colnames(result) <- c("observed_cross_G", 
                        "expected_cross_G")
  
  for (i in seq(length(radii))) {
    cross_G_df <- calculate_cross_G3D(spe,
                                      reference_cell_type,
                                      target_cell_type,
                                      radii[i],
                                      feature_colname)
    
    result[i, ] <- cross_G_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig <- plot_cross_G_gradient3D(result, reference_cell_type, target_cell_type)
    methods::show(fig)
  }
  
  return(result)
}
calculate_cross_G3D <- function(spe,
                                reference_cell_type,
                                target_cell_type,
                                radius,
                                feature_colname = "Cell.Type") {
  
  ### Calculate the observed cross_G
  # Get the number of target cells in the radius around each reference cell
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood3D(spe,
                                                                  reference_cell_type,
                                                                  target_cell_type,
                                                                  radius,
                                                                  feature_colname,
                                                                  show_summary = FALSE,
                                                                  plot_image = FALSE)
  
  reference_target_interactions <- cells_in_neighbourhood_df[[target_cell_type]]
  
  # cross_G: essentially the proportion of reference cells with at least 1 target cell within the chosen radius.
  observed_cross_G <- sum(reference_target_interactions != 0) / length(reference_target_interactions)
  
  ### Calculate the expected cross_G
  # Get rough dimensions of the window the points are in
  spe_coords <- data.frame(spatialCoords(spe))
  
  length <- round(max(spe_coords$Cell.X.Position) - min(spe_coords$Cell.X.Position))
  width  <- round(max(spe_coords$Cell.Y.Position) - min(spe_coords$Cell.Y.Position))
  height <- round(max(spe_coords$Cell.Z.Position) - min(spe_coords$Cell.Z.Position))
  
  # Get volume of the window the cells are in
  volume <- length * width * height
  
  # Get the number of target cells
  n_target_cells <- sum(spe[[feature_colname]] == target_cell_type)
  
  # Get target_cell_type intensity (density)
  target_cell_type_intensity <- n_target_cells / volume
  
  # Apply formula
  expected_cross_G <- 1 - exp(-1 * target_cell_type_intensity * (4 / 3) * pi * radius^3)
  
  result <- data.frame(observed_cross_G = observed_cross_G,
                       expected_cross_G = expected_cross_G)
  
  return(result)
}
calculate_cross_K_gradient3D <- function(spe, 
                                         reference_cell_type, 
                                         target_cell_types, 
                                         radii, 
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = 2 + length(target_cell_types)))
  colnames(result) <- c("reference", "expected", target_cell_types)
  
  for (i in seq(length(radii))) {
    cross_K_df <- calculate_cross_K3D(spe,
                                      reference_cell_type,
                                      target_cell_types,
                                      radii[i],
                                      feature_colname)
    
    result[i, ] <- cross_K_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig1 <- plot_cross_K_gradient3D(result)
    fig2 <- plot_cross_K_gradient_ratio3D(result)
    
    combined_fig <- plot_grid(fig1, fig2, nrow = 2)
    methods::show(combined_fig)
  }
  
  return(result)
}
calculate_cross_K3D <- function(spe, 
                                reference_cell_type, 
                                target_cell_types, 
                                radius, 
                                feature_colname = "Cell.Type") {
  
  if (is.null(spe[[feature_colname]])) stop(paste("No column called", feature_colname, "found in spe object"))
  
  if (is.null(spe[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spe")
    spe$Cell.ID <- paste("Cell", seq(ncol(spe)), sep = "_")
  }  
  
  
  ## Get expected cross K-function
  expected_cross_K <- (4/3) * pi * radius^3
  
  ## For reference_cell_type, check it is found in the spe object
  if (!(reference_cell_type %in% spe[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spe object"))
    result <- data.frame(observed_cross_K = NA,
                         expected_cross_K = expected_cross_K,
                         cross_K_ratio = NA)
    return(result)
  }
  
  ## Get rough dimensions of the window the points are in
  spe_coords <- data.frame(spatialCoords(spe))
  
  length <- round(max(spe_coords$Cell.X.Position) - min(spe_coords$Cell.X.Position))
  width  <- round(max(spe_coords$Cell.Y.Position) - min(spe_coords$Cell.Y.Position))
  height <- round(max(spe_coords$Cell.Z.Position) - min(spe_coords$Cell.Z.Position))
  ## Get volume of the window the cells are in
  volume <- length * width * height
  
  # Number of reference cell types is constant
  n_ref_cells <- sum(spe[[feature_colname]] == reference_cell_type)
  
  # Define result data frame
  result <- data.frame(reference = reference_cell_type, expected = expected_cross_K)
  
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood3D(spe,
                                                                  reference_cell_type,
                                                                  target_cell_types,
                                                                  radius,
                                                                  feature_colname,
                                                                  show_summary = FALSE,
                                                                  plot_image = FALSE)
  
  for (target_cell_type in target_cell_types) {
    
    n_ref_tar_interactions <- sum(cells_in_neighbourhood_df[[target_cell_type]])
    
    n_tar_cells <- sum(spe[[feature_colname]] == target_cell_type)
    
    ## Get observed cross K-function
    if (n_tar_cells == 0) {
      observed_cross_K <- NA
    }
    else {
      observed_cross_K <- (volume * n_ref_tar_interactions) / (n_ref_cells * n_tar_cells)  
    }
    result[[target_cell_type]] <- observed_cross_K
  }
  
  return(result)
}
calculate_cross_L_gradient3D <- function(spe, 
                                         reference_cell_type, 
                                         target_cell_types, 
                                         radii, 
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = 2 + length(target_cell_types)))
  colnames(result) <- c("reference", "expected", target_cell_types)
  
  for (i in seq(length(radii))) {
    cross_L_df <- calculate_cross_L3D(spe,
                                      reference_cell_type,
                                      target_cell_types,
                                      radii[i],
                                      feature_colname)
    
    result[i, ] <- cross_L_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig1 <- plot_cross_L_gradient3D(result)
    fig2 <- plot_cross_L_gradient_ratio3D(result)
    
    combined_fig <- plot_grid(fig1, fig2, nrow = 2)
    methods::show(combined_fig)
  }
  
  return(result)
}
calculate_cross_L3D <- function(spe, 
                                reference_cell_type, 
                                target_cell_types, 
                                radius, 
                                feature_colname = "Cell.Type") {
  
  result <- calculate_cross_K3D(spe = spe,
                                reference_cell_type = reference_cell_type,
                                target_cell_types = target_cell_types,
                                radius = radius,
                                feature_colname = feature_colname)
  
  result[ , c("expected", target_cell_types)] <- (result[ , c("expected", target_cell_types)] / (4 * pi / 3)) ^ (1/3)
  
  return(result)
}
calculate_entropy_background3D <- function(spe,
                                           cell_types_of_interest, 
                                           feature_colname = "Cell.Type") {
  
  # NULL case: entropy is undefined
  if (is.null(cell_types_of_interest)) return(NA)
  
  # One cell type case: entropy is 0
  if (is.character(cell_types_of_interest) && length(cell_types_of_interest) == 1) return(0)
  
  cell_proportions_data <- calculate_cell_proportions3D(spe, cell_types_of_interest, feature_colname, FALSE)
  
  # Calculate entropy of the entire image
  entropy <- -1 * sum(cell_proportions_data$proportion * log(cell_proportions_data$proportion, length(cell_proportions_data$proportion)))
  
  return(entropy) 
}

calculate_entropy_gradient3D <- function(spe,
                                         reference_cell_type,
                                         target_cell_types,
                                         radii,
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = 1))
  colnames(result) <- "entropy"
  
  for (i in seq(length(radii))) {
    entropy_df <- calculate_entropy3D(spe,
                                      reference_cell_type,
                                      target_cell_types,
                                      radii[i],
                                      feature_colname)
    
    if (is.null(entropy_df)) return(NULL)
    
    result[i, "entropy"] <- mean(entropy_df$entropy, na.rm = T)
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    expected_entropy <- calculate_entropy_background3D(spe, target_cell_types, feature_colname)
    fig <- plot_entropy_gradient3D(result, expected_entropy, reference_cell_type, target_cell_types)
    methods::show(fig)
  }
  
  return(result)
}

calculate_entropy_grid_metrics3D <- function(spe, 
                                             n_splits,
                                             cell_types_of_interest,
                                             feature_colname = "Cell.Type",
                                             plot_image = TRUE) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  ## Check cell_types_of_interest are found in the spe object
  unknown_cell_types <- setdiff(cell_types_of_interest, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  ## If cell types have been chosen, check they are found in the spe object
  unknown_cell_types <- setdiff(cell_types_of_interest, unique(spe[[feature_colname]]))
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Add grid metrics to spe
  spe <- get_spe_grid_metrics3D(spe, n_splits, feature_colname)
  
  # Get grid_prism_cell_matrix from spe
  grid_prism_cell_matrix <- spe@metadata$grid_metrics$grid_prism_cell_matrix
  
  ## Define data frame which contains all results
  n_grid_prisms <- n_splits^3
  result <- data.frame(row.names = seq(n_grid_prisms))
  
  for (cell_type in cell_types_of_interest) {
    result[[cell_type]] <- grid_prism_cell_matrix[[cell_type]]
  }
  result$total <- rowSums(result)
  
  ## Get data frame containing proportions for cell_types_of_interest
  df_props <- result[ , cell_types_of_interest] / result$total
  
  ## Use proportion data frame to get entropy
  calculate_entropy <- function(x) {
    entropy <- -1 * sum(x * ifelse(is.infinite(log(x, length(x))), 0, log(x, length(x))))
    return(entropy)
  }
  result$entropy <- apply(df_props, 1, calculate_entropy)
  
  # Add grid_prism_coordinates info to result
  result <- cbind(result, spe@metadata$grid_metrics$grid_prism_coordinates)
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_metrics_continuous3D(result, "entropy")
    methods::show(fig)
  }
  
  return(result)
}
calculate_entropy3D <- function(spe,
                                reference_cell_type,
                                target_cell_types,
                                radius,
                                feature_colname = "Cell.Type") {
  
  # Check target_cell_types
  if (!(is.character(target_cell_types) && length(target_cell_types) >= 2)) {
    stop("`target_cell_types` is not a character vector with at least 2 cell types.")
  }
  
  ## Users should ensure include the reference_cell_type as one of the target_cell_types
  cells_in_neighbourhood_proportion_df <- calculate_cells_in_neighbourhood_proportions3D(spe,
                                                                                         reference_cell_type,
                                                                                         target_cell_types,
                                                                                         radius,
                                                                                         feature_colname)
  
  if (is.null(cells_in_neighbourhood_proportion_df)) return(NULL)
  
  ## Get entropy for each row
  cells_in_neighbourhood_proportion_df$entropy <- apply(cells_in_neighbourhood_proportion_df[ , paste(target_cell_types, "_prop", sep = "")],
                                                        1,
                                                        function(x) -1 * sum(x * log(x, length(target_cell_types))))
  cells_in_neighbourhood_proportion_df$entropy <- ifelse(cells_in_neighbourhood_proportion_df$total > 0 & is.nan(cells_in_neighbourhood_proportion_df$entropy), 
                                                         0,
                                                         cells_in_neighbourhood_proportion_df$entropy)
  
  return(cells_in_neighbourhood_proportion_df)
}
### Start from the grid_prism with the maximum cell proportion.
## Look left, right, forward, back, up and down and see if that grid_prism has at least threshold cell proportion value
## If it does, add it to the answer
## Keep doing this until adjacent grid prisms don't have above threshold, or if you hit a boundary, or it has already been removed
## Return a vector containing all the grid prism numbers which COULD be part of the cluster
calculate_grid_prism_numbers_in_cluster3D <- function(curr_grid_prism_number, 
                                                      grid_prism_cell_proportions, 
                                                      threshold_cell_proportion,
                                                      n_splits,
                                                      answer) {
  
  ## If answer already has curr_grid_prism_number, go back
  if (as.character(curr_grid_prism_number) %in% answer) return(answer)
  
  grid_prism_numbers <- names(grid_prism_cell_proportions)
  
  ## If curr_grid_prism_number has already been removed from grid_prism_numbers, go back
  if (!(as.character(curr_grid_prism_number) %in% grid_prism_numbers)) return(answer)
  
  
  if (grid_prism_cell_proportions[as.character(curr_grid_prism_number)] > threshold_cell_proportion) {
    
    answer <- c(answer, as.character(curr_grid_prism_number))
    
    ### CHECK RIGHT, LEFT, FORWARD, BACKWARD, UP, DOWN
    ## Need to check if going right, left, forward, backward, up or down is possible
    
    # Right
    if (curr_grid_prism_number%%n_splits != 0) {
      answer <- calculate_grid_prism_numbers_in_cluster3D(curr_grid_prism_number + 1,
                                                          grid_prism_cell_proportions,
                                                          threshold_cell_proportion,
                                                          n_splits,
                                                          answer)
    }
    
    # Left
    if (curr_grid_prism_number%%n_splits != 1) {
      answer <- calculate_grid_prism_numbers_in_cluster3D(curr_grid_prism_number - 1,
                                                          grid_prism_cell_proportions,
                                                          threshold_cell_proportion,
                                                          n_splits,
                                                          answer)
    }
    
    # Forward
    if ((curr_grid_prism_number - 1)%%(n_splits^2) < n_splits^2 - n_splits) {
      answer <- calculate_grid_prism_numbers_in_cluster3D(curr_grid_prism_number + n_splits,
                                                          grid_prism_cell_proportions,
                                                          threshold_cell_proportion,
                                                          n_splits,
                                                          answer)
    }
    
    # Backward
    if (curr_grid_prism_number%%(n_splits^2) > n_splits) {
      answer <- calculate_grid_prism_numbers_in_cluster3D(curr_grid_prism_number - n_splits,
                                                          grid_prism_cell_proportions,
                                                          threshold_cell_proportion,
                                                          n_splits,
                                                          answer)
    }
    
    # Up
    if (curr_grid_prism_number <= n_splits^3 - n_splits^2) {
      answer <- calculate_grid_prism_numbers_in_cluster3D(curr_grid_prism_number + n_splits^2,
                                                          grid_prism_cell_proportions,
                                                          threshold_cell_proportion,
                                                          n_splits,
                                                          answer)
    }
    
    # Down
    if (curr_grid_prism_number > n_splits^2) {
      answer <- calculate_grid_prism_numbers_in_cluster3D(curr_grid_prism_number - n_splits^2,
                                                          grid_prism_cell_proportions,
                                                          threshold_cell_proportion,
                                                          n_splits,
                                                          answer)
    }
  }
  
  return(answer)
}#' @title Calculate minimum distances between cell types in 3D spatial data.
#'
#' @description This function calculates the minimum distances between different cell types in a 3D SpatialExperiment object. 
#'    It allows you to specify a subset of cell types to analyse and provides the option to summarise 
#'    the results and plot violin plots of the minimum distances between cell types.
#'
#' @param spe A SpatialExperiment object containing 3D spatial information for the cells.
#' @param cell_types_of_interest A character vector specifying the cell types of interest.
#'   If NULL, all cell types in the `feature_colname` column will be considered.
#' @param feature_colname A string specifying the name of the column in the `colData` slot of the SpatialExperiment
#'    object that contains the cell type information.
#' @param show_summary A logical indicating whether to print a summary of the minimum distances 
#'    for each cell type pair. Defaults to TRUE.
#' @param plot_image A logical indicating whether to plot violin plots of the minimum distances 
#'    between cell type pairs. Defaults to TRUE.
#'
#' @return A data frame containing information about the reference cell, the nearest cell of another type, 
#'    and the distance between them for each cell type pair.
#'
#' @examples
#' minimum_distances <- calculate_minimum_distances_between_cell_types3D(
#'     spe = SPIAT3D::simulated_spe,
#'     cell_types_of_interest = NULL,
#'     feature_colname = "Cell.Type",
#'     show_summary = TRUE,
#'     plot_image = TRUE
#' )
#' 
#' @export


calculate_minimum_distances_between_cell_types3D <- function(spe,
                                                             cell_types_of_interest = NULL,
                                                             feature_colname = "Cell.Type",
                                                             show_summary = TRUE,
                                                             plot_image = TRUE) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (ncol(spe) < 2) {
    stop("There must be at least two cells in spe.")
  }
  if (!(is.null(cell_types_of_interest) || is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(show_summary)) {
    stop("`show_summary` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  if (is.null(spe[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spe")
    spe$Cell.ID <- paste("Cell", seq(ncol(spe)), sep = "_")
  }  
  
  # De-factor feature column in spe object
  spe[[feature_colname]] <- as.character(spe[[feature_colname]])
  
  # Subset spe to only contain the cells of interest
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spe object
    unknown_cell_types <- setdiff(cell_types_of_interest, spe[[feature_colname]])
    if (length(unknown_cell_types) != 0) {
      warning(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                    paste(unknown_cell_types, collapse = ", ")))
    }
    
    spe <- spe[ , spe[[feature_colname]] %in% cell_types_of_interest]
  }
  # If cell_types_of_interest is NULL, use all cells in spe
  else {
    cell_types_of_interest <- unique(spe[[feature_colname]])
  }
  
  # Create a list containing the cell IDs of each cell type
  cell_type_ids <- list()
  for (cell_type in cell_types_of_interest) {
    cell_type_ids[[cell_type]] <- as.character(spe$Cell.ID[spe[[feature_colname]] == cell_type])
  }
  
  # Get spe coords
  spe_coords <- data.frame(spatialCoords(spe))
  
  # Get different possible cell type combinations
  # Each row represents a combination
  # If a row is [1 , 2], then we are comparing cell type 1 and cell type 2
  permu <- gtools::permutations(length(cell_types_of_interest), 2, repeats.allowed = TRUE)
  
  result <- data.frame()
  
  for (i in seq(nrow(permu))) {
    cell_type1 <- cell_types_of_interest[permu[i, 1]]
    cell_type2 <- cell_types_of_interest[permu[i, 2]]
    
    # Don't have one of the cells
    if (sum(spe[[feature_colname]] == cell_type1) == 0 || sum(spe[[feature_colname]] == cell_type2) == 0) {
      result <- rbind(result, data.frame(ref_cell_id = NA, ref_cell_type = cell_type1, nearest_cell_id = NA, nearest_cell_type = cell_type2, distance = NA))
      next
    }
    
    # Get x, y, z coords for all cells of cell_type1 and cell_type2
    cell_type1_coords <- spe_coords[spe[[feature_colname]] == cell_type1, ]
    cell_type2_coords <- spe_coords[spe[[feature_colname]] == cell_type2, ]
    
    # Find all of closest points
    # For each cell of cell_type1, find the closest cell of cell_type2
    if (cell_type1 != cell_type2) {
      nearest_neighbours <- RANN::nn2(data = cell_type2_coords, 
                                      query = cell_type1_coords, 
                                      k = 1)  
    }
    # If we are comparing the same cell_type, and there is only one of this cell type, move on
    else if (nrow(cell_type1_coords) == 1) {
      warning("There is only 1 '", cell_type1, "' cell in your data. It has no nearest neighbour of the same cell type.", sep = "")
      result <- rbind(result, data.frame(ref_cell_id = NA, ref_cell_type = cell_type1, nearest_cell_id = NA, nearest_cell_type = cell_type2, distance = NA))
      next
    }
    # If we are comparing the same cell_type, use the second closest neighbour
    else {
      nearest_neighbours <- RANN::nn2(data = cell_type2_coords, 
                                      query = cell_type1_coords, 
                                      k = 2)
      nearest_neighbours[['nn.idx']] <- nearest_neighbours[['nn.idx']][ , 2]
      nearest_neighbours[['nn.dists']] <- nearest_neighbours[['nn.dists']][ , 2]
    }
    
    # Create the data frame containing the chosen cells and their ids, as well as the nearest cell to them and their ids, and the distance between
    
    df <- data.frame(
      ref_cell_id = cell_type_ids[[cell_type1]],
      ref_cell_type = cell_type1,
      nearest_cell_id = cell_type_ids[[cell_type2]][c(nearest_neighbours$nn.idx)],
      nearest_cell_type = cell_type2,
      distance = nearest_neighbours$nn.dists
    )
    result <- rbind(result, df)
  }
  
  result$pair <- paste(result$ref_cell_type, result$nearest_cell_type,sep = "/")
  
  # Print summary
  if (show_summary) {
    print(summarise_distances_between_cell_types3D(result))  
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_distances_between_cell_types_violin3D(result)
    methods::show(fig)
  }
  
  return(result)
}
calculate_minimum_distances_to_clusters3D <- function(spe, 
                                                      cell_types_inside_cluster, 
                                                      cell_types_outside_cluster, 
                                                      cluster_colname, 
                                                      feature_colname = "Cell.Type", 
                                                      plot_image = T) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (!is.character(cell_types_inside_cluster)) {
    stop("`cell_types_inside_cluster` is not a character vector.")
  }
  if (!is.character(cell_types_outside_cluster)) {
    stop("`cell_types_outside_cluster` is not a character vector.")
  }
  if (!(is.numeric(radius) && length(radius) == 1 && radius > 0)) {
    stop("`radius` is not a positive numeric.")
  }
  if (!is.character(cluster_colname)) {
    stop("`cluster_colname` is not a character. This should be 'alpha_hull_cluster', 'dbscan_cluster', or 'grid_based_cluster', depending on the chosen method.")
  }
  if (is.null(spe[[cluster_colname]])) {
    stop(paste("No column called", cluster_colname, "found in spe object."))
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  ## Add Cell.ID column
  if (is.null(spe[["Cell.ID"]])) {
    warning("Temporarily adding Cell.Id column to your spe")
    spe$Cell.ID <- paste("Cell", seq(ncol(spe)), sep = "_")
  }
  
  ## For each cell type outside clusters, get their set of coords. These exclude cell types in clusters
  spe_coords <- spatialCoords(spe)
  
  # Cells outside cluster have a cluster number of 0 (i.e. they are not in a cluster)
  spe_outside_cluster <- spe[ , spe[[cluster_colname]] == 0]
  
  cell_types_outside_cluster_coords <- list()
  for (cell_type in cell_types_outside_cluster) {
    cell_types_outside_cluster_coords[[cell_type]] <- spatialCoords(spe_outside_cluster)[spe_outside_cluster[[feature_colname]] == cell_type, ]
  }
  
  ## For each cluster, determine the minimum distance of each outside_cell_type  
  result <- vector()
  
  # Get number of clusters
  n_clusters <- max(spe[[cluster_colname]])
  
  for (i in seq(n_clusters)) {
    cluster_coords <- spe_coords[spe[[cluster_colname]] == i & spe[[feature_colname]] %in% cell_types_inside_cluster, ]
    cluster_cell_types <- spe[["Cell.Type"]][spe[[cluster_colname]] == i & spe[[feature_colname]] %in% cell_types_inside_cluster]
    cluster_cell_ids <- spe[["Cell.ID"]][spe[[cluster_colname]] == i & spe[[feature_colname]] %in% cell_types_inside_cluster]
    
    for (outside_cell_type in cell_types_outside_cluster) {
      curr_cell_type_coords <- cell_types_outside_cluster_coords[[outside_cell_type]]
      
      all_closest <- RANN::nn2(data = cluster_coords, 
                               query = curr_cell_type_coords, 
                               k = 1) 
      
      local_dist_mins <- data.frame(
        cluster_number = i,
        outside_cell_id = as.character(spe_outside_cluster$Cell.ID[spe_outside_cluster[["Cell.Type"]] == outside_cell_type]),
        outside_cell_type = outside_cell_type,
        inside_cell_id = cluster_cell_ids[c(all_closest$nn.idx)],
        inside_cell_type = cluster_cell_types[c(all_closest$nn.idx)],
        distance = all_closest$nn.dists
      )
      ## Remove any 0 distance rows
      local_dist_mins <- local_dist_mins[local_dist_mins$distance != 0, ]
      result <- rbind(result, local_dist_mins)
    }
    
    
    ## Plot
    if (plot_image) {
      
      cluster_number_labs <- paste("cluster_", seq(n_clusters), sep = "")
      names(cluster_number_labs) <- seq(n_clusters)
      
      fig <- ggplot(result, aes(x = outside_cell_type, y = distance, fill = outside_cell_type)) + 
        geom_violin() +
        facet_grid(cluster_number~., scales="free_x", labeller = labeller(cluster_number = cluster_number_labs)) +
        theme_bw() +
        theme(axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), legend.position = "none") +
        labs(title="Minimum cell distances to clusters", x = "Cell type", y = "Distance") +
        stat_summary(fun.data = "mean_sdl", fun.args = list(mult= 1), colour = "red")
      
      methods::show(fig)
    }
    
  }
  return(result)
}
calculate_mixing_scores_gradient3D <- function(spe, 
                                               reference_cell_type, 
                                               target_cell_type, 
                                               radii, 
                                               feature_colname = "Cell.Type",
                                               plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = 8))
  colnames(result) <- c("ref_cell_type", 
                        "tar_cell_type", 
                        "n_ref_cells",
                        "n_tar_cells", 
                        "n_ref_tar_interactions",
                        "n_ref_ref_interactions", 
                        "mixing_score", 
                        "normalised_mixing_score")
  
  for (i in seq(length(radii))) {
    mixing_scores <- calculate_mixing_scores3D(spe,
                                               reference_cell_type,
                                               target_cell_type,
                                               radii[i],
                                               feature_colname)
    
    result[i, ] <- mixing_scores
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig1 <- plot_mixing_scores_gradient3D(result, "NMS")
    fig2 <- plot_mixing_scores_gradient3D(result, "MS")
    combined_fig <- plot_grid(fig1, fig2, nrow = 2)
    methods::show(combined_fig)
  }
  
  return(result)
}

calculate_mixing_scores3D <- function(spe, 
                                      reference_cell_types, 
                                      target_cell_types, 
                                      radius, 
                                      feature_colname = "Cell.Type") {
  
  # Define result
  result <- data.frame()
  
  for (reference_cell_type in reference_cell_types) {
    
    for (target_cell_type in target_cell_types) {
      
      # No point getting mixing scores if comparing the same cell type
      if (reference_cell_type == target_cell_type) {
        next
      }
      
      # Get number of reference cells and target cells
      n_ref <- sum(spe[[feature_colname]] == reference_cell_type)
      n_tar <- sum(spe[[feature_colname]] == target_cell_type)
      
      
      # Can't get mixing scores if there are 0 or 1 reference cells
      if (n_ref == 0 || n_ref == 1) {
        result <-  rbind(result, 
                         c(reference_cell_type, 
                           target_cell_type, 
                           n_ref, 
                           n_tar, 
                           0, 
                           0, 
                           NA, 
                           NA))
      }
      
      
      ## Get cells in neighbourhood df
      cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood3D(spe,
                                                                      reference_cell_type,
                                                                      c(reference_cell_type, target_cell_type),
                                                                      radius,
                                                                      feature_colname,
                                                                      FALSE,
                                                                      FALSE)
      
      # Get number of ref-ref interactions
      # Halve it to avoid counting each ref-ref interaction twice
      n_ref_ref_interactions <- 0.5 * sum(cells_in_neighbourhood_df[[reference_cell_type]]) 
      
      # Get number of ref-tar interactions
      n_ref_tar_interactions <- sum(cells_in_neighbourhood_df[[target_cell_type]]) 
      
      
      # Can't get mixing scores if there are no target cells
      if (n_tar == 0) {
        
        result <-  rbind(result, 
                         c(reference_cell_type, 
                           target_cell_type, 
                           n_ref, 
                           0, 
                           0, 
                           n_ref_ref_interactions, 
                           NA, 
                           NA))
      }
      
      # Generic case: We have reference cells and target cells
      else {
        
        if (n_ref_ref_interactions != 0) {
          mixing_score <- n_ref_tar_interactions / n_ref_ref_interactions
          normalised_mixing_score <- 0.5 * mixing_score * n_ref / n_tar
        }
        else {
          mixing_score <- 0
          normalised_mixing_score <- 0
          methods::show(paste("There are no reference to reference interactions for", target_cell_type, "in the specified radius, cannot calculate mixing score"))
        }
        
        result <-  rbind(result, 
                         c(reference_cell_type, 
                           target_cell_type, 
                           n_ref, 
                           n_tar, 
                           n_ref_tar_interactions, 
                           n_ref_ref_interactions, 
                           mixing_score, 
                           normalised_mixing_score))
      }
    }
  }
  
  # Required column names of our output data frame
  colnames(result) <- c("ref_cell_type", 
                        "tar_cell_type", 
                        "n_ref_cells",
                        "n_tar_cells", 
                        "n_ref_tar_interactions",
                        "n_ref_ref_interactions", 
                        "mixing_score", 
                        "normalised_mixing_score")
  
  # Turn numeric data into numeric type
  result[ , 3:8] <- apply(result[ , 3:8], 2, as.numeric)
  
  return(result)
}
calculate_pairwise_distances_between_cell_types3D <- function(spe,
                                                              cell_types_of_interest = NULL,
                                                              feature_colname = "Cell.Type",
                                                              show_summary = TRUE,
                                                              plot_image = TRUE) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (ncol(spe) < 2) {
    stop("There must be at least two cells in spe.")
  }
  if (!(is.null(cell_types_of_interest) && is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(show_summary)) {
    stop("`show_summary` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  if (is.null(spe[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spe")
    spe$Cell.ID <- paste("Cell", seq(ncol(spe)), sep = "_")
  }
  
  # De-factor feature column in spe object
  spe[[feature_colname]] <- as.character(spe[[feature_colname]])
  
  # Subset spe to only contain the cells of interest
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spe object
    unknown_cell_types <- setdiff(cell_types_of_interest, spe[[feature_colname]])
    if (length(unknown_cell_types) != 0) {
      warning(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                    paste(unknown_cell_types, collapse = ", ")))
    }
    
    spe <- spe[ , spe[[feature_colname]] %in% cell_types_of_interest]
  }
  # If cell_types_of_interest is NULL, use all cells in spe
  else {
    cell_types_of_interest <- unique(spe[[feature_colname]])
  }
  
  # Create a list containing the cell IDs of each cell type
  cell_type_ids <- list()
  for (cell_type in cell_types_of_interest) {
    cell_type_ids[[cell_type]] <- as.character(spe$Cell.ID[spe[[feature_colname]] == cell_type])
  }
  
  # Calculate cell to cell distances
  distance_matrix <- -1 * apcluster::negDistMat(spatialCoords(spe))
  rownames(distance_matrix) <- spe$Cell.ID
  colnames(distance_matrix) <- spe$Cell.ID
  
  result <- data.frame()
  
  for (i in seq(length(cell_types_of_interest))) {
    
    for (j in i:length(cell_types_of_interest)) {
      
      # Get current cell types and cell ids
      cell_type1 <- names(cell_type_ids)[i]
      cell_type2 <- names(cell_type_ids)[j]
      
      cell_type1_ids <- cell_type_ids[[cell_type1]]
      cell_type2_ids <- cell_type_ids[[cell_type2]]
      
      ## Don't have a cell type, or the same cell type with only one cell
      if (length(cell_type1_ids) == 0 || length(cell_type2_ids) == 0) {
        result <- rbind(result, data.frame(Var1 = NA, Var2 = NA, value = NA, cell_type1 = cell_type1, cell_type2 = cell_type2, pair = paste(cell_type1, cell_type2, sep="/")))
        next
      }
      
      ## Same cell type only one cell
      if (cell_type1 == cell_type2 && length(cell_type1_ids) == 1) {
        warning("There is only 1 '", cell_type1, "' cell in your data. It has no pair of the same cell type.", sep = "")
        result <- rbind(result, data.frame(Var1 = NA, Var2 = NA, value = NA, cell_type1 = cell_type1, cell_type2 = cell_type2, pair = paste(cell_type1, cell_type2, sep="/")))
        next
      }
      
      # Subset distance_matrix for current cell types
      distance_matrix_subset <- distance_matrix[rownames(distance_matrix) %in% cell_type1_ids, 
                                                colnames(distance_matrix) %in% cell_type2_ids]
      
      ## Different cell types, each only has one cell
      if (length(cell_type1_ids) == 1 && length(cell_type2_ids) == 1) {
        distance_matrix_subset <- as.matrix(distance_matrix_subset)
        rownames(distance_matrix_subset) <- cell_type1_ids
        colnames(distance_matrix_subset) <- cell_type2_ids
      }    
      ## Different cell types, only one cell of cell_type1
      else if (length(cell_type1_ids) == 1) {
        distance_matrix_subset <- as.matrix(distance_matrix_subset)
        colnames(distance_matrix_subset) <- cell_type1_ids
      }
      ## Different cell types, only one cell of cell_type2
      else if (length(cell_type2_ids) == 1) {
        distance_matrix_subset <- as.matrix(distance_matrix_subset)
        colnames(distance_matrix_subset) <- cell_type2_ids
      }
      ## Same cell type, only need part of the matrix (make irrelevant part of matrix equal to NA)
      if (cell_type1 == cell_type2) distance_matrix_subset[upper.tri(distance_matrix_subset, diag = TRUE)] <- NA
      
      # Convert distance_matrix_subset to a data frame
      df <- reshape2::melt(distance_matrix_subset, na.rm = TRUE)
      df$cell_type1 <- cell_type1
      df$cell_type2 <- cell_type2
      df$pair <- paste(cell_type1, cell_type2, sep="/")
      
      result <- rbind(result, df)
    }
  }
  
  # Rearrange columns 
  colnames(result)[c(1, 2, 3)] <- c("cell_type1_id", "cell_type2_id", "distance")
  result <- result[ , c("cell_type1_id", "cell_type1", "cell_type2_id", "cell_type2", "distance", "pair")]
  
  # Print summary
  if (show_summary) {
    print(summarise_distances_between_cell_types3D(result))  
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_distances_between_cell_types_violin3D(result)
    methods::show(fig)
  }
  
  return(result)
}
calculate_prevalence_gradient_AUC3D <- function(prevalence_gradient_df) {
  
  return(sum(prevalence_gradient_df$prevalence) * 0.01)
}
calculate_prevalence_gradient3D <- function(grid_metrics,
                                            metric_colname,
                                            show_AUC = T,
                                            plot_image = T) {
  
  ## Check input parameters
  if (!(is.character(metric_colname))) {
    stop("`metric_colname` is not a character. This should be 'proportion' or 'entropy', depending on the chosen method.")
  }
  if (is.null(grid_metrics[[metric_colname]])) {
    stop("`metric_colname` is not a column in `grid_metrics`.")
  }
  if (!is.logical(show_AUC)) {
    stop("`show_AUC` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Thresholds range from 0 to 1
  thresholds <- seq(0.01, 1, 0.01)
  
  # Define result
  result <- data.frame(threshold = thresholds)
  
  # Get prevalences for each threshold
  result$prevalence <- sapply(thresholds, function(threshold) { 
    calculate_prevalence3D(grid_metrics, metric_colname, threshold) 
  })
  
  # Show AUC of prevalence gradient graph
  if (show_AUC) {
    print(paste("AUC:", round(calculate_prevalence_gradient_AUC3D(result), 2)))
  }
  
  # Plot
  if (plot_image) {
    fig <- ggplot(result, aes(threshold, prevalence)) +
      geom_line() +
      theme_bw() +
      labs(x = "Threshold",
           y = "Prevalence",
           title = paste("Prevalence vs Threshold (", metric_colname, ")", sep = "")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      ylim(0, 100)
    methods::show(fig)
  }
  
  return(result)
}
calculate_prevalence3D <- function(grid_metrics,
                                   metric_colname,
                                   threshold,
                                   above = TRUE) {
  
  ## Check input parameters
  if (!(is.character(metric_colname))) {
    stop("`metric_colname` is not a character. This should be 'proportion' or 'entropy', depending on the chosen method.")
  }
  if (is.null(grid_metrics[[metric_colname]])) {
    stop("`metric_colname` is not a column in `grid_metrics`.")
  }
  if (!(is.numeric(threshold) && length(threshold) == 1 && threshold >= 0 && threshold <= 1)) {
    stop("`threshold` is not a numeric between 0 and 1.")
  }
  if (!is.logical(above)) {
    stop("`above` is not a logical (TRUE or FALSE).")
  }
  
  
  ## Exclude rows with NA values
  grid_metrics <- grid_metrics[!is.na(grid_metrics[[metric_colname]]), ]
  
  if (above) {
    prevalence <- sum(grid_metrics[[metric_colname]] >= threshold) / nrow(grid_metrics) * 100
  }
  else {
    prevalence <- sum(grid_metrics[[metric_colname]] < threshold) / nrow(grid_metrics) * 100    
  }
  
  return(prevalence)
}
calculate_spatial_autocorrelation3D <- function(grid_metrics,
                                                metric_colname,
                                                weight_method = 0.1) {
  
  ## Check input parameters
  if (!(is.character(metric_colname))) {
    stop("`metric_colname` is not a character. This should be 'proportion' or 'entropy', depending on the chosen method.")
  }
  if (is.null(grid_metrics[[metric_colname]])) {
    stop("`metric_colname` is not a column in `grid_metrics`.")
  }
  if (!((is.numeric(weight_method) && length(weight_method) == 1 && weight_method > 0 && weight_method < 1) ||
        (is.character(weight_method) && weight_method %in% c("IDW", "rook", "queen")))) {
    stop("`weight_method` is not a numeric between 0 and 1 or either 'IDW', 'rook' or 'queen'.")
  }
  
  ## Get number of grid prisms
  n_grid_prisms <- nrow(grid_metrics)
  
  ## Get splitting number (should be the cube root of n_grid_prisms)
  n_splits <- (n_grid_prisms)^(1/3)
  
  ## Find the coordinates of each grid prism
  x <- ((seq(n_grid_prisms) - 1) %% n_splits)
  y <- (floor(((seq(n_grid_prisms) - 1) %% (n_splits)^2) / n_splits))
  z <- (floor((seq(n_grid_prisms) - 1) / (n_splits^2)))
  grid_prism_coords <- data.frame(x = x, y = y, z = z)
  
  ## Subset for non NA rows
  grid_prism_coords <- grid_prism_coords[!is.na(grid_metrics[[metric_colname]]), ]
  grid_metrics <- grid_metrics[!is.na(grid_metrics[[metric_colname]]), ]
  
  weight_matrix <- -1 * apcluster::negDistMat(grid_prism_coords)
  ## Use the inverse distance between two points as the weight (IDW is 'inverse distance weighting')
  if (weight_method == "IDW") {
    weight_matrix <- 1 / weight_matrix
  }
  ## Use rook method: adjacent points get a weight of 1, otherwise, weight of 0
  ## Adjacent points are within 1 unit apart. e.g. (0, 0, 0) vs (0, 0, 1)
  else if (weight_method == "rook") {
    weight_matrix <- ifelse(weight_matrix > 1, 0, 1)  
  }
  ## Use queen method: adjacent points get a weight of 1, otherwise, weight of 0
  ## Adjacent points are within sqrt(3) unit apart. e.g. (0, 0, 0) vs (0, 0, 1)
  else if (weight_method == "queen") {
    weight_matrix <- ifelse(weight_matrix > sqrt(3), 0, 1)  
  }
  ## If a number (x) between 0 and 1 is supplied, set a threshold to be x quantile value of c(weight_matrix)
  ## Grid prisms within this specified threshold have a weight of 1, otherwise, weight of 0
  else if (as.numeric(weight_method) && 0 < weight_method && weight_method < 1) {
    threshold <- quantile(c(weight_matrix), weight_method)
    weight_matrix <- ifelse(weight_matrix > threshold, 0, 1)
  }
  
  ## Points along the diagonal are comparing the same point so its weight is zero
  diag(weight_matrix) <- 0
  
  n <- nrow(grid_metrics)
  
  # Center the data
  data_centered <- grid_metrics[[metric_colname]] - mean(grid_metrics[[metric_colname]])
  
  # Calculate numerator using matrix multiplication
  numerator <- sum(data_centered * (weight_matrix %*% data_centered))
  
  # Calculate denominator
  denominator <- sum(data_centered^2) * sum(weight_matrix)
  
  # Moran's I
  I <- (n * numerator) / denominator
  
  return(I)
}
calculate_volume_of_clusters3D <- function(spe, cluster_colname) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (!is.character(cluster_colname)) {
    stop("`cluster_colname` is not a character. This should be 'alpha_hull_cluster', 'dbscan_cluster', or 'grid_based_cluster', depending on the chosen method.")
  }
  if (is.null(spe[[cluster_colname]])) {
    stop(paste("No column called", cluster_colname, "found in spe object."))
  }
  
  # Get number of clusters
  n_clusters <- max(spe[[cluster_colname]])
  
  ### 1. Estimate volume of each cluster by density of the window. ------------
  
  ## For each cluster, determine the number of cells in each cluster of each cluster
  result <- data.frame(matrix(nrow = n_clusters, ncol = 2))
  colnames(result) <- c("cluster_number", "n_cells")
  
  for (i in seq(n_clusters)) {
    result[i, "n_cells"] <- sum(spe[[cluster_colname]] == i)
  }
  result$cluster_number <- as.character(seq(n_clusters))
  
  ## Assume window is a rectangular prism
  spe_coords <- data.frame(spatialCoords(spe))
  
  length <- round(max(spe_coords$Cell.X.Position) - min(spe_coords$Cell.X.Position))
  width  <- round(max(spe_coords$Cell.Y.Position) - min(spe_coords$Cell.Y.Position))
  height <- round(max(spe_coords$Cell.Z.Position) - min(spe_coords$Cell.Z.Position))
  
  window_volume <- length * width * height
  
  result$volume_by_density <- (result$n_cells / ncol(spe)) * window_volume
  
  
  ### 2. If cluster_colname == "alpha_hull_cluster", use the volume method found in the alphashape3d package
  if (cluster_colname == "alpha_hull_cluster") {
    result$volume_by_alpha_hull <- volume_ashape3d(spe@metadata$alpha_hull$ashape3d_object, byComponents = T)
  }
  
  
  ### 3. If cluster_colname == "grid_based_cluster", sum the volume of each grid prism to get volume of each cluster
  if (cluster_colname == "grid_based_cluster") {
    result$volume_by_grid <- 0
    i <- 1
    for (grid_cluster in spe@metadata$grid_prisms) {
      result[i, "volume_by_grid"] <- sum(grid_cluster$l * grid_cluster$w * grid_cluster$h)
      i <- i + 1
    }
  }
  
  return(result)
}
library(dbscan)

dbscan_clustering3D <- function(spe,
                                cell_types_of_interest,
                                radius,
                                minimum_cells_in_radius,
                                minimum_cells_in_cluster,
                                feature_colname = "Cell.Type",
                                plot_image = T) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  ## Check cell types of interst are found in the spe object
  unknown_cell_types <- setdiff(cell_types_of_interest, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  if (!(is.numeric(radius) && length(radius) == 1 && radius > 0)) {
    stop("`radius` is not a positive numeric.")
  }
  if (!(is.integer(minimum_cells_in_radius) && length(minimum_cells_in_radius) == 1 || (is.numeric(minimum_cells_in_radius) && length(minimum_cells_in_radius) == 1 && minimum_cells_in_radius > 0 && minimum_cells_in_radius%%1 == 0))) {
    stop("`minimum_cells_in_radius` is not a positive integer.")
  }
  if (!(is.integer(minimum_cells_in_cluster) && length(minimum_cells_in_cluster) == 1 || (is.numeric(minimum_cells_in_cluster) && length(minimum_cells_in_cluster) == 1 && minimum_cells_in_cluster > 0 && minimum_cells_in_cluster%%1 == 0))) {
    stop("`minimum_cells_in_cluster` is not a positive integer.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  spe_subset <- spe[ , spe[[feature_colname]] %in% cell_types_of_interest]
  spe_subset_coords <- spatialCoords(spe_subset)
  
  db <- dbscan::dbscan(spe_subset_coords, eps = radius, minPts = minimum_cells_in_radius, borderPoints = F)
  n_clusters <- max(db$cluster)
  if (n_clusters == 0) {
    stop("No clusters identified. Consider increasing `radius` and/or decreasing `minimum_cells_in_radius`.")
  }
  
  ## Cell types of interest have a 'cluster' value of 0 if they are noise, 1 if they belong to cluster 1, ...
  ## Check if number of cells in cluster 1, cluster 2, ... is larger than minimum_cells_in_cluster, if they don't, these cells are also assigned a value of 0.
  for (i in seq_len(n_clusters)) {
    if (sum(db$cluster == i) < minimum_cells_in_cluster) {
      db$cluster[db$cluster == i] <- 0
      
      ## Re-number the clusters
      db$cluster[db$cluster > i] <- db$cluster[db$cluster > i] - 1
    }
  }
  n_clusters <- max(db$cluster)
  if (n_clusters == 0) {
    stop("All clusters identified do not meet the `minimum_cells_in_cluster` threshold. Consider lowering the `minimum_cells_in_cluster` parameter.")
  }
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(spe), colData(spe))
  
  df_cell_types_of_interest <- df[df[[feature_colname]] %in% cell_types_of_interest, ]
  df_other_cell_types <- df[!(df[[feature_colname]] %in% cell_types_of_interest), ]
  
  df_cell_types_of_interest$dbscan_cluster <- db$cluster
  df_other_cell_types$dbscan_cluster <- 0
  
  ## Convert data frame to spe object
  df <- rbind(df_cell_types_of_interest, df_other_cell_types)
  
  spe <- SpatialExperiment(
    assay = matrix(data = NA, nrow = nrow(df), ncol = nrow(df)),
    colData = df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = spe@metadata)
  
  ## Plot
  if (plot_image) {
    df$dbscan_cluster <- ifelse(df$dbscan_cluster == 0, "non_cluster", paste("cluster_", df$dbscan_cluster, sep = ""))
    
    fig <- plot_ly(df,
                   type = "scatter3d",
                   mode = 'markers',
                   x = ~Cell.X.Position,
                   y = ~Cell.Y.Position,
                   z = ~Cell.Z.Position,
                   color = ~dbscan_cluster,
                   colors = rainbow(length(unique(df$dbscan_cluster))),
                   marker = list(size = 2)) %>% 
      layout(scene = list(xaxis = list(title = 'x'),
                          yaxis = list(title = 'y'),
                          zaxis = list(title = 'z')))
    
    methods::show(fig)
  }
  
  return(spe)
}
get_spe_grid_metrics3D <- function(spe, 
                                   n_splits, 
                                   feature_colname = "Cell.Type") {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  
  spe_coords <- spatialCoords(spe)
  
  ## Get dimensions of the window
  min_x <- min(spe_coords[ , "Cell.X.Position"])
  min_y <- min(spe_coords[ , "Cell.Y.Position"])
  min_z <- min(spe_coords[ , "Cell.Z.Position"])
  
  max_x <- max(spe_coords[ , "Cell.X.Position"])
  max_y <- max(spe_coords[ , "Cell.Y.Position"])
  max_z <- max(spe_coords[ , "Cell.Z.Position"])
  
  length <- round(max_x - min_x)
  width  <- round(max_y - min_y)
  height <- round(max_z - min_z)
  
  ## Get distance of row, col and lay
  d_row <- length / n_splits
  d_col <- width / n_splits
  d_lay <- height / n_splits
  
  # Shift spe_coords so they begin at the origin
  spe_coords[, "Cell.X.Position"] <- spe_coords[, "Cell.X.Position"] - min_x
  spe_coords[, "Cell.Y.Position"] <- spe_coords[, "Cell.Y.Position"] - min_y
  spe_coords[, "Cell.Z.Position"] <- spe_coords[, "Cell.Z.Position"] - min_z
  
  ## Figure out which 'grid prism number' each cell is inside
  spe$grid_prism_num <- floor(spe_coords[ , "Cell.X.Position"] / d_row) +
    floor(spe_coords[ , "Cell.Y.Position"] / d_col) * n_splits + 
    floor(spe_coords[ , "Cell.Z.Position"] / d_lay) * n_splits^2 + 1
  
  ## Determine the cell types found in each grid prism
  n_grid_prisms <- n_splits^3
  grid_prism_cell_matrix <- as.data.frame.matrix(table(spe[[feature_colname]], factor(spe$grid_prism_num, levels = seq(n_grid_prisms))))
  grid_prism_cell_matrix <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                       t(grid_prism_cell_matrix))
  
  ## Determine centre coordinates of each grid prism
  grid_prism_coordinates <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                       x_coord = ((seq(n_grid_prisms) - 1) %% n_splits + 0.5) * d_row + round(min_x),
                                       y_coord = (floor(((seq(n_grid_prisms) - 1) %% (n_splits)^2) / n_splits) + 0.5) * d_col + round(min_y),
                                       z_coord = (floor((seq(n_grid_prisms) - 1) / (n_splits^2)) + 0.5) * d_lay + round(min_z))
  
  spe@metadata[["grid_metrics"]] <- list("grid_prism_cell_matrix" = grid_prism_cell_matrix,
                                         "grid_prism_coordinates" = grid_prism_coordinates)
  
  return(spe)
}
grid_based_cluster_recursion3D <- function(df,  # Using a df is much faster than using a spe
                                           cell_types_of_interest,
                                           threshold_cell_proportion,
                                           x, y, z, l, w, h,
                                           feature_colname,
                                           answer) {
  
  # Look at cells only in the current grid prism
  df <- df[df$Cell.X.Position >= x &
             df$Cell.X.Position < (x + l) &
             df$Cell.Y.Position >= y &
             df$Cell.Y.Position < (y + w) &
             df$Cell.Z.Position >= z &
             df$Cell.Z.Position < (z + h), ]
  
  # Get cell types from spe grid prism
  cell_types <- df[[feature_colname]]
  
  # Number of cells in prism is getting too small
  if (length(cell_types) <= 2) return(data.frame())
  
  # Get total cell proportion for chosen cell_types_of_interest
  cell_proportion <- mean(cell_types %in% cell_types_of_interest)
  
  # Keep grid prism if cell proportion is above the threshold cell proportion
  if (cell_proportion >= threshold_cell_proportion) {
    return(data.frame(x, y, z, l, w, h))
  }
  
  # some cell_types_of_interest still in the grid prism, check sub-grid prisms (8 to check)
  else if (cell_proportion > 0) {
    # (0, 0, 0)
    answer <- rbind(answer, grid_based_cluster_recursion3D(df,
                                                           cell_types_of_interest,
                                                           threshold_cell_proportion,
                                                           x, y, z, l/2, w/2, h/2,
                                                           feature_colname,
                                                           data.frame()))
    
    # (0.5, 0, 0)
    answer <- rbind(answer, grid_based_cluster_recursion3D(df,
                                                           cell_types_of_interest,
                                                           threshold_cell_proportion,
                                                           x + l/2, y, z, l/2, w/2, h/2,
                                                           feature_colname,
                                                           data.frame()))
    
    # (0, 0.5, 0)
    answer <- rbind(answer, grid_based_cluster_recursion3D(df,
                                                           cell_types_of_interest,
                                                           threshold_cell_proportion,
                                                           x, y + w/2, z, l/2, w/2, h/2,
                                                           feature_colname,
                                                           data.frame()))
    # (0.5, 0.5, 0)
    answer <- rbind(answer, grid_based_cluster_recursion3D(df,
                                                           cell_types_of_interest,
                                                           threshold_cell_proportion,
                                                           x + l/2, y + w/2, z, l/2, w/2, h/2,
                                                           feature_colname,
                                                           data.frame()))
    
    # (0, 0, 0.5)
    answer <- rbind(answer, grid_based_cluster_recursion3D(df,
                                                           cell_types_of_interest,
                                                           threshold_cell_proportion,
                                                           x, y, z + h/2, l/2, w/2, h/2,
                                                           feature_colname,
                                                           data.frame()))
    
    # (0.5, 0, 0.5)
    answer <- rbind(answer, grid_based_cluster_recursion3D(df,
                                                           cell_types_of_interest,
                                                           threshold_cell_proportion,
                                                           x + l/2, y, z + h/2, l/2, w/2, h/2,
                                                           feature_colname,
                                                           data.frame()))
    
    # (0, 0.5, 0.5)
    answer <- rbind(answer, grid_based_cluster_recursion3D(df,
                                                           cell_types_of_interest,
                                                           threshold_cell_proportion,
                                                           x, y + w/2, z + h/2, l/2, w/2, h/2,
                                                           feature_colname,
                                                           data.frame()))
    # (0.5, 0.5, 0.5)
    answer <- rbind(answer, grid_based_cluster_recursion3D(df,
                                                           cell_types_of_interest,
                                                           threshold_cell_proportion,
                                                           x + l/2, y + w/2, z + h/2, l/2, w/2, h/2,
                                                           feature_colname,
                                                           data.frame()))
    
    return(answer)
  }
  
  # cell proportion is zero
  else {
    return(data.frame())
  }
}
grid_based_clustering3D <- function(spe,
                                    cell_types_of_interest,
                                    n_splits,
                                    minimum_cells_in_cluster,
                                    feature_colname = "Cell.Type",
                                    plot_image = TRUE) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  ## Check cell types of interst are found in the spe object
  unknown_cell_types <- setdiff(cell_types_of_interest, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  if (!(is.integer(minimum_cells_in_cluster) && length(minimum_cells_in_cluster) == 1 || (is.numeric(minimum_cells_in_cluster) && length(minimum_cells_in_cluster) == 1 && minimum_cells_in_cluster > 0 && minimum_cells_in_cluster%%1 == 0))) {
    stop("`minimum_cells_in_cluster` is not a positive integer.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Add grid metrics to spe
  spe <- get_spe_grid_metrics3D(spe, n_splits, feature_colname)
  
  # Get grid_prism_cell_matrix from spe
  grid_prism_cell_matrix <- spe@metadata$grid_metrics$grid_prism_cell_matrix
  
  ## Calculate proportions for each grid prism
  if (length(cell_types_of_interest) == 1) {
    grid_prism_cell_proportions <- grid_prism_cell_matrix[ , cell_types_of_interest]
  }
  else {
    grid_prism_cell_proportions <- rowSums(grid_prism_cell_matrix[ , cell_types_of_interest])
  }
  grid_prism_cell_proportions <- grid_prism_cell_proportions / rowSums(grid_prism_cell_matrix[ , unique(spe[[feature_colname]])])
  n_grid_prisms <- n_splits^3
  names(grid_prism_cell_proportions) <- seq(n_grid_prisms)
  
  
  ## Create template for final result
  result <- list()
  n_clusters <- 1
  
  ## Get dimensions of the window
  spe_coords <- data.frame(spatialCoords(spe))
  
  min_x <- min(spe_coords$Cell.X.Position)
  min_y <- min(spe_coords$Cell.Y.Position)
  min_z <- min(spe_coords$Cell.Z.Position)
  
  max_x <- max(spe_coords$Cell.X.Position)
  max_y <- max(spe_coords$Cell.Y.Position)
  max_z <- max(spe_coords$Cell.Z.Position)
  
  length <- round(max_x - min_x)
  width  <- round(max_y - min_y)
  height <- round(max_z - min_z)
  
  ## Get distance of row, col and lay
  d_row <- length / n_splits
  d_col <- width / n_splits
  d_lay <- height / n_splits
  
  
  ### CLUSTER DETECTION RECURSIVE ALGORITHM LOOP ###
  
  # First, remove all 0s and NANs from grid_prism_cell_proportions
  grid_prism_cell_proportions <- grid_prism_cell_proportions[grid_prism_cell_proportions != 0 & !is.nan(grid_prism_cell_proportions)]
  
  while (length(grid_prism_cell_proportions) != 0) {
    # Get the maximum cell proportion and its corresponding grid prism number
    maximum_cell_proportion <- max(grid_prism_cell_proportions)
    maximum_cell_proportion_prism_number <- as.numeric(names(which.max(grid_prism_cell_proportions)))
    
    # Break out the loop if maximum cell proportion is less than 0.5
    if (maximum_cell_proportion < 0.5) break 
    
    # Else, find all the grid prisms adjacent to the maximum cell proportion grid prism. 
    # These are potentially apart of the cluster
    # Adjacent grid prisms must have cell proportion > 0.25 * max cell proportion
    grid_prisms_in_cluster <- calculate_grid_prism_numbers_in_cluster3D(maximum_cell_proportion_prism_number,
                                                                        grid_prism_cell_proportions,
                                                                        0.25 * maximum_cell_proportion,
                                                                        n_splits,
                                                                        c())
    
    # Perform the recursive algorithm on each grid prism potentially apart of the cluster to get a more precise shape of each cluster
    # Create data frame with spatial coords and cell types as columns. Use this as input
    result[[n_clusters]] <- data.frame()
    df <- spe_coords
    df[[feature_colname]] <- spe[[feature_colname]] 
    for (grid_prism in as.numeric(grid_prisms_in_cluster)) {
      result[[n_clusters]] <- rbind(result[[n_clusters]],
                                    grid_based_cluster_recursion3D(df,
                                                                   cell_types_of_interest,
                                                                   0.75 * maximum_cell_proportion,
                                                                   ((grid_prism - 1) %% n_splits) * d_row + round(min_x),
                                                                   (floor(((grid_prism - 1) %% n_splits^2) / n_splits)) * d_col + round(min_y),
                                                                   (floor((grid_prism - 1) / n_splits^2)) * d_lay + round(min_z),
                                                                   d_row, d_col, d_lay,
                                                                   feature_colname,
                                                                   data.frame()))
      
      
    }
    colnames(result[[n_clusters]]) <- c("x", "y", "z", "l", "w", "h")
    n_clusters <- n_clusters + 1
    
    # Remove grid prisms which have just been examined
    grid_prism_cell_proportions <- grid_prism_cell_proportions[setdiff((names(grid_prism_cell_proportions)), 
                                                                       grid_prisms_in_cluster)]
    
  }
  # Name each grid_based cluster
  names(result) <- paste("cluster", seq_len(length(result)), sep = "_")
  
  ## Add grid_based_cluster column to spe, indicating which cluster each cell belongs to
  spe$grid_based_cluster <- 0
  cluster_number <- 1
  
  for (i in seq_len(length(result))) {
    cluster_info <- result[[paste("cluster", i, sep = "_")]]
    for (j in seq(nrow(cluster_info))) {
      x <- cluster_info$x[j]
      y <- cluster_info$y[j]
      z <- cluster_info$z[j]
      l <- cluster_info$l[j]
      w <- cluster_info$w[j]
      h <- cluster_info$h[j]
      
      spe$grid_based_cluster <- ifelse(spe_coords$Cell.X.Position >= x &
                                         spe_coords$Cell.X.Position < (x + l) &
                                         spe_coords$Cell.Y.Position >= y &
                                         spe_coords$Cell.Y.Position < (y + w) &
                                         spe_coords$Cell.Z.Position >= z &
                                         spe_coords$Cell.Z.Position < (z + h) &
                                         spe[[feature_colname]] %in% cell_types_of_interest, 
                                       cluster_number, 
                                       spe$grid_based_cluster)
      
    }
    # Check if current cluster surpasses the minimum_cells_in_cluster threshold
    if (sum(spe$grid_based_cluster == cluster_number) < minimum_cells_in_cluster) {
      spe$grid_based_cluster[spe$grid_based_cluster == cluster_number] <- 0
      result[[paste("cluster", i, sep = "_")]] <- NULL
      
    }
    else {
      cluster_number <- cluster_number + 1 
    }
  }
  
  n_clusters <- max(spe$grid_based_cluster)
  if (n_clusters == 0) {
    stop("All clusters identified do not meet the `minimum_cells_in_cluster` threshold. Consider lowering the `minimum_cells_in_cluster` parameter.")
  }
  
  # re-name each grid_based cluster
  names(result) <- paste("cluster", seq_len(length(result)), sep = "_")
  
  # Add grid_clustering result to spe metadata
  spe@metadata[["grid_prisms"]] <- result
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_based_clusters3D(spe, feature_colname = feature_colname)
    methods::show(fig)
  }
  
  return(spe)
}
plot_alpha_hull_clusters3D <- function(spe_with_alpha_hull, 
                                       plot_cell_types = NULL,
                                       plot_colours = NULL,
                                       feature_colname = "Cell.Type") {
  
  # Check input parameters
  if (class(spe_with_alpha_hull) != "SpatialExperiment") {
    stop("`spe_with_alpha_hull` is not a SpatialExperiment object.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe_with_alpha_hull[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  
  ## If no cell types chosen, use all cell types found in data frame
  if (is.null(plot_cell_types)) plot_cell_types <- unique(spe_with_alpha_hull[[feature_colname]])
  
  ## If cell types have been chosen, check they are found in the spe object
  unknown_cell_types <- setdiff(plot_cell_types, spe_with_alpha_hull[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following plot_cell_types are not found in the spe object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  
  ## If no colours inputted, use rainbow palette
  if (is.null(plot_colours)) {
    plot_colours <- rainbow(length(plot_cell_types))
  }
  
  ## User inputs mismatching cell types and colours
  if (length(plot_cell_types) != length(plot_colours)) {
    stop("Length of plot_cell_types is not equal to length of plot_colours")
  }
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(spe_with_alpha_hull), "Cell.Type" = spe_with_alpha_hull[[feature_colname]])
  
  ## Factor for feature column
  df[["Cell.Type"]] <- factor(df[, "Cell.Type"],
                              levels = plot_cell_types)
  
  ## Add points to fig
  fig <- plot_ly() %>%
    add_trace(
      data = df,
      type = "scatter3d",
      mode = 'markers',
      x = ~Cell.X.Position,
      y = ~Cell.Y.Position,
      z = ~Cell.Z.Position,
      marker = list(size = 2),
      color = ~Cell.Type,
      colors = plot_colours
    ) %>% 
    layout(scene = list(xaxis = list(title = 'x'),
                        yaxis = list(title = 'y'),
                        zaxis = list(title = 'z')))
  
  
  ## Get alpha hull numbers (ignoring 0)
  alpha_hull_clusters <- spe_with_alpha_hull$alpha_hull_cluster[spe_with_alpha_hull$alpha_hull_cluster != 0]
  
  # Get number of alpha hulls
  n_alpha_hulls <- length(unique(alpha_hull_clusters))
  
  vertices <- spe_with_alpha_hull@metadata$alpha_hull$vertices
  faces <- data.frame(spe_with_alpha_hull@metadata$alpha_hull$faces)
  alpha_hull_colours <- rainbow(n_alpha_hulls)
  
  ## Add alpha hulls to fig, one by one  
  for (i in seq(n_alpha_hulls)) {
    faces_temp <- faces[faces[ , 1] %in% which(alpha_hull_clusters == i) , ]
    
    ## Ignore the weird cases where some cells represent clusters, but no faces are associated with them??
    if (nrow(faces_temp) == 0) next
    
    opacity_level <- 0.20
    
    fig <- fig %>%
      add_trace(
        type = 'mesh3d',
        x = vertices[, 1], 
        y = vertices[, 2], 
        z = vertices[, 3],
        i = faces_temp[, 1] - 1, 
        j = faces_temp[, 2] - 1, 
        k = faces_temp[, 3] - 1,
        opacity = opacity_level,
        facecolor = rep(alpha_hull_colours[i], nrow(faces_temp))
      )
  }
  
  return(fig)
}
plot_cells_in_neighbourhood_gradient3D <- function(cells_in_neighbourhood_gradient_df, reference_cell_type = NULL) {
  
  plot_result <- reshape2::melt(cells_in_neighbourhood_gradient_df, "radius")
  
  fig <- ggplot(plot_result, aes(radius, value, color = variable)) + 
    geom_line() + 
    labs(title = "Average cells in neighbourhood gradient", x = "Radius", y = "Average cells in neighbourhood") + 
    scale_color_discrete(name = "Cell type") +
    theme_bw()
  
  if (!is.null(reference_cell_type)) {
    fig <- fig + labs(subtitle = paste("Reference: ", reference_cell_type, ", Target: ", paste(colnames(cells_in_neighbourhood_gradient_df)[seq(ncol(cells_in_neighbourhood_gradient_df) - 1)], collapse = ", "), sep = ""))
  }
  
  return(fig)
}
plot_cells_in_neighbourhood_proportions_gradient3D <- function(cells_in_neighbourhood_proportions_gradient_df, reference_cell_type = NULL) {
  
  plot_result <- reshape2::melt(cells_in_neighbourhood_proportions_gradient_df, id.vars = c("radius"))
  fig <- ggplot(plot_result, aes(radius, value, color = variable)) +
    geom_point() +
    geom_line() +
    labs(title = "Average cells in neighbourhood proportions gradient", x = "Radius", y = "Cell proportion", color = "Cell type") +
    theme_bw() +
    ylim(0, 1)
  
  if (!is.null(reference_cell_type)) {
    fig <- fig + labs(subtitle = paste("Reference: ", reference_cell_type, ", Target: ", paste(colnames(cells_in_neighbourhood_proportions_gradient_df)[seq(ncol(cells_in_neighbourhood_proportions_gradient_df) - 1)], collapse = ", "), sep = ""))
  }
  
  return(fig)
}
## For scales parameter, use "free_x" or "free". "free_y" looks silly
plot_cells_in_neighbourhood_violin3D <- function(cells_in_neighbourhood_df, reference_cell_type, scales = "free_x") {
  
  ## Target cell types will be all the columns except the first column
  target_cell_types <- colnames(cells_in_neighbourhood_df)[c(-1)]
  
  df <- reshape2::melt(cells_in_neighbourhood_df, measure.vars = target_cell_types)
  colnames(df) <- c("ref_cell_id", "tar_cell_type", "count")
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  tar_cell_type <- count <- NULL
  
  fig <- ggplot(df, aes(x = tar_cell_type, y = count)) + 
    geom_violin() +
    facet_wrap(~tar_cell_type, scales=scales, strip.position="bottom") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(title=paste("Cells in neighbourhood of", reference_cell_type, "cells"), x = "Target cell type", y = "Number of cells") +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult= 1), colour = "red")
  
  message("Plots show mean ± sd")
  
  return(fig)
}
plot_cells3D <- function(spe,
                         plot_cell_types = NULL,
                         plot_colours = NULL,
                         feature_colname = "Cell.Type") {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (!is.null(plot_cell_types)) {
    if(!is.character(plot_cell_types)) {
      stop("`plot_cell_types` is not a character vector.")
    }
  } 
  if (!is.null(plot_colours)) {
    non_colours <- plot_colours[which(!(sapply(plot_colours, function(X) {
      tryCatch(is.matrix(col2rgb(X)), 
               error = function(e) FALSE)
    })))]
    if (length(non_colours) > 0) {
      stop(paste("The following plot_colours are not colours:\n   ",
                 paste(non_colours, collapse = ", ")))
    } 
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste(feature_colname, "is not a valid column in your spe object."))
  }
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(spe), "Cell.Type" = spe[[feature_colname]])
  
  ## If no cell types chosen, use all cell types found in data frame
  if (is.null(plot_cell_types)) {
    warning("plot_cell_types not specified, all cell types found in the spe object will be used.")
    plot_cell_types <- unique(df[["Cell.Type"]])
  }
  ## If no colours inputted, use rainbow palette
  if (is.null(plot_colours)) {
    warning("plot_colours not specified, rainbow palette will be used.")
    plot_colours <- rainbow(length(plot_cell_types))
  }
  ## User inputs mismatching cell types and colours
  if (length(plot_cell_types) != length(plot_colours)) {
    stop("Length of plot_cell_types is not equal to length of plot_colours")
  }
  
  ## If cell types have been chosen, check they are found in the spe object
  spe_cell_types <- unique(spe[[feature_colname]])
  unknown_cell_types <- setdiff(plot_cell_types, spe_cell_types)
  
  if (length(unknown_cell_types) == length(plot_cell_types)) {
    stop("None of the plot_cell_types are found in the spe object")
  }
  
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following plot_cell_types are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    plot_colours <- plot_colours[which(plot_cell_types %in% spe_cell_types)]
    plot_cell_types <- intersect(plot_cell_types, spe_cell_types)
  }
  
  ## Factor for feature column
  df[, "Cell.Type"] <- factor(df[, "Cell.Type"],
                              levels = plot_cell_types)
  
  ## Plot
  fig <- plot_ly(df,
                 type = "scatter3d",
                 mode = 'markers',
                 x = ~Cell.X.Position,
                 y = ~Cell.Y.Position,
                 z = ~Cell.Z.Position,
                 color = ~Cell.Type,
                 colors = plot_colours,
                 marker = list(size = 2))
  
  fig <- fig %>% layout(scene = list(xaxis = list(title = 'x', showgrid = T, showaxeslabels = F, showticklabels = T, gridwidth = 5, 
                                                  titlefont = list(size = 20), tickfont = list(size = 15)),
                                     yaxis = list(title = 'y', showgrid = T, showaxeslabels = F, showticklabels = T, gridwidth = 5,
                                                  titlefont = list(size = 20), tickfont = list(size = 15)),
                                     zaxis = list(title = 'z', showgrid = T, showaxeslabels = F, showticklabels = T, gridwidth = 5,
                                                  titlefont = list(size = 20), tickfont = list(size = 15))))
  
  return(fig)
}
plot_co_occurrence_gradient3D <- function(co_occurrence_gradient_df) {
  
  target_cell_types <- colnames(co_occurrence_gradient_df)
  target_cell_types <- target_cell_types[!target_cell_types %in% c("reference", "radius")]
  
  co_occurrence_gradient_df$expected <- 1
  
  plot_result <- reshape2::melt(co_occurrence_gradient_df, "radius", c(target_cell_types, "expected"))
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Co-occurrence gradient", x = "Radius", y = "Co-occurrence value") +
    scale_colour_discrete(name = "") +
    theme_bw()
  
  return(fig) 
}
plot_cross_G_gradient3D <- function(cross_G_gradient_df, reference_cell_type = NULL, target_cell_type = NULL) {
  
  plot_result <- reshape2::melt(cross_G_gradient_df, "radius", c("observed_cross_G", "expected_cross_G"))
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Cross G-function gradient", x = "Radius", y = "Cross G-function value") +
    scale_colour_discrete(name = "", labels = c("Observed cross G", "Expected CSR cross G")) +
    theme_bw()
  
  if (!is.null(reference_cell_type) && !is.null(target_cell_type)) {
    fig <- fig + labs(subtitle = paste("Reference: ", reference_cell_type, ", Target: ", target_cell_type, sep = ""))
  }
  
  return(fig) 
}
plot_cross_K_gradient_ratio3D <- function(cross_K_gradient_df) {
  
  target_cell_types <- colnames(cross_K_gradient_df)[!colnames(cross_K_gradient_df) %in% c("reference", "expected", "radius")]
  
  # Normalize columns by 'expected'
  for (target_cell_type in target_cell_types) {
    cross_K_gradient_df[[target_cell_type]] <- cross_K_gradient_df[[target_cell_type]] / cross_K_gradient_df$expected
  }
  cross_K_gradient_df$expected <- 1
  
  plot_result <- reshape2::melt(cross_K_gradient_df, "radius", c(target_cell_types, "expected"))
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Cross K-function gradient ratio", x = "Radius", y = "Cross K-function ratio") +
    scale_colour_discrete(name = "") +
    theme_bw()
  
  return(fig) 
}
plot_cross_K_gradient3D <- function(cross_K_gradient_df) {
  
  target_cell_types <- colnames(cross_K_gradient_df)[!colnames(cross_K_gradient_df) %in% c("reference", "expected", "radius")]
  
  plot_result <- reshape2::melt(cross_K_gradient_df, "radius", c(target_cell_types, "expected"))
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Cross K-function gradient", x = "Radius", y = "Cross K-function value") +
    scale_colour_discrete(name = "") +
    theme_bw()
  
  return(fig) 
}
plot_cross_L_gradient_ratio3D <- function(cross_L_gradient_df) {
  
  target_cell_types <- colnames(cross_L_gradient_df)[!colnames(cross_L_gradient_df) %in% c("reference", "expected", "radius")]
  
  # Normalize columns by 'expected'
  for (target_cell_type in target_cell_types) {
    cross_L_gradient_df[[target_cell_type]] <- cross_L_gradient_df[[target_cell_type]] / cross_L_gradient_df$expected
  }
  cross_L_gradient_df$expected <- 1
  
  plot_result <- reshape2::melt(cross_L_gradient_df, "radius", c(target_cell_types, "expected"))
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Cross L-function gradient ratio", x = "Radius", y = "Cross L-function ratio") +
    scale_colour_discrete(name = "") +
    theme_bw()
  
  return(fig) 
}
plot_cross_L_gradient3D <- function(cross_L_gradient_df) {
  
  target_cell_types <- colnames(cross_L_gradient_df)[!colnames(cross_L_gradient_df) %in% c("reference", "expected", "radius")]
  
  plot_result <- reshape2::melt(cross_L_gradient_df, "radius", c(target_cell_types, "expected"))
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Cross L-function gradient", x = "Radius", y = "Cross L-function value") +
    scale_colour_discrete(name = "") +
    theme_bw()
  
  return(fig) 
}
## For scales parameter, use "free_x" or "free". "free_y" looks silly
plot_distances_between_cell_types_violin3D <- function(distances_df, scales = "free_x") {
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  pair <- distance <- NULL
  
  fig <- ggplot(distances_df, aes(x = pair, y = distance)) + 
    geom_violin() +
    facet_wrap(~pair, scales=scales, strip.position="bottom") +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
    labs(title="Cell distances", x = "Reference/Target pair", y = "Distance") +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult= 1), colour = "red")
  
  message("Plots show mean ± sd")
  
  return(fig)
}
plot_entropy_gradient3D <- function(entropy_gradient_df, expected_entropy = NULL, reference_cell_type = NULL, target_cell_types = NULL) {
  
  plot_result <- entropy_gradient_df
  
  if (!is.null(expected_entropy)) {
    if (!is.numeric(expected_entropy) || length(expected_entropy) != 1) stop("Please enter a single number for expected_entropy")
    plot_result$expected_entropy <- expected_entropy
    plot_result <- reshape2::melt(plot_result, "radius", c("entropy", "expected_entropy"))
    labels <- c("Observed entropy", "Expected CSR entropy")
  }
  else {
    plot_result <- reshape2::melt(plot_result, "radius", c("entropy"))
    labels <- c("Observed entropy")
  }
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Average entropy gradient", x = "Radius", y = "Entropy") +
    scale_colour_discrete(name = "", labels = labels) +
    theme_bw()
  
  if (!is.null(reference_cell_type) && !is.null(target_cell_types)) {
    fig <- fig + labs(subtitle = paste("Reference: ", reference_cell_type, ", Target: ", target_cell_types, sep = ""))
  }
  
  return(fig)
}
plot_grid_based_clusters3D <- function(spe_with_grid, 
                                       plot_cell_types = NULL,
                                       plot_colours = NULL,
                                       feature_colname = "Cell.Type") {
  
  # Check input parameters
  if (class(spe_with_grid) != "SpatialExperiment") {
    stop("`spe_with_grid` is not a SpatialExperiment object.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe_with_grid[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  
  ## If no cell types chosen, use all cell types found in data frame
  if (is.null(plot_cell_types)) plot_cell_types <- unique(spe_with_grid[[feature_colname]])
  
  ## If cell types have been chosen, check they are found in the spe object
  unknown_cell_types <- setdiff(plot_cell_types, spe_with_grid[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following plot_cell_types are not found in the spe object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  
  ## If no colours inputted, use rainbow palette
  if (is.null(plot_colours)) plot_colours <- rainbow(length(plot_cell_types))
  
  ## User inputs mismatching cell types and colours
  if (length(plot_cell_types) != length(plot_colours)) stop("Length of plot_cell_types is not equal to length of plot_colours")
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(spe_with_grid), colData(spe_with_grid))
  
  ## Factor for feature column
  df[[feature_colname]] <- factor(df[[feature_colname]], levels = plot_cell_types)
  
  ## Add points to fig
  fig <- plot_ly() %>%
    add_trace(
      data = df,
      type = "scatter3d",
      mode = 'markers',
      x = ~Cell.X.Position,
      y = ~Cell.Y.Position,
      z = ~Cell.Z.Position,
      marker = list(size = 2),
      color = ~.data[[feature_colname]],
      colors = plot_colours
    ) %>% 
    layout(scene = list(xaxis = list(title = 'x'),
                        yaxis = list(title = 'y'),
                        zaxis = list(title = 'z')))
  
  # Get number of grid-based clusters
  n_grid_based_clusters <- length(spe_with_grid@metadata[["grid_prisms"]])
  
  faces <- data.frame(edge1 = c(1, 1, 1, 1, 1, 1, 8, 8, 8, 8, 8, 8),
                      edge2 = c(2, 5, 2, 3, 3, 5, 6, 4 ,7, 6, 7, 4),
                      edge3 = c(6, 6, 4, 4, 7, 7, 2, 2, 5, 5, 3, 3))
  grid_based_colours <- rainbow(n_grid_based_clusters)
  
  ## Add grid-based clusters to fig, one by one  
  for (i in seq(n_grid_based_clusters)) {
    
    grid_based_cluster <- spe_with_grid@metadata[["grid_prisms"]][[i]]
    
    for (j in seq(nrow(grid_based_cluster))) {
      
      x <- grid_based_cluster$x[j]
      y <- grid_based_cluster$y[j]
      z <- grid_based_cluster$z[j]
      l <- grid_based_cluster$l[j]
      w <- grid_based_cluster$w[j]
      h <- grid_based_cluster$h[j]
      vertices <- data.frame(x = c(x, x + l, x, x + l, x, x + l, x, x + l),
                             y = c(y, y, y + w, y + w, y, y, y + w, y + w),
                             z = c(z, z, z, z, z + h, z + h, z + h, z + h))
      
      fig <- fig %>%
        add_trace(
          type = 'mesh3d',
          x = vertices[, 1], 
          y = vertices[, 2], 
          z = vertices[, 3],
          i = faces[, 1] - 1, 
          j = faces[, 2] - 1, 
          k = faces[, 3] - 1,
          opacity = 0.2,
          facecolor = rep(grid_based_colours[i], 12) # Always 12 faces per grid prism
        )      
    }
  }
  
  return(fig)
}
plot_grid_metrics_continuous3D <- function(grid_metrics, metric_colname) {
  
  ## Check input parameters
  if (!(is.character(metric_colname) && metric_colname %in% c("proportion", "entropy"))) {
    stop("`metric_colname` is not 'proportion' or 'entropy'.")
  }
  if (is.null(grid_metrics[[metric_colname]])) {
    stop("`metric_colname` is not a column in `grid_metrics`.")
  }
  
  ## Color of each dot is related to its entropy
  pal <- colorRampPalette(hcl.colors(n = 5, palette = "Red-Blue", rev = TRUE))
  
  ## Add size column and for NA entropy values, make the size small
  grid_metrics$size <- ifelse(is.na(grid_metrics[[metric_colname]]), 3, 10)
  
  fig <- plot_ly(grid_metrics,
                 type = "scatter3d",
                 mode = 'markers',
                 x = ~x_coord,
                 y = ~y_coord,
                 z = ~z_coord,
                 color = as.formula(paste0('~', metric_colname)),
                 colors = pal(nrow(grid_metrics)),
                 marker = list(size = ~size),
                 symbol = 1,
                 symbols = "square")
  
  fig <- fig %>% layout(scene = list(xaxis = list(title = 'x'),
                                     yaxis = list(title = 'y'),
                                     zaxis = list(title = 'z')))
  
  return(fig)
}
plot_grid_metrics_discrete3D <- function(grid_metrics, metric_colname) {
  
  ## Check input parameters
  if (!(is.character(metric_colname) && metric_colname %in% c("proportion", "entropy"))) {
    stop("`metric_colname` is not 'proportion' or 'entropy'.")
  }
  if (is.null(grid_metrics[[metric_colname]])) {
    stop("`metric_colname` is not a column in `grid_metrics`.")
  }
  
  ## Define low, medium and high categories
  # Low: between 0 and 1/3
  # Medium: between 1/3 and 2/3
  # High: between 2/3 and 1
  
  grid_metrics$rank <- ifelse(is.na(grid_metrics[[metric_colname]]), "na",
                              ifelse(grid_metrics[[metric_colname]] < 1/3, "low",
                                     ifelse(grid_metrics[[metric_colname]] < 2/3, "medium", "high")))
  grid_metrics$rank <- factor(grid_metrics$rank, c("low", "medium", "high", "na"))
  
  fig <- plot_ly(grid_metrics,
                 type = "scatter3d",
                 mode = 'markers',
                 x = ~x_coord,
                 y = ~y_coord,
                 z = ~z_coord,
                 color = ~rank,
                 colors = c("#AEB6E5", "#BC6EB9", "#A93154", "gray"),
                 symbol = 1,
                 symbols = "square",
                 marker = list(size = 4))
  
  fig <- fig %>% layout(scene = list(xaxis = list(title = 'x'),
                                     yaxis = list(title = 'y'),
                                     zaxis = list(title = 'z')))
  return(fig)
}
plot_mixing_scores_gradient3D <- function(mixing_scores_gradient_df, metric = "MS") {
  
  if (!metric %in% c("MS", "NMS")) {
    stop("'metric' should be 'MS' or 'NMS', for mixing score and normalised mixing score respectively.")
  }
  
  if (metric == "NMS") {
    plot_result <- mixing_scores_gradient_df
    plot_result$expected_normalised_mixing_score <- 1
    plot_result <- reshape2::melt(plot_result, "radius", c("normalised_mixing_score", "expected_normalised_mixing_score"))
    
    fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
      geom_line() +
      labs(title = "Normalised mixing score (NMS) gradient", 
           subtitle = paste("Reference: ", mixing_scores_gradient_df$ref_cell_type[1], ", Target: ", mixing_scores_gradient_df$tar_cell_type[1], sep = ""), 
           x = "Radius", y = "NMS") +
      scale_colour_discrete(name = "", labels = c("Observed NMS", "Expected CSR NMS")) +
      theme_bw() 
  }
  else if (metric == "MS") {
    plot_result <- mixing_scores_gradient_df
    n_tar_cells <- plot_result$n_tar_cells[1]
    n_ref_cells <- plot_result$n_ref_cells[1]
    plot_result$expected_mixing_score <- n_tar_cells * n_ref_cells / ((n_ref_cells - 1) * n_ref_cells / 2)
    plot_result <- reshape2::melt(plot_result, "radius", c("mixing_score", "expected_mixing_score"))
    
    fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
      geom_line() +
      labs(title = "Mixing score (MS) gradient", 
           subtitle = paste("Reference: ", mixing_scores_gradient_df$ref_cell_type[1], ", Target: ", mixing_scores_gradient_df$tar_cell_type[1], sep = ""), 
           x = "Radius", y = "MS") +
      scale_colour_discrete(name = "", labels = c("Observed MS", "Expected CSR MS  ")) +
      theme_bw()  
  }
  return(fig)
}

summarise_cells_in_neighbourhood3D <- function(cells_in_neighbourhood_df) {
  
  ## Target cell types will be all the columns except the first column
  target_cell_types <- colnames(cells_in_neighbourhood_df)[c(-1)]
  
  ## Set up data frame for summarised_results list
  df <- data.frame(row.names = c("mean", "min", "max", "median", "st_dev"))
  
  for (target_cell_type in target_cell_types) {
    
    ## Get statistical measures for each target cell type
    target_cell_type_values <- cells_in_neighbourhood_df[[target_cell_type]]
    df[[target_cell_type]] <- c(mean(target_cell_type_values),
                                min(target_cell_type_values),
                                max(target_cell_type_values),
                                median(target_cell_type_values),
                                sd(target_cell_type_values))
    
  }
  
  return(data.frame(t(df)))
}
summarise_distances_between_cell_types3D <- function(distances_df) {
  
  pair <- distance <- NULL
  
  # summarise the results
  distances_df_summarised <- distances_df %>% 
    dplyr::group_by(pair) %>%
    dplyr::summarise(mean(distance), 
                     min(distance), 
                     max(distance),
                     stats::median(distance), 
                     stats::sd(distance))
  
  distances_df_summarised <- data.frame(distances_df_summarised)
  
  colnames(distances_df_summarised) <- c("pair", 
                                         "mean", 
                                         "min", 
                                         "max", 
                                         "median", 
                                         "std_dev")
  
  for (i in seq(nrow(distances_df_summarised))) {
    # Get cell_types for each pair
    cell_types <- strsplit(distances_df_summarised[i,"pair"], "/")[[1]]
    
    distances_df_summarised[i, "reference"] <- cell_types[1]
    distances_df_summarised[i, "target"] <- cell_types[2]
  }
  
  return(distances_df_summarised)
}


calculate_all_gradient_cc_metrics2D <- function(spe, 
                                                reference_cell_type, 
                                                target_cell_types, 
                                                radii, 
                                                feature_colname = "Cell.Type", 
                                                plot_image = T) {
  
  # Define constants
  cross_K_df_colnames <- c("reference",
                           "expected",
                           target_cell_types)
  mixing_score_df_colnames <- c("ref_cell_type", 
                                "tar_cell_type", 
                                "n_ref_cells",
                                "n_tar_cells", 
                                "n_ref_tar_interactions",
                                "n_ref_ref_interactions", 
                                "mixing_score", 
                                "normalised_mixing_score")
  cross_G_df_colnames <- c("observed_cross_G",
                           "expected_cross_G")
  co_occurrence_df_colnames <- c("reference",
                                 target_cell_types)
  
  ## Define result
  result <- list("mixing_score" = list(),
                 "cells_in_neighbourhood" = data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types))),
                 "cells_in_neighbourhood_proportion" = data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types))),
                 "entropy" = data.frame(matrix(nrow = length(radii), ncol = 1)),
                 "cross_K" = data.frame(matrix(nrow = length(radii), ncol = length(cross_K_df_colnames))),
                 "cross_L" = data.frame(matrix(nrow = length(radii), ncol = length(cross_K_df_colnames))),
                 "cross_G" = list(),
                 "co_occurrence" = data.frame(matrix(nrow = length(radii), ncol = length(co_occurrence_df_colnames))))
  colnames(result[["cells_in_neighbourhood"]]) <- target_cell_types
  colnames(result[["cells_in_neighbourhood_proportion"]]) <- target_cell_types
  colnames(result[["entropy"]]) <- "entropy"
  colnames(result[["cross_K"]]) <- cross_K_df_colnames
  colnames(result[["cross_L"]]) <- cross_K_df_colnames
  colnames(result[["co_occurrence"]]) <- co_occurrence_df_colnames
  
  # Define individual data frames for mixing_score and cross_G
  for (target_cell_type in target_cell_types) {
    if (reference_cell_type != target_cell_type) {
      result[["mixing_score"]][[target_cell_type]] <- data.frame(matrix(nrow = length(radii), ncol = length(mixing_score_df_colnames)))
      colnames(result[["mixing_score"]][[target_cell_type]]) <- mixing_score_df_colnames
    }
    result[["cross_G"]][[target_cell_type]] <- data.frame(matrix(nrow = length(radii), ncol = length(cross_G_df_colnames)))
    colnames(result[["cross_G"]][[target_cell_type]]) <- cross_G_df_colnames
  }
  
  # Get gradient results for each metric
  for (i in seq(length(radii))) {
    df <- calculate_all_single_radius_cc_metrics2D(spe,
                                                   reference_cell_type,
                                                   target_cell_types,
                                                   radii[i],
                                                   feature_colname)
    
    if (is.null(df)) return(NULL)
    
    df[["cells_in_neighbourhood"]]$ref_cell_id <- NULL
    
    result[["cells_in_neighbourhood"]][i, ] <- apply(df[["cells_in_neighbourhood"]], 2, mean)
    result[["cells_in_neighbourhood_proportion"]][i, ] <- apply(df[["cells_in_neighbourhood_proportion"]][ , paste(target_cell_types, "_prop", sep = "")], 2, mean, na.rm = T)
    result[["entropy"]][i, "entropy"] <- mean(df[["entropy"]]$entropy, na.rm = T)
    result[["cross_K"]][i, ] <- df[["cross_K"]]
    result[["cross_L"]][i, ] <- df[["cross_L"]]
    result[["co_occurrence"]][i, ] <- df[["co_occurrence"]]
    
    for (target_cell_type in names(df[["mixing_score"]])) {
      result[["mixing_score"]][[target_cell_type]][i, ] <- df[["mixing_score"]][[target_cell_type]]
    }
    for (target_cell_type in names(df[["cross_G"]])) {
      result[["cross_G"]][[target_cell_type]][i, ] <- df[["cross_G"]][[target_cell_type]]
    }
  }
  
  # Add radius column to each data frame
  result[["cells_in_neighbourhood"]]$radius <- radii
  result[["cells_in_neighbourhood_proportion"]]$radius <- radii
  result[["entropy"]]$radius <- radii
  result[["cross_K"]]$radius <- radii
  result[["cross_L"]]$radius <- radii
  result[["co_occurrence"]]$radius <- radii
  for (target_cell_type in names(df[["mixing_score"]])) {
    result[["mixing_score"]][[target_cell_type]]$radius <- radii
  }
  for (target_cell_type in names(df[["cross_G"]])) {
    result[["cross_G"]][[target_cell_type]]$radius <- radii
  }
  
  
  ## Plot
  if (plot_image) {
    fig_ACIN <- plot_cells_in_neighbourhood_gradient2D(result[["cells_in_neighbourhood"]], reference_cell_type)
    methods::show(fig_ACIN)
    
    fig_ACINP <- plot_cells_in_neighbourhood_proportions_gradient2D(result[["cells_in_neighbourhood_proportion"]], reference_cell_type)
    methods::show(fig_ACINP)
    
    expected_entropy <- calculate_entropy_background2D(spe, target_cell_types, feature_colname)
    fig_AE <- plot_entropy_gradient2D(result[["entropy"]], expected_entropy, reference_cell_type, target_cell_types)
    methods::show(fig_AE)
    
    for (target_cell_type in names(result[["mixing_score"]])) {
      fig_NMS <- plot_mixing_scores_gradient2D(result[["mixing_score"]][[target_cell_type]], "NMS")
      fig_MS <- plot_mixing_scores_gradient2D(result[["mixing_score"]][[target_cell_type]], "MS")
      fig_NMS_MS <- plot_grid(fig_NMS, fig_MS, nrow = 2)
      methods::show(fig_NMS_MS)
    }
    fig_CK <- plot_cross_K_gradient2D(result[["cross_K"]])
    fig_CKR <- plot_cross_K_gradient_ratio2D(result[["cross_K"]])
    fig_CK_CKR <- plot_grid(fig_CK, fig_CKR, nrow = 2)
    methods::show(fig_CK_CKR)
    
    fig_CL <- plot_cross_L_gradient2D(result[["cross_L"]])
    fig_CLR <- plot_cross_L_gradient_ratio2D(result[["cross_L"]])
    fig_CL_CLR <- plot_grid(fig_CL, fig_CLR, nrow = 2)
    methods::show(fig_CL_CLR)
    
    for (target_cell_type in names(result[["cross_G"]])) {
      fig_CG <- plot_cross_G_gradient2D(result[["cross_G"]][[target_cell_type]], reference_cell_type, target_cell_type)
      methods::show(fig_CG)
    }
    
    fig_co_occ <- plot_co_occurrence_gradient2D(result[["co_occurrence"]])
    methods::show(fig_co_occ)
  }
  
  return(result)
}


### Calculate all single radius cell-colocalisation metrics
# If a function only requires one target cell type, iterate through each cell type in target_cell_types, else use all target_cell_types

calculate_all_single_radius_cc_metrics2D <- function(spe, 
                                                     reference_cell_type, 
                                                     target_cell_types, 
                                                     radius, 
                                                     feature_colname = "Cell.Type") {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (!(is.character(reference_cell_type) && length(reference_cell_type) == 1)) {
    stop("`reference_cell_type` is not a character.")
  }
  if (!is.character(target_cell_types)) {
    stop("`target_cell_types` is not a character vector.")
  }
  if (!(is.numeric(radius) && length(radius) == 1 && radius > 0)) {
    stop(paste(radius, " is not a positive numeric."))
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  
  ## For reference_cell_type, check it is found in the spe object
  if (!(reference_cell_type %in% spe[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spe object"))
    return(NULL)
  }
  ## For target_cell_types, check they are found in the spe object
  unknown_cell_types <- setdiff(target_cell_types, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
  }
  
  # Define result
  result <- list("cells_in_neighbourhood" = list(),
                 "cells_in_neighbourhood_proportion" = list(),
                 "entropy" = list(),
                 "mixing_score" = list(),
                 "cross_K" = list(),
                 "cross_L" = list(),
                 "cross_G" = list(),
                 "co_occurrence" = list())
  
  # Define other constants
  mixing_score_df_colnames <- c("ref_cell_type", 
                                "tar_cell_type", 
                                "n_ref_cells",
                                "n_tar_cells", 
                                "n_ref_tar_interactions",
                                "n_ref_ref_interactions", 
                                "mixing_score", 
                                "normalised_mixing_score")
  
  cross_K_df_colnames <- c("reference",
                           "expected",
                           target_cell_types)
  cross_G_df_colnames <- c("observed_cross_G",
                           "expected_cross_G")
  co_occurrence_df_colnames <- c("reference",
                                 target_cell_types)
  
  # Get rough dimensions of window for cross_K
  spe_coords <- data.frame(spatialCoords(spe))
  length <- round(max(spe_coords$Cell.X.Position) - min(spe_coords$Cell.X.Position))
  width  <- round(max(spe_coords$Cell.Y.Position) - min(spe_coords$Cell.Y.Position))
  ## Get volume of the window the cells are in
  volume <- length * width
  
  
  
  # All single radius cc metrics stem from calculate_entropy2D function
  entropy_df <- calculate_entropy2D(spe, 
                                    reference_cell_type, 
                                    target_cell_types, 
                                    radius, 
                                    feature_colname)  
  
  ## Cells in neighbourhood ----------
  result[["cells_in_neighbourhood"]] <- entropy_df[ , c("ref_cell_id", target_cell_types)]
  
  ## Cells in neighbourhood proportion ----------
  result[["cells_in_neighbourhood_proportion"]] <- entropy_df[ , c("ref_cell_id", target_cell_types, paste(target_cell_types, "_prop", sep = ""))]
  
  ## Entropy --------------
  result[["entropy"]] <- entropy_df
  
  ## Mixing score -----------------
  for (target_cell_type in target_cell_types) {
    mixing_score_df <- data.frame(matrix(nrow = 1, ncol = length(mixing_score_df_colnames)))
    colnames(mixing_score_df) <- mixing_score_df_colnames
    mixing_score_df$ref_cell_type <- reference_cell_type
    
    # No need to fill in mixing_score_df if the reference and target cell is the same
    if (reference_cell_type != target_cell_type) {
      mixing_score_df$tar_cell_type <- target_cell_type
      mixing_score_df$n_ref_cells <- sum(spe[[feature_colname]] == reference_cell_type)
      mixing_score_df$n_tar_cells <- sum(spe[[feature_colname]] == target_cell_type)
      mixing_score_df$n_ref_tar_interactions <- sum(entropy_df[[target_cell_type]])
      mixing_score_df$n_ref_ref_interactions <- sum(entropy_df[[reference_cell_type]])
      mixing_score_df$mixing_score <- mixing_score_df$n_ref_tar_interactions / (0.5 * mixing_score_df$n_ref_ref_interactions)
      mixing_score_df$normalised_mixing_score <- 0.5 * mixing_score_df$mixing_score * mixing_score_df$n_ref_cells / mixing_score_df$n_tar_cell
      if (is.infinite(mixing_score_df$mixing_score)) mixing_score_df$mixing_score <- NA
      if (is.infinite(mixing_score_df$normalised_mixing_score)) mixing_score_df$normalised_mixing_score <- NA
      result[["mixing_score"]][[target_cell_type]] <- mixing_score_df
    }
  }
  
  ## Cross_K ---------------------
  cross_K_df <- data.frame(matrix(nrow = 1, ncol = length(cross_K_df_colnames)))
  colnames(cross_K_df) <- cross_K_df_colnames
  cross_K_df$reference <- reference_cell_type
  cross_K_df$expected <- (4/3) * pi * radius^3
  
  for (target_cell_type in target_cell_types) {
    cross_K_df[[target_cell_type]] <- (((volume * sum(entropy_df[[target_cell_type]])) / sum(spe[[feature_colname]] == reference_cell_type)) / sum(spe[[feature_colname]] == target_cell_type)) 
  }
  result[["cross_K"]] <- cross_K_df
  
  ## Cross_L ---------------------
  cross_L_df <- cross_K_df
  cross_L_df[ , c("expected", target_cell_types)] <- (cross_L_df[ , c("expected", target_cell_types)] / (4 * pi / 3)) ^ (1/3)
  result[["cross_L"]] <- cross_L_df
  
  ## Cross_G ---------------------
  for (target_cell_type in target_cell_types) {
    cross_G_df <- data.frame(matrix(nrow = 1, ncol = length(cross_G_df_colnames)))
    colnames(cross_G_df) <- cross_G_df_colnames
    
    reference_target_interactions <- entropy_df[[target_cell_type]]
    n_target_cells <- sum(spe[[feature_colname]] == target_cell_type)
    target_cell_type_intensity <- n_target_cells / volume
    observed_cross_G <- sum(reference_target_interactions != 0) / length(reference_target_interactions)
    expected_cross_G <- 1 - exp(-1 * target_cell_type_intensity * (4 / 3) * pi * radius^3)
    
    cross_G_df$observed_cross_G <- observed_cross_G
    cross_G_df$expected_cross_G <- expected_cross_G
    result[["cross_G"]][[target_cell_type]] <- cross_G_df
  }
  
  
  ## Co_occurrence ---------------
  all_cell_types <- unique(spe[[feature_colname]])
  cells_in_neighbourhood_proportions_df <- calculate_cells_in_neighbourhood_proportions2D(spe,
                                                                                          reference_cell_type,
                                                                                          all_cell_types,
                                                                                          radius,
                                                                                          feature_colname)
  
  co_occurrence_df <- data.frame(matrix(nrow = 1, ncol = length(co_occurrence_df_colnames)))
  colnames(co_occurrence_df) <- co_occurrence_df_colnames
  co_occurrence_df$reference <- reference_cell_type
  
  n_cells_in_spe <- length(spe[[feature_colname]])
  n_cells_in_reference_cell_type_radius <- sum(cells_in_neighbourhood_proportions_df$total)
  
  for (target_cell_type in target_cell_types) {
    n_target_cells_in_reference_cell_type_radius <- sum(cells_in_neighbourhood_proportions_df[[target_cell_type]])
    target_cell_type_proportion_in_reference_cell_type_radius <- n_target_cells_in_reference_cell_type_radius / n_cells_in_reference_cell_type_radius
    n_target_cells_in_spe <- sum(spe[[feature_colname]] == target_cell_type)
    target_cell_type_proportion_in_spe <- n_target_cells_in_spe / n_cells_in_spe
    target_cell_type_co_occurrence <- target_cell_type_proportion_in_reference_cell_type_radius / target_cell_type_proportion_in_spe
    
    co_occurrence_df[[target_cell_type]] <- target_cell_type_co_occurrence
  }
  result[["co_occurrence"]] <- co_occurrence_df
  
  return(result)
}

calculate_cell_proportion_grid_metrics2D <- function(spe, 
                                                     n_splits,
                                                     reference_cell_types,
                                                     target_cell_types,
                                                     feature_colname = "Cell.Type",
                                                     plot_image = TRUE) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  ## Check reference_cell_types are found in the spe object
  unknown_cell_types <- setdiff(reference_cell_types, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in reference_cell_types are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  ## Check target_cell_types are found in the spe object
  unknown_cell_types <- setdiff(target_cell_types, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  # Check if there is intersection between reference_cell_types and target_cell_types
  if (length(intersect(reference_cell_types, target_cell_types)) > 0) {
    stop("Cannot have same cells in both reference_cell_types and target_cell_types")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Add grid metrics to spe
  spe <- get_spe_grid_metrics2D(spe, n_splits, feature_colname)
  
  # Get grid_prism_cell_matrix from spe
  grid_prism_cell_matrix <- spe@metadata$grid_metrics$grid_prism_cell_matrix
  
  ## Define data frame which contains all results
  n_grid_prisms <- n_splits^3
  result <- data.frame(row.names = seq(n_grid_prisms))
  
  # Fill in the result data frame
  if (length(reference_cell_types) == 1) {
    result$reference <- grid_prism_cell_matrix[[reference_cell_types]]
  }
  else {
    result$reference <- rowSums(grid_prism_cell_matrix[ , reference_cell_types])
  }
  if (length(target_cell_types) == 1) {
    result$target <- grid_prism_cell_matrix[[target_cell_types]]
  }
  else {
    result$target <- rowSums(grid_prism_cell_matrix[ , target_cell_types])
  }
  result$total <- result$reference + result$target
  result$proportion <- result$target / result$total
  
  # Add grid_prism_coordinates info to result
  result <- cbind(result, spe@metadata$grid_metrics$grid_prism_coordinates)
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_metrics_continuous2D(result, "proportion")
    methods::show(fig)
  }
  
  return(result)
}

calculate_cell_proportions2D <- function(spe,
                                         cell_types_of_interest = NULL, 
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (ncol(spe) == 0) {
    stop("No cells found for calculating cell proportions.")
  }
  if (!(is.null(cell_types_of_interest) || is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Creates frequency/bar plot of all cell types in the entire image
  cell_proportions <- data.frame(table(spe[[feature_colname]]))
  names(cell_proportions) <- c("cell_type", 'frequency')
  
  # Only include cell types the user has chosen
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spe object
    unknown_cell_types <- setdiff(cell_types_of_interest, cell_proportions$cell_type)
    if (length(unknown_cell_types) != 0) {
      stop(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                 paste(unknown_cell_types, collapse = ", ")))
    }
    
    # Subset for cell types chosen by user
    cell_proportions <- cell_proportions[(cell_proportions$cell_type %in% cell_types_of_interest), ]
    
  }
  
  # Get frequency total for all cells
  cell_type_frequency_total <- sum(cell_proportions$frequency)
  
  # Get proportions and percentages
  cell_proportions$proportion <- cell_proportions$frequency / cell_type_frequency_total
  cell_proportions$percentage <- cell_proportions$proportion * 100
  
  # Order the cell types by proportion (highest cell proportion is first)
  cell_proportions <- cell_proportions[rev(order(cell_proportions$proportion)), ]
  rownames(cell_proportions) <- seq(nrow(cell_proportions))
  
  # Plot
  if (plot_image) {
    
    labels <- paste(round(cell_proportions$percentage, 1), "%", sep = "")
    
    fig <- ggplot(cell_proportions, aes(x = factor(cell_type, cell_type), y = percentage, fill = cell_type)) +
      geom_bar(stat='identity') + 
      theme_bw() +
      labs(title="Cell proportions", x = "Cell type", y = "Percentage") +
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "none") +
      geom_text(aes(label = labels), vjust = 0)
    
    methods::show(fig)
  }
  
  return(cell_proportions)
}
calculate_cells_in_neighbourhood_gradient2D <- function(spe, 
                                                        reference_cell_type, 
                                                        target_cell_types, 
                                                        radii, 
                                                        feature_colname = "Cell.Type",
                                                        plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types)))
  colnames(result) <- target_cell_types
  
  for (i in seq(length(radii))) {
    cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood2D(spe,
                                                                    reference_cell_type,
                                                                    target_cell_types,
                                                                    radii[i],
                                                                    feature_colname,
                                                                    FALSE,
                                                                    FALSE)
    
    if (is.null(cells_in_neighbourhood_df)) return(NULL)
    
    cells_in_neighbourhood_df$ref_cell_id <- NULL
    result[i, ] <- apply(cells_in_neighbourhood_df, 2, mean)
  }
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig <- plot_cells_in_neighbourhood_gradient2D(result, reference_cell_type)
    methods::show(fig)
  }
  
  return(result)
}

calculate_cells_in_neighbourhood_proportions_gradient2D <- function(spe, 
                                                                    reference_cell_type, 
                                                                    target_cell_types, 
                                                                    radii, 
                                                                    feature_colname = "Cell.Type",
                                                                    plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types)))
  colnames(result) <- target_cell_types
  
  for (i in seq(length(radii))) {
    cell_proportions_neighbourhood_proportions_df <- calculate_cells_in_neighbourhood_proportions2D(spe,
                                                                                                    reference_cell_type,
                                                                                                    target_cell_types,
                                                                                                    radii[i],
                                                                                                    feature_colname)
    
    if (is.null(cell_proportions_neighbourhood_proportions_df)) return(NULL)
    
    result[i, ] <- apply(cell_proportions_neighbourhood_proportions_df[ , paste(target_cell_types, "_prop", sep = "")], 2, mean, na.rm = T)
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  # Plot
  if (plot_image) {
    fig <- plot_cells_in_neighbourhood_proportions_gradient2D(result, reference_cell_type)
    methods::show(fig)
  }
  
  return(result)
}

calculate_cells_in_neighbourhood_proportions2D <- function(spe, 
                                                           reference_cell_type, 
                                                           target_cell_types, 
                                                           radius, 
                                                           feature_colname = "Cell.Type") {
  
  ## Get cells in neighbourhood df
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood2D(spe,
                                                                  reference_cell_type,
                                                                  target_cell_types,
                                                                  radius,
                                                                  feature_colname,
                                                                  FALSE,
                                                                  FALSE)
  
  if (is.null(cells_in_neighbourhood_df)) return(NULL)
  
  ## Get total number of target cells for each row (first column is the reference cell id column, so we exclude it)
  cells_in_neighbourhood_df$total <- apply(cells_in_neighbourhood_df[ , c(-1)], 1, sum)
  
  cells_in_neighbourhood_df[ , paste(target_cell_types, "_prop", sep = "")] <- cells_in_neighbourhood_df[ , target_cell_types] / cells_in_neighbourhood_df$total
  
  return(cells_in_neighbourhood_df)
}
calculate_cells_in_neighbourhood2D <- function(spe, 
                                               reference_cell_type, 
                                               target_cell_types, 
                                               radius, 
                                               feature_colname = "Cell.Type",
                                               show_summary = TRUE,
                                               plot_image = TRUE) {
  
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (!(is.character(reference_cell_type) && length(reference_cell_type) == 1)) {
    stop("`reference_cell_type` is not a character.")
  }
  if (!is.character(target_cell_types)) {
    stop("`target_cell_types` is not a character vector.")
  }
  if (!(is.numeric(radius) && length(radius) == 1 && radius > 0)) {
    stop("`radius` is not a positive numeric.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(show_summary)) {
    stop("`show_summary` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  ## For reference_cell_type, check it is found in the spe object
  if (!(reference_cell_type %in% spe[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spe object"))
    return(NULL)
  }
  ## For target_cell_types, check they are found in the spe object
  unknown_cell_types <- setdiff(target_cell_types, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
  }
  
  if (is.null(spe[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spe")
    spe$Cell.ID <- paste("Cell", seq(ncol(spe)), sep = "_")
  }  
  
  # Get spe coords
  spe_coords <- data.frame(spatialCoords(spe))
  
  # Get reference_cell_type coords
  reference_cell_type_coords <- spe_coords[spe[[feature_colname]] == reference_cell_type, ]
  
  result <- data.frame(matrix(nrow = nrow(reference_cell_type_coords), ncol = 0))
  
  for (target_cell_type in target_cell_types) {
    
    if (sum(spe[[feature_colname]] == target_cell_type) == 0) {
      result[[target_cell_type]] <- NA
      next
    }
    
    ## Get target_cell_type coords
    target_cell_type_coords <- spe_coords[spe[[feature_colname]] == target_cell_type, ]
    
    ## Determine number of target cells specified distance for each reference cell
    ref_tar_result <- dbscan::frNN(target_cell_type_coords, 
                                   eps = radius,
                                   query = reference_cell_type_coords, 
                                   sort = FALSE)
    
    n_targets <- rapply(ref_tar_result$id, length)
    
    
    # Don't want to include the reference cell as one of the target cells
    if (reference_cell_type == target_cell_type) n_targets <- n_targets - 1
    
    ## Add to data frame
    result[[target_cell_type]] <- n_targets
  }
  
  result <- data.frame(ref_cell_id = spe$Cell.ID[spe[[feature_colname]] == reference_cell_type], result)
  
  ## Print summary
  if (show_summary) {
    print(summarise_cells_in_neighbourhood2D(result))    
  }
  
  ## Plot
  if (plot_image) {
    fig <- plot_cells_in_neighbourhood_violin2D(result, reference_cell_type)
    methods::show(fig)
  }
  
  return(result)
}

calculate_co_occurrence_gradient2D <- function(spe, 
                                               reference_cell_type, 
                                               target_cell_types, 
                                               radii, 
                                               feature_colname = "Cell.Type",
                                               plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types) + 1))
  colnames(result) <- c("reference", target_cell_types)
  
  for (i in seq(length(radii))) {
    co_occurrence_df <- calculate_co_occurrence2D(spe,
                                                  reference_cell_type,
                                                  target_cell_types,
                                                  radii[i],
                                                  feature_colname)
    
    result[i, ] <- co_occurrence_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig <- plot_co_occurrence_gradient2D(result)
    methods::show(fig)
  }
  
  return(result)
}
calculate_co_occurrence2D <- function(spe, 
                                      reference_cell_type, 
                                      target_cell_types, 
                                      radius, 
                                      feature_colname = "Cell.Type") {
  
  # Get all cell types in spe
  all_cell_types <- unique(spe[[feature_colname]])
  
  cells_in_neighbourhood_proportions_df <- calculate_cells_in_neighbourhood_proportions2D(spe,
                                                                                          reference_cell_type,
                                                                                          all_cell_types,
                                                                                          radius,
                                                                                          feature_colname)
  
  result <- data.frame(reference = reference_cell_type)
  
  # Get total number of cells in spe
  n_cells_in_spe <- length(spe[[feature_colname]])
  
  # Get total number of cells in radius around reference cell type
  n_cells_in_reference_cell_type_radius <- sum(cells_in_neighbourhood_proportions_df$total)
  
  for (target_cell_type in target_cell_types) {
    
    # Get total number of target cells in radius around reference cell type
    n_target_cells_in_reference_cell_type_radius <- sum(cells_in_neighbourhood_proportions_df[[target_cell_type]])
    
    # Get proportion of target cells in radius around reference cell type
    target_cell_type_proportion_in_reference_cell_type_radius <- n_target_cells_in_reference_cell_type_radius / n_cells_in_reference_cell_type_radius
    
    # Get proportion of target cell type in spe
    n_target_cells_in_spe <- sum(spe[[feature_colname]] == target_cell_type)
    target_cell_type_proportion_in_spe <- n_target_cells_in_spe / n_cells_in_spe
    
    # Get co-occurence value for taget cell type
    target_cell_type_co_occurrence <- target_cell_type_proportion_in_reference_cell_type_radius / target_cell_type_proportion_in_spe
    
    # Add to result data frame
    result[[target_cell_type]] <- target_cell_type_co_occurrence
  }
  
  return(result)
}
calculate_cross_G_gradient2D <- function(spe, 
                                         reference_cell_type, 
                                         target_cell_type, 
                                         radii, 
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = 2))
  colnames(result) <- c("observed_cross_G", 
                        "expected_cross_G")
  
  for (i in seq(length(radii))) {
    cross_G_df <- calculate_cross_G2D(spe,
                                      reference_cell_type,
                                      target_cell_type,
                                      radii[i],
                                      feature_colname)
    
    result[i, ] <- cross_G_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig <- plot_cross_G_gradient2D(result, reference_cell_type, target_cell_type)
    methods::show(fig)
  }
  
  return(result)
}
calculate_cross_G2D <- function(spe,
                                reference_cell_type,
                                target_cell_type,
                                radius,
                                feature_colname = "Cell.Type") {
  
  ### Calculate the observed cross_G
  # Get the number of target cells in the radius around each reference cell
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood2D(spe,
                                                                  reference_cell_type,
                                                                  target_cell_type,
                                                                  radius,
                                                                  feature_colname,
                                                                  show_summary = FALSE,
                                                                  plot_image = FALSE)
  
  reference_target_interactions <- cells_in_neighbourhood_df[[target_cell_type]]
  
  # cross_G: essentially the proportion of reference cells with at least 1 target cell within the chosen radius.
  observed_cross_G <- sum(reference_target_interactions != 0) / length(reference_target_interactions)
  
  ### Calculate the expected cross_G
  # Get rough dimensions of the window the points are in
  spe_coords <- data.frame(spatialCoords(spe))
  
  length <- round(max(spe_coords$Cell.X.Position) - min(spe_coords$Cell.X.Position))
  width  <- round(max(spe_coords$Cell.Y.Position) - min(spe_coords$Cell.Y.Position))
  
  # Get volume of the window the cells are in
  volume <- length * width
  
  # Get the number of target cells
  n_target_cells <- sum(spe[[feature_colname]] == target_cell_type)
  
  # Get target_cell_type intensity (density)
  target_cell_type_intensity <- n_target_cells / volume
  
  # Apply formula
  expected_cross_G <- 1 - exp(-1 * target_cell_type_intensity * pi * radius^2)
  
  result <- data.frame(observed_cross_G = observed_cross_G,
                       expected_cross_G = expected_cross_G)
  
  return(result)
}
calculate_cross_K_gradient2D <- function(spe, 
                                         reference_cell_type, 
                                         target_cell_types, 
                                         radii, 
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = 2 + length(target_cell_types)))
  colnames(result) <- c("reference", "expected", target_cell_types)
  
  for (i in seq(length(radii))) {
    cross_K_df <- calculate_cross_K2D(spe,
                                      reference_cell_type,
                                      target_cell_types,
                                      radii[i],
                                      feature_colname)
    
    result[i, ] <- cross_K_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig1 <- plot_cross_K_gradient2D(result)
    fig2 <- plot_cross_K_gradient_ratio2D(result)
    
    combined_fig <- plot_grid(fig1, fig2, nrow = 2)
    methods::show(combined_fig)
  }
  
  return(result)
}
calculate_cross_K2D <- function(spe, 
                                reference_cell_type, 
                                target_cell_types, 
                                radius, 
                                feature_colname = "Cell.Type") {
  
  if (is.null(spe[[feature_colname]])) stop(paste("No column called", feature_colname, "found in spe object"))
  
  if (is.null(spe[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spe")
    spe$Cell.ID <- paste("Cell", seq(ncol(spe)), sep = "_")
  }  
  
  
  ## Get expected cross K-function
  expected_cross_K <- pi * radius^2
  
  ## For reference_cell_type, check it is found in the spe object
  if (!(reference_cell_type %in% spe[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spe object"))
    result <- data.frame(observed_cross_K = NA,
                         expected_cross_K = expected_cross_K,
                         cross_K_ratio = NA)
    return(result)
  }
  
  ## Get rough dimensions of the window the points are in
  spe_coords <- data.frame(spatialCoords(spe))
  
  length <- round(max(spe_coords$Cell.X.Position) - min(spe_coords$Cell.X.Position))
  width  <- round(max(spe_coords$Cell.Y.Position) - min(spe_coords$Cell.Y.Position))
  
  ## Get volume of the window the cells are in
  volume <- length * width
  
  # Number of reference cell types is constant
  n_ref_cells <- sum(spe[[feature_colname]] == reference_cell_type)
  
  # Define result data frame
  result <- data.frame(reference = reference_cell_type, expected = expected_cross_K)
  
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood2D(spe,
                                                                  reference_cell_type,
                                                                  target_cell_types,
                                                                  radius,
                                                                  feature_colname,
                                                                  show_summary = FALSE,
                                                                  plot_image = FALSE)
  
  for (target_cell_type in target_cell_types) {
    
    n_ref_tar_interactions <- sum(cells_in_neighbourhood_df[[target_cell_type]])
    
    n_tar_cells <- sum(spe[[feature_colname]] == target_cell_type)
    
    ## Get observed cross K-function
    if (n_tar_cells == 0) {
      observed_cross_K <- NA
    }
    else {
      observed_cross_K <- (volume * n_ref_tar_interactions) / (n_ref_cells * n_tar_cells)  
    }
    result[[target_cell_type]] <- observed_cross_K
  }
  
  return(result)
}
calculate_cross_L_gradient2D <- function(spe, 
                                         reference_cell_type, 
                                         target_cell_types, 
                                         radii, 
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = 2 + length(target_cell_types)))
  colnames(result) <- c("reference", "expected", target_cell_types)
  
  for (i in seq(length(radii))) {
    cross_L_df <- calculate_cross_L2D(spe,
                                      reference_cell_type,
                                      target_cell_types,
                                      radii[i],
                                      feature_colname)
    
    result[i, ] <- cross_L_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig1 <- plot_cross_L_gradient2D(result)
    fig2 <- plot_cross_L_gradient_ratio2D(result)
    
    combined_fig <- plot_grid(fig1, fig2, nrow = 2)
    methods::show(combined_fig)
  }
  
  return(result)
}
calculate_cross_L2D <- function(spe, 
                                reference_cell_type, 
                                target_cell_types, 
                                radius, 
                                feature_colname = "Cell.Type") {
  
  result <- calculate_cross_K2D(spe = spe,
                                reference_cell_type = reference_cell_type,
                                target_cell_types = target_cell_types,
                                radius = radius,
                                feature_colname = feature_colname)
  
  result[ , c("expected", target_cell_types)] <- (result[ , c("expected", target_cell_types)] / (pi)) ^ (1/2)
  
  return(result)
}
calculate_entropy_background2D <- function(spe,
                                           cell_types_of_interest, 
                                           feature_colname = "Cell.Type") {
  
  # NULL case: entropy is undefined
  if (is.null(cell_types_of_interest)) return(NA)
  
  # One cell type case: entropy is 0
  if (is.character(cell_types_of_interest) && length(cell_types_of_interest) == 1) return(0)
  
  cell_proportions_data <- calculate_cell_proportions2D(spe, cell_types_of_interest, feature_colname, FALSE)
  
  # Calculate entropy of the entire image
  entropy <- -1 * sum(cell_proportions_data$proportion * log(cell_proportions_data$proportion, length(cell_proportions_data$proportion)))
  
  return(entropy) 
}

calculate_entropy_gradient2D <- function(spe,
                                         reference_cell_type,
                                         target_cell_types,
                                         radii,
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = 1))
  colnames(result) <- "entropy"
  
  for (i in seq(length(radii))) {
    entropy_df <- calculate_entropy2D(spe,
                                      reference_cell_type,
                                      target_cell_types,
                                      radii[i],
                                      feature_colname)
    
    if (is.null(entropy_df)) return(NULL)
    
    result[i, "entropy"] <- mean(entropy_df$entropy, na.rm = T)
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    expected_entropy <- calculate_entropy_background2D(spe, target_cell_types, feature_colname)
    fig <- plot_entropy_gradient2D(result, expected_entropy, reference_cell_type, target_cell_types)
    methods::show(fig)
  }
  
  return(result)
}

calculate_entropy_grid_metrics2D <- function(spe, 
                                             n_splits,
                                             cell_types_of_interest,
                                             feature_colname = "Cell.Type",
                                             plot_image = TRUE) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  ## Check cell_types_of_interest are found in the spe object
  unknown_cell_types <- setdiff(cell_types_of_interest, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  ## If cell types have been chosen, check they are found in the spe object
  unknown_cell_types <- setdiff(cell_types_of_interest, unique(spe[[feature_colname]]))
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Add grid metrics to spe
  spe <- get_spe_grid_metrics2D(spe, n_splits, feature_colname)
  
  # Get grid_prism_cell_matrix from spe
  grid_prism_cell_matrix <- spe@metadata$grid_metrics$grid_prism_cell_matrix
  
  ## Define data frame which contains all results
  n_grid_prisms <- n_splits^3
  result <- data.frame(row.names = seq(n_grid_prisms))
  
  for (cell_type in cell_types_of_interest) {
    result[[cell_type]] <- grid_prism_cell_matrix[[cell_type]]
  }
  result$total <- rowSums(result)
  
  ## Get data frame containing proportions for cell_types_of_interest
  df_props <- result[ , cell_types_of_interest] / result$total
  
  ## Use proportion data frame to get entropy
  calculate_entropy <- function(x) {
    entropy <- -1 * sum(x * ifelse(is.infinite(log(x, length(x))), 0, log(x, length(x))))
    return(entropy)
  }
  result$entropy <- apply(df_props, 1, calculate_entropy)
  
  # Add grid_prism_coordinates info to result
  result <- cbind(result, spe@metadata$grid_metrics$grid_prism_coordinates)
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_metrics_continuous2D(result, "entropy")
    methods::show(fig)
  }
  
  return(result)
}
calculate_entropy2D <- function(spe,
                                reference_cell_type,
                                target_cell_types,
                                radius,
                                feature_colname = "Cell.Type") {
  
  # Check target_cell_types
  if (!(is.character(target_cell_types) && length(target_cell_types) >= 2)) {
    stop("`target_cell_types` is not a character vector with at least 2 cell types.")
  }
  
  ## Users should ensure include the reference_cell_type as one of the target_cell_types
  cells_in_neighbourhood_proportion_df <- calculate_cells_in_neighbourhood_proportions2D(spe,
                                                                                         reference_cell_type,
                                                                                         target_cell_types,
                                                                                         radius,
                                                                                         feature_colname)
  
  if (is.null(cells_in_neighbourhood_proportion_df)) return(NULL)
  
  ## Get entropy for each row
  cells_in_neighbourhood_proportion_df$entropy <- apply(cells_in_neighbourhood_proportion_df[ , paste(target_cell_types, "_prop", sep = "")],
                                                        1,
                                                        function(x) -1 * sum(x * log(x, length(target_cell_types))))
  cells_in_neighbourhood_proportion_df$entropy <- ifelse(cells_in_neighbourhood_proportion_df$total > 0 & is.nan(cells_in_neighbourhood_proportion_df$entropy), 
                                                         0,
                                                         cells_in_neighbourhood_proportion_df$entropy)
  
  return(cells_in_neighbourhood_proportion_df)
}

calculate_minimum_distances_between_cell_types2D <- function(spe,
                                                             cell_types_of_interest = NULL,
                                                             feature_colname = "Cell.Type",
                                                             show_summary = TRUE,
                                                             plot_image = TRUE) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (ncol(spe) < 2) {
    stop("There must be at least two cells in spe.")
  }
  if (!(is.null(cell_types_of_interest) || is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(show_summary)) {
    stop("`show_summary` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  if (is.null(spe[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spe")
    spe$Cell.ID <- paste("Cell", seq(ncol(spe)), sep = "_")
  }  
  
  # De-factor feature column in spe object
  spe[[feature_colname]] <- as.character(spe[[feature_colname]])
  
  # Subset spe to only contain the cells of interest
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spe object
    unknown_cell_types <- setdiff(cell_types_of_interest, spe[[feature_colname]])
    if (length(unknown_cell_types) != 0) {
      warning(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                    paste(unknown_cell_types, collapse = ", ")))
    }
    
    spe <- spe[ , spe[[feature_colname]] %in% cell_types_of_interest]
  }
  # If cell_types_of_interest is NULL, use all cells in spe
  else {
    cell_types_of_interest <- unique(spe[[feature_colname]])
  }
  
  # Create a list containing the cell IDs of each cell type
  cell_type_ids <- list()
  for (cell_type in cell_types_of_interest) {
    cell_type_ids[[cell_type]] <- as.character(spe$Cell.ID[spe[[feature_colname]] == cell_type])
  }
  
  # Get spe coords
  spe_coords <- data.frame(spatialCoords(spe))
  
  # Get different possible cell type combinations
  # Each row represents a combination
  # If a row is [1 , 2], then we are comparing cell type 1 and cell type 2
  permu <- gtools::permutations(length(cell_types_of_interest), 2, repeats.allowed = TRUE)
  
  result <- data.frame()
  
  for (i in seq(nrow(permu))) {
    cell_type1 <- cell_types_of_interest[permu[i, 1]]
    cell_type2 <- cell_types_of_interest[permu[i, 2]]
    
    # Don't have one of the cells
    if (sum(spe[[feature_colname]] == cell_type1) == 0 || sum(spe[[feature_colname]] == cell_type2) == 0) {
      result <- rbind(result, data.frame(ref_cell_id = NA, ref_cell_type = cell_type1, nearest_cell_id = NA, nearest_cell_type = cell_type2, distance = NA))
      next
    }
    
    # Get x, y, z coords for all cells of cell_type1 and cell_type2
    cell_type1_coords <- spe_coords[spe[[feature_colname]] == cell_type1, ]
    cell_type2_coords <- spe_coords[spe[[feature_colname]] == cell_type2, ]
    
    # Find all of closest points
    # For each cell of cell_type1, find the closest cell of cell_type2
    if (cell_type1 != cell_type2) {
      nearest_neighbours <- RANN::nn2(data = cell_type2_coords, 
                                      query = cell_type1_coords, 
                                      k = 1)  
    }
    # If we are comparing the same cell_type, and there is only one of this cell type, move on
    else if (nrow(cell_type1_coords) == 1) {
      warning("There is only 1 '", cell_type1, "' cell in your data. It has no nearest neighbour of the same cell type.", sep = "")
      result <- rbind(result, data.frame(ref_cell_id = NA, ref_cell_type = cell_type1, nearest_cell_id = NA, nearest_cell_type = cell_type2, distance = NA))
      next
    }
    # If we are comparing the same cell_type, use the second closest neighbour
    else {
      nearest_neighbours <- RANN::nn2(data = cell_type2_coords, 
                                      query = cell_type1_coords, 
                                      k = 2)
      nearest_neighbours[['nn.idx']] <- nearest_neighbours[['nn.idx']][ , 2]
      nearest_neighbours[['nn.dists']] <- nearest_neighbours[['nn.dists']][ , 2]
    }
    
    # Create the data frame containing the chosen cells and their ids, as well as the nearest cell to them and their ids, and the distance between
    
    df <- data.frame(
      ref_cell_id = cell_type_ids[[cell_type1]],
      ref_cell_type = cell_type1,
      nearest_cell_id = cell_type_ids[[cell_type2]][c(nearest_neighbours$nn.idx)],
      nearest_cell_type = cell_type2,
      distance = nearest_neighbours$nn.dists
    )
    result <- rbind(result, df)
  }
  
  result$pair <- paste(result$ref_cell_type, result$nearest_cell_type,sep = "/")
  
  # Print summary
  if (show_summary) {
    print(summarise_distances_between_cell_types2D(result))  
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_distances_between_cell_types_violin2D(result)
    methods::show(fig)
  }
  
  return(result)
}
calculate_minimum_distances_to_clusters2D <- function(spe, 
                                                      cell_types_inside_cluster, 
                                                      cell_types_outside_cluster, 
                                                      cluster_colname, 
                                                      feature_colname = "Cell.Type", 
                                                      plot_image = T) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (!is.character(cell_types_inside_cluster)) {
    stop("`cell_types_inside_cluster` is not a character vector.")
  }
  if (!is.character(cell_types_outside_cluster)) {
    stop("`cell_types_outside_cluster` is not a character vector.")
  }
  if (!(is.numeric(radius) && length(radius) == 1 && radius > 0)) {
    stop("`radius` is not a positive numeric.")
  }
  if (!is.character(cluster_colname)) {
    stop("`cluster_colname` is not a character. This should be 'alpha_hull_cluster', 'dbscan_cluster', or 'grid_based_cluster', depending on the chosen method.")
  }
  if (is.null(spe[[cluster_colname]])) {
    stop(paste("No column called", cluster_colname, "found in spe object."))
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  ## Add Cell.ID column
  if (is.null(spe[["Cell.ID"]])) {
    warning("Temporarily adding Cell.Id column to your spe")
    spe$Cell.ID <- paste("Cell", seq(ncol(spe)), sep = "_")
  }
  
  ## For each cell type outside clusters, get their set of coords. These exclude cell types in clusters
  spe_coords <- spatialCoords(spe)
  
  # Cells outside cluster have a cluster number of 0 (i.e. they are not in a cluster)
  spe_outside_cluster <- spe[ , spe[[cluster_colname]] == 0]
  
  cell_types_outside_cluster_coords <- list()
  for (cell_type in cell_types_outside_cluster) {
    cell_types_outside_cluster_coords[[cell_type]] <- spatialCoords(spe_outside_cluster)[spe_outside_cluster[[feature_colname]] == cell_type, ]
  }
  
  ## For each cluster, determine the minimum distance of each outside_cell_type  
  result <- vector()
  
  # Get number of clusters
  n_clusters <- max(spe[[cluster_colname]])
  
  for (i in seq(n_clusters)) {
    cluster_coords <- spe_coords[spe[[cluster_colname]] == i & spe[[feature_colname]] %in% cell_types_inside_cluster, ]
    cluster_cell_types <- spe[["Cell.Type"]][spe[[cluster_colname]] == i & spe[[feature_colname]] %in% cell_types_inside_cluster]
    cluster_cell_ids <- spe[["Cell.ID"]][spe[[cluster_colname]] == i & spe[[feature_colname]] %in% cell_types_inside_cluster]
    
    for (outside_cell_type in cell_types_outside_cluster) {
      curr_cell_type_coords <- cell_types_outside_cluster_coords[[outside_cell_type]]
      
      all_closest <- RANN::nn2(data = cluster_coords, 
                               query = curr_cell_type_coords, 
                               k = 1) 
      
      local_dist_mins <- data.frame(
        cluster_number = i,
        outside_cell_id = as.character(spe_outside_cluster$Cell.ID[spe_outside_cluster[["Cell.Type"]] == outside_cell_type]),
        outside_cell_type = outside_cell_type,
        inside_cell_id = cluster_cell_ids[c(all_closest$nn.idx)],
        inside_cell_type = cluster_cell_types[c(all_closest$nn.idx)],
        distance = all_closest$nn.dists
      )
      ## Remove any 0 distance rows
      local_dist_mins <- local_dist_mins[local_dist_mins$distance != 0, ]
      result <- rbind(result, local_dist_mins)
    }
    
    
    ## Plot
    if (plot_image) {
      
      cluster_number_labs <- paste("cluster_", seq(n_clusters), sep = "")
      names(cluster_number_labs) <- seq(n_clusters)
      
      fig <- ggplot(result, aes(x = outside_cell_type, y = distance, fill = outside_cell_type)) + 
        geom_violin() +
        facet_grid(cluster_number~., scales="free_x", labeller = labeller(cluster_number = cluster_number_labs)) +
        theme_bw() +
        theme(axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), legend.position = "none") +
        labs(title="Minimum cell distances to clusters", x = "Cell type", y = "Distance") +
        stat_summary(fun.data = "mean_sdl", fun.args = list(mult= 1), colour = "red")
      
      methods::show(fig)
    }
    
  }
  return(result)
}
calculate_mixing_scores_gradient2D <- function(spe, 
                                               reference_cell_type, 
                                               target_cell_type, 
                                               radii, 
                                               feature_colname = "Cell.Type",
                                               plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = 8))
  colnames(result) <- c("ref_cell_type", 
                        "tar_cell_type", 
                        "n_ref_cells",
                        "n_tar_cells", 
                        "n_ref_tar_interactions",
                        "n_ref_ref_interactions", 
                        "mixing_score", 
                        "normalised_mixing_score")
  
  for (i in seq(length(radii))) {
    mixing_scores <- calculate_mixing_scores2D(spe,
                                               reference_cell_type,
                                               target_cell_type,
                                               radii[i],
                                               feature_colname)
    
    result[i, ] <- mixing_scores
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig1 <- plot_mixing_scores_gradient2D(result, "NMS")
    fig2 <- plot_mixing_scores_gradient2D(result, "MS")
    combined_fig <- plot_grid(fig1, fig2, nrow = 2)
    methods::show(combined_fig)
  }
  
  return(result)
}

calculate_mixing_scores2D <- function(spe, 
                                      reference_cell_types, 
                                      target_cell_types, 
                                      radius, 
                                      feature_colname = "Cell.Type") {
  
  # Define result
  result <- data.frame()
  
  for (reference_cell_type in reference_cell_types) {
    
    for (target_cell_type in target_cell_types) {
      
      # No point getting mixing scores if comparing the same cell type
      if (reference_cell_type == target_cell_type) {
        next
      }
      
      # Get number of reference cells and target cells
      n_ref <- sum(spe[[feature_colname]] == reference_cell_type)
      n_tar <- sum(spe[[feature_colname]] == target_cell_type)
      
      
      # Can't get mixing scores if there are 0 or 1 reference cells
      if (n_ref == 0 || n_ref == 1) {
        result <-  rbind(result, 
                         c(reference_cell_type, 
                           target_cell_type, 
                           n_ref, 
                           n_tar, 
                           0, 
                           0, 
                           NA, 
                           NA))
      }
      
      
      ## Get cells in neighbourhood df
      cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood2D(spe,
                                                                      reference_cell_type,
                                                                      c(reference_cell_type, target_cell_type),
                                                                      radius,
                                                                      feature_colname,
                                                                      FALSE,
                                                                      FALSE)
      
      # Get number of ref-ref interactions
      # Halve it to avoid counting each ref-ref interaction twice
      n_ref_ref_interactions <- 0.5 * sum(cells_in_neighbourhood_df[[reference_cell_type]]) 
      
      # Get number of ref-tar interactions
      n_ref_tar_interactions <- sum(cells_in_neighbourhood_df[[target_cell_type]]) 
      
      
      # Can't get mixing scores if there are no target cells
      if (n_tar == 0) {
        
        result <-  rbind(result, 
                         c(reference_cell_type, 
                           target_cell_type, 
                           n_ref, 
                           0, 
                           0, 
                           n_ref_ref_interactions, 
                           NA, 
                           NA))
      }
      
      # Generic case: We have reference cells and target cells
      else {
        
        if (n_ref_ref_interactions != 0) {
          mixing_score <- n_ref_tar_interactions / n_ref_ref_interactions
          normalised_mixing_score <- 0.5 * mixing_score * n_ref / n_tar
        }
        else {
          mixing_score <- 0
          normalised_mixing_score <- 0
          methods::show(paste("There are no reference to reference interactions for", target_cell_type, "in the specified radius, cannot calculate mixing score"))
        }
        
        result <-  rbind(result, 
                         c(reference_cell_type, 
                           target_cell_type, 
                           n_ref, 
                           n_tar, 
                           n_ref_tar_interactions, 
                           n_ref_ref_interactions, 
                           mixing_score, 
                           normalised_mixing_score))
      }
    }
  }
  
  # Required column names of our output data frame
  colnames(result) <- c("ref_cell_type", 
                        "tar_cell_type", 
                        "n_ref_cells",
                        "n_tar_cells", 
                        "n_ref_tar_interactions",
                        "n_ref_ref_interactions", 
                        "mixing_score", 
                        "normalised_mixing_score")
  
  # Turn numeric data into numeric type
  result[ , 3:8] <- apply(result[ , 3:8], 2, as.numeric)
  
  return(result)
}
calculate_pairwise_distances_between_cell_types2D <- function(spe,
                                                              cell_types_of_interest = NULL,
                                                              feature_colname = "Cell.Type",
                                                              show_summary = TRUE,
                                                              plot_image = TRUE) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (ncol(spe) < 2) {
    stop("There must be at least two cells in spe.")
  }
  if (!(is.null(cell_types_of_interest) && is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(show_summary)) {
    stop("`show_summary` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  if (is.null(spe[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spe")
    spe$Cell.ID <- paste("Cell", seq(ncol(spe)), sep = "_")
  }
  
  # De-factor feature column in spe object
  spe[[feature_colname]] <- as.character(spe[[feature_colname]])
  
  # Subset spe to only contain the cells of interest
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spe object
    unknown_cell_types <- setdiff(cell_types_of_interest, spe[[feature_colname]])
    if (length(unknown_cell_types) != 0) {
      warning(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                    paste(unknown_cell_types, collapse = ", ")))
    }
    
    spe <- spe[ , spe[[feature_colname]] %in% cell_types_of_interest]
  }
  # If cell_types_of_interest is NULL, use all cells in spe
  else {
    cell_types_of_interest <- unique(spe[[feature_colname]])
  }
  
  # Create a list containing the cell IDs of each cell type
  cell_type_ids <- list()
  for (cell_type in cell_types_of_interest) {
    cell_type_ids[[cell_type]] <- as.character(spe$Cell.ID[spe[[feature_colname]] == cell_type])
  }
  
  # Calculate cell to cell distances
  distance_matrix <- -1 * apcluster::negDistMat(spatialCoords(spe))
  rownames(distance_matrix) <- spe$Cell.ID
  colnames(distance_matrix) <- spe$Cell.ID
  
  result <- data.frame()
  
  for (i in seq(length(cell_types_of_interest))) {
    
    for (j in i:length(cell_types_of_interest)) {
      
      # Get current cell types and cell ids
      cell_type1 <- names(cell_type_ids)[i]
      cell_type2 <- names(cell_type_ids)[j]
      
      cell_type1_ids <- cell_type_ids[[cell_type1]]
      cell_type2_ids <- cell_type_ids[[cell_type2]]
      
      ## Don't have a cell type, or the same cell type with only one cell
      if (length(cell_type1_ids) == 0 || length(cell_type2_ids) == 0) {
        result <- rbind(result, data.frame(Var1 = NA, Var2 = NA, value = NA, cell_type1 = cell_type1, cell_type2 = cell_type2, pair = paste(cell_type1, cell_type2, sep="/")))
        next
      }
      
      ## Same cell type only one cell
      if (cell_type1 == cell_type2 && length(cell_type1_ids) == 1) {
        warning("There is only 1 '", cell_type1, "' cell in your data. It has no pair of the same cell type.", sep = "")
        result <- rbind(result, data.frame(Var1 = NA, Var2 = NA, value = NA, cell_type1 = cell_type1, cell_type2 = cell_type2, pair = paste(cell_type1, cell_type2, sep="/")))
        next
      }
      
      # Subset distance_matrix for current cell types
      distance_matrix_subset <- distance_matrix[rownames(distance_matrix) %in% cell_type1_ids, 
                                                colnames(distance_matrix) %in% cell_type2_ids]
      
      ## Different cell types, each only has one cell
      if (length(cell_type1_ids) == 1 && length(cell_type2_ids) == 1) {
        distance_matrix_subset <- as.matrix(distance_matrix_subset)
        rownames(distance_matrix_subset) <- cell_type1_ids
        colnames(distance_matrix_subset) <- cell_type2_ids
      }    
      ## Different cell types, only one cell of cell_type1
      else if (length(cell_type1_ids) == 1) {
        distance_matrix_subset <- as.matrix(distance_matrix_subset)
        colnames(distance_matrix_subset) <- cell_type1_ids
      }
      ## Different cell types, only one cell of cell_type2
      else if (length(cell_type2_ids) == 1) {
        distance_matrix_subset <- as.matrix(distance_matrix_subset)
        colnames(distance_matrix_subset) <- cell_type2_ids
      }
      ## Same cell type, only need part of the matrix (make irrelevant part of matrix equal to NA)
      if (cell_type1 == cell_type2) distance_matrix_subset[upper.tri(distance_matrix_subset, diag = TRUE)] <- NA
      
      # Convert distance_matrix_subset to a data frame
      df <- reshape2::melt(distance_matrix_subset, na.rm = TRUE)
      df$cell_type1 <- cell_type1
      df$cell_type2 <- cell_type2
      df$pair <- paste(cell_type1, cell_type2, sep="/")
      
      result <- rbind(result, df)
    }
  }
  
  # Rearrange columns 
  colnames(result)[c(1, 2, 3)] <- c("cell_type1_id", "cell_type2_id", "distance")
  result <- result[ , c("cell_type1_id", "cell_type1", "cell_type2_id", "cell_type2", "distance", "pair")]
  
  # Print summary
  if (show_summary) {
    print(summarise_distances_between_cell_types2D(result))  
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_distances_between_cell_types_violin2D(result)
    methods::show(fig)
  }
  
  return(result)
}
calculate_prevalence_gradient_AUC2D <- function(prevalence_gradient_df) {
  
  return(sum(prevalence_gradient_df$prevalence) * 0.01)
}
calculate_prevalence_gradient2D <- function(grid_metrics,
                                            metric_colname,
                                            show_AUC = T,
                                            plot_image = T) {
  
  ## Check input parameters
  if (!(is.character(metric_colname))) {
    stop("`metric_colname` is not a character. This should be 'proportion' or 'entropy', depending on the chosen method.")
  }
  if (is.null(grid_metrics[[metric_colname]])) {
    stop("`metric_colname` is not a column in `grid_metrics`.")
  }
  if (!is.logical(show_AUC)) {
    stop("`show_AUC` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Thresholds range from 0 to 1
  thresholds <- seq(0.01, 1, 0.01)
  
  # Define result
  result <- data.frame(threshold = thresholds)
  
  # Get prevalences for each threshold
  result$prevalence <- sapply(thresholds, function(threshold) { 
    calculate_prevalence2D(grid_metrics, metric_colname, threshold) 
  })
  
  # Show AUC of prevalence gradient graph
  if (show_AUC) {
    print(paste("AUC:", round(calculate_prevalence_gradient_AUC2D(result), 2)))
  }
  
  # Plot
  if (plot_image) {
    fig <- ggplot(result, aes(threshold, prevalence)) +
      geom_line() +
      theme_bw() +
      labs(x = "Threshold",
           y = "Prevalence",
           title = paste("Prevalence vs Threshold (", metric_colname, ")", sep = "")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      ylim(0, 100)
    methods::show(fig)
  }
  
  return(result)
}
calculate_prevalence2D <- function(grid_metrics,
                                   metric_colname,
                                   threshold,
                                   above = TRUE) {
  
  ## Check input parameters
  if (!(is.character(metric_colname))) {
    stop("`metric_colname` is not a character. This should be 'proportion' or 'entropy', depending on the chosen method.")
  }
  if (is.null(grid_metrics[[metric_colname]])) {
    stop("`metric_colname` is not a column in `grid_metrics`.")
  }
  if (!(is.numeric(threshold) && length(threshold) == 1 && threshold >= 0 && threshold <= 1)) {
    stop("`threshold` is not a numeric between 0 and 1.")
  }
  if (!is.logical(above)) {
    stop("`above` is not a logical (TRUE or FALSE).")
  }
  
  
  ## Exclude rows with NA values
  grid_metrics <- grid_metrics[!is.na(grid_metrics[[metric_colname]]), ]
  
  if (above) {
    prevalence <- sum(grid_metrics[[metric_colname]] >= threshold) / nrow(grid_metrics) * 100
  }
  else {
    prevalence <- sum(grid_metrics[[metric_colname]] < threshold) / nrow(grid_metrics) * 100    
  }
  
  return(prevalence)
}
calculate_spatial_autocorrelation2D <- function(grid_metrics,
                                                metric_colname,
                                                weight_method = 0.1) {
  
  ## Check input parameters
  if (!(is.character(metric_colname))) {
    stop("`metric_colname` is not a character. This should be 'proportion' or 'entropy', depending on the chosen method.")
  }
  if (is.null(grid_metrics[[metric_colname]])) {
    stop("`metric_colname` is not a column in `grid_metrics`.")
  }
  if (!((is.numeric(weight_method) && length(weight_method) == 1 && weight_method > 0 && weight_method < 1) ||
        (is.character(weight_method) && weight_method %in% c("IDW", "rook", "queen")))) {
    stop("`weight_method` is not a numeric between 0 and 1 or either 'IDW', 'rook' or 'queen'.")
  }
  
  ## Get number of grid prisms
  n_grid_prisms <- nrow(grid_metrics)
  
  ## Get splitting number (should be the SQUARE root of n_grid_prisms)
  n_splits <- (n_grid_prisms)^(1/2)
  
  ## Find the coordinates of each grid prism
  x <- ((seq(n_grid_prisms) - 1) %% n_splits)
  y <- (floor(((seq(n_grid_prisms) - 1) %% (n_splits)^2) / n_splits))
  grid_prism_coords <- data.frame(x = x, y = y)
  
  ## Subset for non NA rows
  grid_prism_coords <- grid_prism_coords[!is.na(grid_metrics[[metric_colname]]), ]
  grid_metrics <- grid_metrics[!is.na(grid_metrics[[metric_colname]]), ]
  
  weight_matrix <- -1 * apcluster::negDistMat(grid_prism_coords)
  ## Use the inverse distance between two points as the weight (IDW is 'inverse distance weighting')
  if (weight_method == "IDW") {
    weight_matrix <- 1 / weight_matrix
  }
  ## Use rook method: adjacent points get a weight of 1, otherwise, weight of 0
  ## Adjacent points are within 1 unit apart. e.g. (0, 0, 0) vs (0, 0, 1)
  else if (weight_method == "rook") {
    weight_matrix <- ifelse(weight_matrix > 1, 0, 1)  
  }
  ## Use queen method: adjacent points get a weight of 1, otherwise, weight of 0
  ## Adjacent points are within sqrt(2) unit apart. e.g. (0, 0, 0) vs (0, 0, 1)
  else if (weight_method == "queen") {
    weight_matrix <- ifelse(weight_matrix > sqrt(2), 0, 1)  
  }
  ## If a number (x) between 0 and 1 is supplied, set a threshold to be x quantile value of c(weight_matrix)
  ## Grid prisms within this specified threshold have a weight of 1, otherwise, weight of 0
  else if (as.numeric(weight_method) && 0 < weight_method && weight_method < 1) {
    threshold <- quantile(c(weight_matrix), weight_method)
    weight_matrix <- ifelse(weight_matrix > threshold, 0, 1)
  }
  
  ## Points along the diagonal are comparing the same point so its weight is zero
  diag(weight_matrix) <- 0
  
  n <- nrow(grid_metrics)
  
  # Center the data
  data_centered <- grid_metrics[[metric_colname]] - mean(grid_metrics[[metric_colname]])
  
  # Calculate numerator using matrix multiplication
  numerator <- sum(data_centered * (weight_matrix %*% data_centered))
  
  # Calculate denominator
  denominator <- sum(data_centered^2) * sum(weight_matrix)
  
  # Moran's I
  I <- (n * numerator) / denominator
  
  return(I)
}



get_spe_grid_metrics2D <- function(spe, 
                                   n_splits, 
                                   feature_colname = "Cell.Type") {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  
  spe_coords <- spatialCoords(spe)
  
  ## Get dimensions of the window
  min_x <- min(spe_coords[ , "Cell.X.Position"])
  min_y <- min(spe_coords[ , "Cell.Y.Position"])
  
  max_x <- max(spe_coords[ , "Cell.X.Position"])
  max_y <- max(spe_coords[ , "Cell.Y.Position"])
  
  length <- round(max_x - min_x)
  width  <- round(max_y - min_y)
  
  ## Get distance of row, col and lay
  d_row <- length / n_splits
  d_col <- width / n_splits
  
  # Shift spe_coords so they begin at the origin
  spe_coords[, "Cell.X.Position"] <- spe_coords[, "Cell.X.Position"] - min_x
  spe_coords[, "Cell.Y.Position"] <- spe_coords[, "Cell.Y.Position"] - min_y
  
  ## Figure out which 'grid prism number' each cell is inside
  spe$grid_prism_num <- floor(spe_coords[ , "Cell.X.Position"] / d_row) +
    floor(spe_coords[ , "Cell.Y.Position"] / d_col) * n_splits
  
  ## Determine the cell types found in each grid prism
  n_grid_prisms <- n_splits^2
  grid_prism_cell_matrix <- as.data.frame.matrix(table(spe[[feature_colname]], factor(spe$grid_prism_num, levels = seq(n_grid_prisms))))
  grid_prism_cell_matrix <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                       t(grid_prism_cell_matrix))
  
  ## Determine centre coordinates of each grid prism
  grid_prism_coordinates <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                       x_coord = ((seq(n_grid_prisms) - 1) %% n_splits + 0.5) * d_row + round(min_x),
                                       y_coord = (floor(((seq(n_grid_prisms) - 1) %% (n_splits)^2) / n_splits) + 0.5) * d_col + round(min_y))
  
  spe@metadata[["grid_metrics"]] <- list("grid_prism_cell_matrix" = grid_prism_cell_matrix,
                                         "grid_prism_coordinates" = grid_prism_coordinates)
  
  return(spe)
}

plot_cells_in_neighbourhood_gradient2D <- function(cells_in_neighbourhood_gradient_df, reference_cell_type = NULL) {
  
  plot_result <- reshape2::melt(cells_in_neighbourhood_gradient_df, "radius")
  
  fig <- ggplot(plot_result, aes(radius, value, color = variable)) + 
    geom_line() + 
    labs(title = "Average cells in neighbourhood gradient", x = "Radius", y = "Average cells in neighbourhood") + 
    scale_color_discrete(name = "Cell type") +
    theme_bw()
  
  if (!is.null(reference_cell_type)) {
    fig <- fig + labs(subtitle = paste("Reference: ", reference_cell_type, ", Target: ", paste(colnames(cells_in_neighbourhood_gradient_df)[seq(ncol(cells_in_neighbourhood_gradient_df) - 1)], collapse = ", "), sep = ""))
  }
  
  return(fig)
}
plot_cells_in_neighbourhood_proportions_gradient2D <- function(cells_in_neighbourhood_proportions_gradient_df, reference_cell_type = NULL) {
  
  plot_result <- reshape2::melt(cells_in_neighbourhood_proportions_gradient_df, id.vars = c("radius"))
  fig <- ggplot(plot_result, aes(radius, value, color = variable)) +
    geom_point() +
    geom_line() +
    labs(title = "Average cells in neighbourhood proportions gradient", x = "Radius", y = "Cell proportion", color = "Cell type") +
    theme_bw() +
    ylim(0, 1)
  
  if (!is.null(reference_cell_type)) {
    fig <- fig + labs(subtitle = paste("Reference: ", reference_cell_type, ", Target: ", paste(colnames(cells_in_neighbourhood_proportions_gradient_df)[seq(ncol(cells_in_neighbourhood_proportions_gradient_df) - 1)], collapse = ", "), sep = ""))
  }
  
  return(fig)
}
## For scales parameter, use "free_x" or "free". "free_y" looks silly
plot_cells_in_neighbourhood_violin2D <- function(cells_in_neighbourhood_df, reference_cell_type, scales = "free_x") {
  
  ## Target cell types will be all the columns except the first column
  target_cell_types <- colnames(cells_in_neighbourhood_df)[c(-1)]
  
  df <- reshape2::melt(cells_in_neighbourhood_df, measure.vars = target_cell_types)
  colnames(df) <- c("ref_cell_id", "tar_cell_type", "count")
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  tar_cell_type <- count <- NULL
  
  fig <- ggplot(df, aes(x = tar_cell_type, y = count)) + 
    geom_violin() +
    facet_wrap(~tar_cell_type, scales=scales, strip.position="bottom") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(title=paste("Cells in neighbourhood of", reference_cell_type, "cells"), x = "Target cell type", y = "Number of cells") +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult= 1), colour = "red")
  
  message("Plots show mean ± sd")
  
  return(fig)
}

plot_co_occurrence_gradient2D <- function(co_occurrence_gradient_df) {
  
  target_cell_types <- colnames(co_occurrence_gradient_df)
  target_cell_types <- target_cell_types[!target_cell_types %in% c("reference", "radius")]
  
  co_occurrence_gradient_df$expected <- 1
  
  plot_result <- reshape2::melt(co_occurrence_gradient_df, "radius", c(target_cell_types, "expected"))
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Co-occurrence gradient", x = "Radius", y = "Co-occurrence value") +
    scale_colour_discrete(name = "") +
    theme_bw()
  
  return(fig) 
}
plot_cross_G_gradient2D <- function(cross_G_gradient_df, reference_cell_type = NULL, target_cell_type = NULL) {
  
  plot_result <- reshape2::melt(cross_G_gradient_df, "radius", c("observed_cross_G", "expected_cross_G"))
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Cross G-function gradient", x = "Radius", y = "Cross G-function value") +
    scale_colour_discrete(name = "", labels = c("Observed cross G", "Expected CSR cross G")) +
    theme_bw()
  
  if (!is.null(reference_cell_type) && !is.null(target_cell_type)) {
    fig <- fig + labs(subtitle = paste("Reference: ", reference_cell_type, ", Target: ", target_cell_type, sep = ""))
  }
  
  return(fig) 
}
plot_cross_K_gradient_ratio2D <- function(cross_K_gradient_df) {
  
  target_cell_types <- colnames(cross_K_gradient_df)[!colnames(cross_K_gradient_df) %in% c("reference", "expected", "radius")]
  
  # Normalize columns by 'expected'
  for (target_cell_type in target_cell_types) {
    cross_K_gradient_df[[target_cell_type]] <- cross_K_gradient_df[[target_cell_type]] / cross_K_gradient_df$expected
  }
  cross_K_gradient_df$expected <- 1
  
  plot_result <- reshape2::melt(cross_K_gradient_df, "radius", c(target_cell_types, "expected"))
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Cross K-function gradient ratio", x = "Radius", y = "Cross K-function ratio") +
    scale_colour_discrete(name = "") +
    theme_bw()
  
  return(fig) 
}
plot_cross_K_gradient2D <- function(cross_K_gradient_df) {
  
  target_cell_types <- colnames(cross_K_gradient_df)[!colnames(cross_K_gradient_df) %in% c("reference", "expected", "radius")]
  
  plot_result <- reshape2::melt(cross_K_gradient_df, "radius", c(target_cell_types, "expected"))
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Cross K-function gradient", x = "Radius", y = "Cross K-function value") +
    scale_colour_discrete(name = "") +
    theme_bw()
  
  return(fig) 
}
plot_cross_L_gradient_ratio2D <- function(cross_L_gradient_df) {
  
  target_cell_types <- colnames(cross_L_gradient_df)[!colnames(cross_L_gradient_df) %in% c("reference", "expected", "radius")]
  
  # Normalize columns by 'expected'
  for (target_cell_type in target_cell_types) {
    cross_L_gradient_df[[target_cell_type]] <- cross_L_gradient_df[[target_cell_type]] / cross_L_gradient_df$expected
  }
  cross_L_gradient_df$expected <- 1
  
  plot_result <- reshape2::melt(cross_L_gradient_df, "radius", c(target_cell_types, "expected"))
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Cross L-function gradient ratio", x = "Radius", y = "Cross L-function ratio") +
    scale_colour_discrete(name = "") +
    theme_bw()
  
  return(fig) 
}
plot_cross_L_gradient2D <- function(cross_L_gradient_df) {
  
  target_cell_types <- colnames(cross_L_gradient_df)[!colnames(cross_L_gradient_df) %in% c("reference", "expected", "radius")]
  
  plot_result <- reshape2::melt(cross_L_gradient_df, "radius", c(target_cell_types, "expected"))
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Cross L-function gradient", x = "Radius", y = "Cross L-function value") +
    scale_colour_discrete(name = "") +
    theme_bw()
  
  return(fig) 
}
## For scales parameter, use "free_x" or "free". "free_y" looks silly
plot_distances_between_cell_types_violin2D <- function(distances_df, scales = "free_x") {
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  pair <- distance <- NULL
  
  fig <- ggplot(distances_df, aes(x = pair, y = distance)) + 
    geom_violin() +
    facet_wrap(~pair, scales=scales, strip.position="bottom") +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
    labs(title="Cell distances", x = "Reference/Target pair", y = "Distance") +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult= 1), colour = "red")
  
  message("Plots show mean ± sd")
  
  return(fig)
}
plot_entropy_gradient2D <- function(entropy_gradient_df, expected_entropy = NULL, reference_cell_type = NULL, target_cell_types = NULL) {
  
  plot_result <- entropy_gradient_df
  
  if (!is.null(expected_entropy)) {
    if (!is.numeric(expected_entropy) || length(expected_entropy) != 1) stop("Please enter a single number for expected_entropy")
    plot_result$expected_entropy <- expected_entropy
    plot_result <- reshape2::melt(plot_result, "radius", c("entropy", "expected_entropy"))
    labels <- c("Observed entropy", "Expected CSR entropy")
  }
  else {
    plot_result <- reshape2::melt(plot_result, "radius", c("entropy"))
    labels <- c("Observed entropy")
  }
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Average entropy gradient", x = "Radius", y = "Entropy") +
    scale_colour_discrete(name = "", labels = labels) +
    theme_bw()
  
  if (!is.null(reference_cell_type) && !is.null(target_cell_types)) {
    fig <- fig + labs(subtitle = paste("Reference: ", reference_cell_type, ", Target: ", target_cell_types, sep = ""))
  }
  
  return(fig)
}

plot_mixing_scores_gradient2D <- function(mixing_scores_gradient_df, metric = "MS") {
  
  if (!metric %in% c("MS", "NMS")) {
    stop("'metric' should be 'MS' or 'NMS', for mixing score and normalised mixing score respectively.")
  }
  
  if (metric == "NMS") {
    plot_result <- mixing_scores_gradient_df
    plot_result$expected_normalised_mixing_score <- 1
    plot_result <- reshape2::melt(plot_result, "radius", c("normalised_mixing_score", "expected_normalised_mixing_score"))
    
    fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
      geom_line() +
      labs(title = "Normalised mixing score (NMS) gradient", 
           subtitle = paste("Reference: ", mixing_scores_gradient_df$ref_cell_type[1], ", Target: ", mixing_scores_gradient_df$tar_cell_type[1], sep = ""), 
           x = "Radius", y = "NMS") +
      scale_colour_discrete(name = "", labels = c("Observed NMS", "Expected CSR NMS")) +
      theme_bw() 
  }
  else if (metric == "MS") {
    plot_result <- mixing_scores_gradient_df
    n_tar_cells <- plot_result$n_tar_cells[1]
    n_ref_cells <- plot_result$n_ref_cells[1]
    plot_result$expected_mixing_score <- n_tar_cells * n_ref_cells / ((n_ref_cells - 1) * n_ref_cells / 2)
    plot_result <- reshape2::melt(plot_result, "radius", c("mixing_score", "expected_mixing_score"))
    
    fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
      geom_line() +
      labs(title = "Mixing score (MS) gradient", 
           subtitle = paste("Reference: ", mixing_scores_gradient_df$ref_cell_type[1], ", Target: ", mixing_scores_gradient_df$tar_cell_type[1], sep = ""), 
           x = "Radius", y = "MS") +
      scale_colour_discrete(name = "", labels = c("Observed MS", "Expected CSR MS  ")) +
      theme_bw()  
  }
  return(fig)
}

summarise_cells_in_neighbourhood2D <- function(cells_in_neighbourhood_df) {
  
  ## Target cell types will be all the columns except the first column
  target_cell_types <- colnames(cells_in_neighbourhood_df)[c(-1)]
  
  ## Set up data frame for summarised_results list
  df <- data.frame(row.names = c("mean", "min", "max", "median", "st_dev"))
  
  for (target_cell_type in target_cell_types) {
    
    ## Get statistical measures for each target cell type
    target_cell_type_values <- cells_in_neighbourhood_df[[target_cell_type]]
    df[[target_cell_type]] <- c(mean(target_cell_type_values),
                                min(target_cell_type_values),
                                max(target_cell_type_values),
                                median(target_cell_type_values),
                                sd(target_cell_type_values))
    
  }
  
  return(data.frame(t(df)))
}
summarise_distances_between_cell_types2D <- function(distances_df) {
  
  pair <- distance <- NULL
  
  # summarise the results
  distances_df_summarised <- distances_df %>% 
    dplyr::group_by(pair) %>%
    dplyr::summarise(mean(distance), 
                     min(distance), 
                     max(distance),
                     stats::median(distance), 
                     stats::sd(distance))
  
  distances_df_summarised <- data.frame(distances_df_summarised)
  
  colnames(distances_df_summarised) <- c("pair", 
                                         "mean", 
                                         "min", 
                                         "max", 
                                         "median", 
                                         "std_dev")
  
  for (i in seq(nrow(distances_df_summarised))) {
    # Get cell_types for each pair
    cell_types <- strsplit(distances_df_summarised[i,"pair"], "/")[[1]]
    
    distances_df_summarised[i, "reference"] <- cell_types[1]
    distances_df_summarised[i, "target"] <- cell_types[2]
  }
  
  return(distances_df_summarised)
}

### Slicing function -------------------------------------------------

# Function to get slices from spe
get_spe_slices_list <- function(spe, slice_bottom_z_coords, slice_top_z_coords) {
  
  spe_slices_list <- list()
  
  if (length(slice_bottom_z_coords) != length(slice_top_z_coords)) stop("Lengths of slice_bottom_z_coords and slice_top_z_coords should be equal.")
  n_slices <- length(slice_bottom_z_coords)
  
  for (i in seq(n_slices)) {
    bottom_z_coord <- slice_bottom_z_coords[i]
    top_z_coord <- slice_top_z_coords[i]
    spe_z_coords <- spatialCoords(spe)[ , "Cell.Z.Position"]
    spe_slice <- spe[, bottom_z_coord < spe_z_coords & spe_z_coords < top_z_coord]
    spatialCoords(spe_slice) <- spatialCoords(spe_slice)[ , c("Cell.X.Position", "Cell.Y.Position")]
    
    spe_slices_list[[i]] <- spe_slice
  }
  return(spe_slices_list)
}


### Spes 3D setup --------------------------------------------

setwd("~/R/SPIAT-3D_benchmarking/simulations_and_analysis_S1/S1_data")
spes_metadata <- readRDS("spes_metadata.RDS")

# Get number of spe groups
n_spe_groups <- length(spes_metadata)

# Get number of spes in each group (should be the same)
n_spes <- length(spes_metadata[[1]])

# Define AMD data frames as well as constants
cell_types <- c("A", "B")

AMD_pairs <- c("A/A", "A/B", "B/A", "B/B")
AMD_colnames <- c("spe", "reference", "target", "AMD")
AMD_df <- data.frame(matrix(nrow = n_spes * length(AMD_pairs), ncol = 4))
colnames(AMD_df) <- AMD_colnames


# Define MS, NMS, ACIN, ACINP, CKR, AE data frames as well as constants
radii <- seq(20, 100, 10)
radii_colnames <- paste("r", radii, sep = "")

MS_colnames <- c("spe", "reference", "target", radii_colnames)
MS_df <- data.frame(matrix(nrow = n_spes * length(cell_types), ncol = length(MS_colnames)))
colnames(MS_df) <- MS_colnames

# NMS has same data frame output as MS
NMS_df <- MS_df

# Only choose prop(A) as prop(A) = 1 - prop(B) always
ACINP_colnames <- c("spe", "reference", "target", radii_colnames)
ACINP_df <- data.frame(matrix(nrow = n_spes * length(cell_types), ncol = length(ACINP_colnames)))
colnames(ACINP_df) <- ACINP_colnames

# AE has same data frame output as ACINP
AE_df <- ACINP_df

## ACIN and CKR are twice as large
# (ref A and tar A or B) OR (ref B and tar B or A)
ACIN_colnames <- c("spe", "reference", "target", radii_colnames)
ACIN_df <- data.frame(matrix(nrow = n_spes * length(cell_types)^2, ncol = length(ACIN_colnames)))
colnames(ACIN_df) <- ACIN_colnames

# CKR, CLR, COO, CGR have same data frame ouptut as ACIN
CKR_df <- ACIN_df
CLR_df <- ACIN_df
COO_df <- ACIN_df
CGR_df <- ACIN_df

# Define SAC and prevalence data frames as well as constants
n_splits <- 10
thresholds <- seq(0.01, 1, 0.01)
thresholds_colnames <- paste("t", thresholds, sep = "")

prop_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))

PBSAC_df_colnames <- c("spe", "reference", "target", "PBSAC")
PBSAC_df <- data.frame(matrix(nrow = n_spes * nrow(prop_cell_types), ncol = length(PBSAC_df_colnames)))
colnames(PBSAC_df) <- PBSAC_df_colnames

PBP_df_colnames <- c("spe", "reference", "target", thresholds_colnames)
PBP_df <- data.frame(matrix(nrow = n_spes * nrow(prop_cell_types), ncol = length(PBP_df_colnames)))
colnames(PBP_df) <- PBP_df_colnames


entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))

EBSAC_df_colnames <- c("spe", "cell_types", "EBSAC")
EBSAC_df <- data.frame(matrix(nrow = n_spes * nrow(entropy_cell_types), ncol = length(EBSAC_df_colnames)))
colnames(EBSAC_df) <- EBSAC_df_colnames

EBP_df_colnames <- c("spe", "cell_types", thresholds_colnames)
EBP_df <- data.frame(matrix(nrow = n_spes * nrow(entropy_cell_types), ncol = length(EBP_df_colnames)))
colnames(EBP_df) <- EBP_df_colnames


# Add all to list:
metric_df_list3D <- list(AMD = AMD_df,
                         MS = MS_df,
                         NMS = NMS_df,
                         ACINP = ACINP_df,
                         AE = AE_df,
                         ACIN = ACIN_df,
                         CKR = CKR_df,
                         CLR = CLR_df,
                         COO = COO_df,
                         CGR = CGR_df,
                         PBSAC = PBSAC_df,
                         PBP = PBP_df,
                         EBSAC = EBSAC_df,
                         EBP = EBP_df)

metric_df_lists3D <- list(mixed_ellipsoid = metric_df_list3D,
                          mixed_network = metric_df_list3D,
                          ringed_ellipsoid = metric_df_list3D,
                          ringed_network = metric_df_list3D,
                          separated_ellipsoid = metric_df_list3D,
                          separated_network = metric_df_list3D)



### Spes 2D setup --------------------------------------------

slice_bottom_z_coords <- c(145, 175, 205)
slice_top_z_coords <- slice_bottom_z_coords + 10
n_slices <- length(slice_bottom_z_coords)

metric_df_list2D <- metric_df_list3D

for (df in metric_df_list2D) {
  df <- cbind(df[1], slice = NA, df[-1])
  df <- df[rep(1:nrow(df), 3), ]
}

metric_df_lists2D <- list(mixed_ellipsoid = metric_df_list2D,
                          mixed_network = metric_df_list2D,
                          ringed_ellipsoid = metric_df_list2D,
                          ringed_network = metric_df_list2D,
                          separated_ellipsoid = metric_df_list2D,
                          separated_network = metric_df_list2D)






### spe 3D and 2D analysis ------------------------------------------
arrangements <- c("separated")
shapes <- c("network")

for (arrangement in arrangements) {
  
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    for (i in seq_len(n_spes)) {
      print(i)
      ### 3D analysis -----------------------------
      spe <- simulate_spe_metadata3D(spes_metadata[[spes_metadata_index]][[i]], plot_image = F)
      spe_name <- paste("spe_", i, sep = "")  
      
      minimum_distance_data <- calculate_minimum_distances_between_cell_types3D(spe,
                                                                                cell_types,
                                                                                show_summary = F,
                                                                                plot_image = F)
      
      minimum_distance_data_summary <- summarise_distances_between_cell_types3D(minimum_distance_data)
      ## Fill in 4 rows at a time for AMD df (as we have A/A, A/B, B/A, B/B)
      index <- 4 * (i - 1) + 1 # index is 1, 5, 9, 13
      
      metric_df_lists3D[[spes_metadata_index]][["AMD"]][index:(index + 3), "spe"] <- spe_name
      metric_df_lists3D[[spes_metadata_index]][["AMD"]][index:(index + 3), "reference"] <- minimum_distance_data_summary$reference
      metric_df_lists3D[[spes_metadata_index]][["AMD"]][index:(index + 3), "target"] <- minimum_distance_data_summary$target
      metric_df_lists3D[[spes_metadata_index]][["AMD"]][index:(index + 3), "AMD"] <- minimum_distance_data_summary$mean
      
      
      
      
      
      index1 <- 2 * (i - 1) + 1 # index1 is 1, 3, 5, ...
      index2 <- 4 * (i - 1) + 1 # index2 is 1, 5, 9, 13...
      for (reference_cell_type in cell_types) {
        gradient_data <- calculate_all_gradient_cc_metrics3D(spe,
                                                             reference_cell_type,
                                                             cell_types,
                                                             radii,
                                                             plot_image = F)
        
        target_cell_type <- setdiff(cell_types, reference_cell_type)
        
        metric_df_lists3D[[spes_metadata_index]][["MS"]][index1, c("spe", "reference", "target")] <- c(spe_name, reference_cell_type, target_cell_type)
        metric_df_lists3D[[spes_metadata_index]][["MS"]][index1, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$mixing_score
        
        metric_df_lists3D[[spes_metadata_index]][["NMS"]][index1, c("spe", "reference", "target")] <- c(spe_name, reference_cell_type, target_cell_type)
        metric_df_lists3D[[spes_metadata_index]][["NMS"]][index1, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$normalised_mixing_score
        
        metric_df_lists3D[[spes_metadata_index]][["ACINP"]][index1, c("spe", "reference", "target")] <- c(spe_name, reference_cell_type, "B")
        metric_df_lists3D[[spes_metadata_index]][["ACINP"]][index1, radii_colnames] <- gradient_data[["cells_in_neighbourhood_proportion"]][["B"]]
        
        metric_df_lists3D[[spes_metadata_index]][["AE"]][index1, c("spe", "reference", "target")] <- c(spe_name, reference_cell_type, "A,B")
        metric_df_lists3D[[spes_metadata_index]][["AE"]][index1, radii_colnames] <- gradient_data[["entropy"]]$entropy
        
        index1 <- index1 + 1
        
        for (target_cell_type in cell_types) {
          ## Calculate ACIN, CKR, CKL, COO as target cell type can also be the reference cell type
          
          # ACIN
          metric_df_lists3D[[spes_metadata_index]][["ACIN"]][index2, c("spe", "reference", "target")] <- c(spe_name, reference_cell_type, target_cell_type)
          metric_df_lists3D[[spes_metadata_index]][["ACIN"]][index2, radii_colnames] <- gradient_data[["cells_in_neighbourhood"]][[target_cell_type]]
          
          # CKR
          metric_df_lists3D[[spes_metadata_index]][["CKR"]][index2, c("spe", "reference", "target")] <- c(spe_name, reference_cell_type, target_cell_type)
          metric_df_lists3D[[spes_metadata_index]][["CKR"]][index2, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]] / gradient_data[["cross_K"]][["expected"]]
          
          # CLR
          metric_df_lists3D[[spes_metadata_index]][["CLR"]][index2, c("spe", "reference", "target")] <- c(spe_name, reference_cell_type, target_cell_type)
          metric_df_lists3D[[spes_metadata_index]][["CLR"]][index2, radii_colnames] <- gradient_data[["cross_L"]][[target_cell_type]] / gradient_data[["cross_L"]][["expected"]]
          
          # COO
          metric_df_lists3D[[spes_metadata_index]][["COO"]][index2, c("spe", "reference", "target")] <- c(spe_name, reference_cell_type, target_cell_type)
          metric_df_lists3D[[spes_metadata_index]][["COO"]][index2, radii_colnames] <- gradient_data[["co_occurrence"]][[target_cell_type]]
          
          # CGR
          metric_df_lists3D[[spes_metadata_index]][["CGR"]][index2, c("spe", "reference", "target")] <- c(spe_name, reference_cell_type, target_cell_type)
          metric_df_lists3D[[spes_metadata_index]][["CGR"]][index2, radii_colnames] <- gradient_data[["cross_G"]][[target_cell_type]][["observed_cross_G"]] / gradient_data[["cross_G"]][[target_cell_type]][["expected_cross_G"]]
          
          index2 <- index2 + 1
        }
      }
      
      
      # Get proportion grid metrics
      for (j in seq_len(nrow(prop_cell_types))) {
        proportion_grid_metrics <- calculate_cell_proportion_grid_metrics3D(spe, 
                                                                            n_splits,
                                                                            strsplit(prop_cell_types$ref[j], ",")[[1]], 
                                                                            strsplit(prop_cell_types$tar[j], ",")[[1]],
                                                                            plot_image = F)
        
        PBSAC <- calculate_spatial_autocorrelation3D(proportion_grid_metrics, 
                                                     "proportion",
                                                     weight_method = 0.1)
        
        PBP_df <- calculate_prevalence_gradient3D(proportion_grid_metrics,
                                                  "proportion",
                                                  show_AUC = F,
                                                  plot_image = F)
        
        index <- nrow(prop_cell_types) * (i - 1) + j
        metric_df_lists3D[[spes_metadata_index]][["PBSAC"]][index, c("spe", "reference", "target")] <- c(spe_name, prop_cell_types$ref[j], prop_cell_types$tar[j])
        metric_df_lists3D[[spes_metadata_index]][["PBSAC"]][index, "PBSAC"] <- PBSAC
        
        metric_df_lists3D[[spes_metadata_index]][["PBP"]][index, c("spe", "reference", "target")] <- c(spe_name, prop_cell_types$ref[j], prop_cell_types$tar[j])
        metric_df_lists3D[[spes_metadata_index]][["PBP"]][index, thresholds_colnames] <- PBP_df$prevalence
      }
      
      # Get entropy grid metrics
      for (j in seq_len(nrow(entropy_cell_types))) {
        entropy_grid_metrics <- calculate_entropy_grid_metrics3D(spe, 
                                                                 n_splits,
                                                                 strsplit(entropy_cell_types$cell_types[j], ",")[[1]], 
                                                                 plot_image = F)
        
        EBSAC <- calculate_spatial_autocorrelation3D(entropy_grid_metrics, 
                                                     "entropy",
                                                     weight_method = 0.1)
        
        EBP_df <- calculate_prevalence_gradient3D(entropy_grid_metrics,
                                                  "entropy",
                                                  show_AUC = F,
                                                  plot_image = F)
        
        index <- nrow(entropy_cell_types) * (i - 1) + j
        metric_df_lists3D[[spes_metadata_index]][["EBSAC"]][index, c("spe", "cell_types")] <- c(spe_name, entropy_cell_types$cell_types[j])
        metric_df_lists3D[[spes_metadata_index]][["EBSAC"]][index, "EBSAC"] <- EBSAC
        
        metric_df_lists3D[[spes_metadata_index]][["EBP"]][index, c("spe", "cell_types")] <- c(spe_name, entropy_cell_types$cell_types[j])
        metric_df_lists3D[[spes_metadata_index]][["EBP"]][index, thresholds_colnames] <- EBP_df$prevalence
      }  
      
      
      ### 2D slicing analysis -------------------------
      spe_slices <- get_spe_slices_list(spe, slice_bottom_z_coords, slice_top_z_coords)
      for (slice_index in seq(n_slices)) {
        spe_slice <- spe_slices[[slice_index]]
        
        minimum_distance_data <- calculate_minimum_distances_between_cell_types2D(spe_slice,
                                                                                  cell_types,
                                                                                  show_summary = F,
                                                                                  plot_image = F)
        
        minimum_distance_data_summary <- summarise_distances_between_cell_types2D(minimum_distance_data)
        ## Fill in 4 rows at a time for AMD df (as we have A/A, A/B, B/A, B/B)
        index <- n_slices * 4 * (i - 1) + 4 * (slice_index - 1) + 1
        metric_df_lists2D[[spes_metadata_index]][["AMD"]][index:(index + 3), "spe"] <- spe_name
        metric_df_lists2D[[spes_metadata_index]][["AMD"]][index:(index + 3), "slice"] <- slice_index
        metric_df_lists2D[[spes_metadata_index]][["AMD"]][index:(index + 3), "reference"] <- minimum_distance_data_summary$reference
        metric_df_lists2D[[spes_metadata_index]][["AMD"]][index:(index + 3), "target"] <- minimum_distance_data_summary$target
        metric_df_lists2D[[spes_metadata_index]][["AMD"]][index:(index + 3), "AMD"] <- minimum_distance_data_summary$mean
        
        
        index1 <- n_slices * 2 * (i - 1) + 2 * (slice_index - 1) + 1
        index2 <- n_slices * 4 * (i - 1) + 4 * (slice_index - 1) + 1
        for (reference_cell_type in cell_types) {
          gradient_data <- calculate_all_gradient_cc_metrics2D(spe_slice,
                                                               reference_cell_type,
                                                               cell_types,
                                                               radii,
                                                               plot_image = F)
          
          target_cell_type <- setdiff(cell_types, reference_cell_type)
          
          metric_df_lists2D[[spes_metadata_index]][["MS"]][index1, c("spe", "slice", "reference", "target")] <- c(spe_name, slice_index, reference_cell_type, target_cell_type)
          metric_df_lists2D[[spes_metadata_index]][["NMS"]][index1, c("spe", "slice", "reference", "target")] <- c(spe_name, slice_index, reference_cell_type, target_cell_type)
          metric_df_lists2D[[spes_metadata_index]][["ACINP"]][index1, c("spe", "slice","reference", "target")] <- c(spe_name, slice_index, reference_cell_type, "B")
          metric_df_lists2D[[spes_metadata_index]][["AE"]][index1, c("spe", "slice","reference", "target")] <- c(spe_name, slice_index, reference_cell_type, "A,B")
          
          if (!is.null(gradient_data)) {
            metric_df_lists2D[[spes_metadata_index]][["MS"]][index1, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$mixing_score
            metric_df_lists2D[[spes_metadata_index]][["NMS"]][index1, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$normalised_mixing_score
            metric_df_lists2D[[spes_metadata_index]][["ACINP"]][index1, radii_colnames] <- gradient_data[["cells_in_neighbourhood_proportion"]][["B"]]
            metric_df_lists2D[[spes_metadata_index]][["AE"]][index1, radii_colnames] <- gradient_data[["entropy"]]$entropy        
          }
          else {
            metric_df_lists2D[[spes_metadata_index]][["MS"]][index1, radii_colnames] <- NA
            metric_df_lists2D[[spes_metadata_index]][["NMS"]][index1, radii_colnames] <- NA
            metric_df_lists2D[[spes_metadata_index]][["ACINP"]][index1, radii_colnames] <- NA
            metric_df_lists2D[[spes_metadata_index]][["AE"]][index1, radii_colnames] <- NA
          }
          
          
          index1 <- index1 + 1
          
          for (target_cell_type in cell_types) {
            ## Calculate ACIN, CKR, CLR, COO as target cell type can also be the reference cell type
            
            metric_df_lists2D[[spes_metadata_index]][["ACIN"]][index2, c("spe", "slice", "reference", "target")] <- c(spe_name, slice_index, reference_cell_type, target_cell_type)
            metric_df_lists2D[[spes_metadata_index]][["CKR"]][index2, c("spe", "slice", "reference", "target")] <- c(spe_name, slice_index, reference_cell_type, target_cell_type)
            metric_df_lists2D[[spes_metadata_index]][["CLR"]][index2, c("spe", "slice", "reference", "target")] <- c(spe_name, slice_index, reference_cell_type, target_cell_type)
            metric_df_lists2D[[spes_metadata_index]][["COO"]][index2, c("spe", "slice", "reference", "target")] <- c(spe_name, slice_index, reference_cell_type, target_cell_type)
            metric_df_lists2D[[spes_metadata_index]][["CGR"]][index2, c("spe", "slice", "reference", "target")] <- c(spe_name, slice_index, reference_cell_type, target_cell_type)
            
            if (!is.null(gradient_data)) {
              metric_df_lists2D[[spes_metadata_index]][["ACIN"]][index2, radii_colnames] <- gradient_data[["cells_in_neighbourhood"]][[target_cell_type]]
              metric_df_lists2D[[spes_metadata_index]][["CKR"]][index2, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]] / gradient_data[["cross_K"]][["expected"]]
              metric_df_lists2D[[spes_metadata_index]][["CLR"]][index2, radii_colnames] <- gradient_data[["cross_L"]][[target_cell_type]] / gradient_data[["cross_L"]][["expected"]]
              metric_df_lists2D[[spes_metadata_index]][["COO"]][index2, radii_colnames] <- gradient_data[["co_occurrence"]][[target_cell_type]]
              metric_df_lists2D[[spes_metadata_index]][["CGR"]][index2, radii_colnames] <- gradient_data[["cross_G"]][[target_cell_type]][["observed_cross_G"]] / gradient_data[["cross_G"]][[target_cell_type]][["expected_cross_G"]]
            }
            else {
              metric_df_lists2D[[spes_metadata_index]][["ACIN"]][index2, radii_colnames] <- NA
              metric_df_lists2D[[spes_metadata_index]][["CKR"]][index2, radii_colnames] <- NA
              metric_df_lists2D[[spes_metadata_index]][["CLR"]][index2, radii_colnames] <- NA
              metric_df_lists2D[[spes_metadata_index]][["COO"]][index2, radii_colnames] <- NA
              metric_df_lists2D[[spes_metadata_index]][["CGR"]][index2, radii_colnames] <- NA
            }
            
            index2 <- index2 + 1
          }
        }
        
        
        # Get proportion grid metrics
        for (j in seq_len(nrow(prop_cell_types))) {
          proportion_grid_metrics <- calculate_cell_proportion_grid_metrics2D(spe_slice, 
                                                                              n_splits,
                                                                              strsplit(prop_cell_types$ref[j], ",")[[1]], 
                                                                              strsplit(prop_cell_types$tar[j], ",")[[1]],
                                                                              plot_image = F)
          
          if (!is.null(proportion_grid_metrics)) {
            PBSAC <- calculate_spatial_autocorrelation2D(proportion_grid_metrics, 
                                                         "proportion",
                                                         weight_method = 0.1)
            
            PBP_df <- calculate_prevalence_gradient2D(proportion_grid_metrics,
                                                      "proportion",
                                                      show_AUC = F,
                                                      plot_image = F)
          }
          else {
            PBSAC <- NA
            PBP_df <- data.frame(threshold = seq(0.01, 1, 0.01), prevalence = NA)
          }
          
          
          
          index <- n_slices * nrow(prop_cell_types) * (i - 1) + nrow(prop_cell_types) * (slice_index - 1) + j
          
          metric_df_lists2D[[spes_metadata_index]][["PBSAC"]][index, c("spe", "slice", "reference", "target")] <- c(spe_name, slice_index, prop_cell_types$ref[j], prop_cell_types$tar[j])
          metric_df_lists2D[[spes_metadata_index]][["PBSAC"]][index, "PBSAC"] <- PBSAC
          
          metric_df_lists2D[[spes_metadata_index]][["PBP"]][index, c("spe", "slice", "reference", "target")] <- c(spe_name, slice_index, prop_cell_types$ref[j], prop_cell_types$tar[j])
          metric_df_lists2D[[spes_metadata_index]][["PBP"]][index, thresholds_colnames] <- PBP_df$prevalence
        }
        
        # Get entropy grid metrics
        for (j in seq_len(nrow(entropy_cell_types))) {
          entropy_grid_metrics <- calculate_entropy_grid_metrics2D(spe_slice, 
                                                                   n_splits,
                                                                   strsplit(entropy_cell_types$cell_types[j], ",")[[1]], 
                                                                   plot_image = F)
          
          if (!is.null(entropy_grid_metrics)) {
            EBSAC <- calculate_spatial_autocorrelation2D(entropy_grid_metrics, 
                                                         "entropy",
                                                         weight_method = 0.1)
            
            EBP_df <- calculate_prevalence_gradient2D(entropy_grid_metrics,
                                                      "entropy",
                                                      show_AUC = F,
                                                      plot_image = F)
          }
          else {
            EBSAC <- NA
            EBP_df <- data.frame(threshold = seq(0.01, 1, 0.01), prevalence = NA)
          }
          
          index <- n_slices * nrow(entropy_cell_types) * (i - 1) + nrow(entropy_cell_types) * (slice_index - 1) + j
          
          metric_df_lists2D[[spes_metadata_index]][["EBSAC"]][index, c("spe", "slice", "cell_types")] <- c(spe_name, slice_index, entropy_cell_types$cell_types[j])
          metric_df_lists2D[[spes_metadata_index]][["EBSAC"]][index, "EBSAC"] <- EBSAC
          
          metric_df_lists2D[[spes_metadata_index]][["EBP"]][index, c("spe", "slice", "cell_types")] <- c(spe_name, slice_index, entropy_cell_types$cell_types[j])
          metric_df_lists2D[[spes_metadata_index]][["EBP"]][index, thresholds_colnames] <- EBP_df$prevalence
        }  
      }
    }
  }
}

setwd("~/R/SPIAT-3D_benchmarking/simulations_and_analysis_S1/S1_data")
saveRDS(metric_df_lists3D, "SN_metric_df_lists3D.RDS")
saveRDS(metric_df_lists2D, "SN_metric_df_lists2D.RDS")


