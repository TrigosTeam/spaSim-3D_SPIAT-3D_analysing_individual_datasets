# Data entry (THE ONLY PART YOU NEED TO CHANGE) -----

setwd("***directory to your data***")
data3D <- read.csv("***your data.csv***")

cell_types <- "***whatever you want***"

save_directory <- "***Where ever you want to save***" # e.g. "~/R/SPIAT-3D_benchmarking/public_3D_data_analysis/metric_df_lists"
file_name <- "*** file save name.RDS ***" # e.g. openST_human_metastatic_lymph_node_metric_df_list.RDS

# SPIAT-3D functions -----
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
library(alphashape3d)

alpha_hull_clustering3D <- function(spatial_df, 
                                    cell_types_of_interest, 
                                    alpha, 
                                    minimum_cells_in_cluster,
                                    feature_colname = "Cell.Type", 
                                    plot_image = T) {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
  }
  ## Check cell types of interest are found in the spatial_df object
  unknown_cell_types <- setdiff(cell_types_of_interest, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
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
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  ## Subset for the chosen cell_types_of_interest
  spatial_df_subset <- spatial_df[ , spatial_df[[feature_colname]] %in% cell_types_of_interest]
  spatial_df_subset_coords <- spatial_df_subset[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
  
  ## Get the alpha hull
  alpha_hull <- ashape3d(as.matrix(spatial_df_subset_coords), alpha = alpha)
  
  if (sum(alpha_hull$triang[, 9]) == 0) stop("alpha value is too small? No alpha hulls identified")
  
  ## Determine which alpha hull cluster each cell_type_of_interest belongs to
  alpha_hull_clusters <- components_ashape3d(alpha_hull)
  
  df_cell_types_of_interest <- spatial_df[spatial_df[[feature_colname]] %in% cell_types_of_interest, ]
  df_other_cell_types <- spatial_df[!(spatial_df[[feature_colname]] %in% cell_types_of_interest), ]
  
  df_cell_types_of_interest$alpha_hull_cluster <- alpha_hull_clusters
  df_other_cell_types$alpha_hull_cluster <- 0
  
  ## Ignore cell_types_of_interest which belong to an alpha hull cluster with less than minimum_cells_in_cluster
  alpha_hull_clusters_table <- table(alpha_hull_clusters)
  maximium_alpha_hull_cluster <- Position(function(x) x < minimum_cells_in_cluster, alpha_hull_clusters_table)
  maximium_alpha_hull_cluster <- as.numeric(names(alpha_hull_clusters_table[maximium_alpha_hull_cluster]))
  
  if (!is.na(maximium_alpha_hull_cluster) && maximium_alpha_hull_cluster != -1) {
    spatial_df_subset_coords <- spatial_df_subset_coords[alpha_hull_clusters >= 1 & alpha_hull_clusters < maximium_alpha_hull_cluster, ]
    
    df_cell_types_of_interest$alpha_hull_cluster <- ifelse(alpha_hull_clusters >= 1 & alpha_hull_clusters < maximium_alpha_hull_cluster, 
                                                           alpha_hull_clusters, 0)
    
    ## Get the alpha hull again...
    alpha_hull <- ashape3d(as.matrix(spatial_df_subset_coords), alpha = alpha)
  }
  
  ## Convert data frame to spatial_df object
  df <- rbind(df_cell_types_of_interest, df_other_cell_types)
  
  spatial_df <- SpatialExperiment(
    assay = matrix(data = NA, nrow = nrow(df), ncol = nrow(df)),
    colData = df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = spatial_df@metadata)
  
  ## Get the information of the vertices and faces of the alpha hull (what 3 vertices make up each face triangle?)
  vertices <- alpha_hull$x
  faces <- alpha_hull$triang[alpha_hull$triang[, 9] == 2, c("tr1", "tr2", "tr3")]
  spatial_df@metadata$alpha_hull <- list(vertices = vertices, faces = faces, ashape3d_object = alpha_hull)
  
  ## Plot
  if (plot_image) {
    fig <- plot_alpha_hull_clusters3D(spatial_df, feature_colname = feature_colname)
    methods::show(fig)
  }
  
  return(spatial_df)
}
calculate_all_gradient_cc_metrics3D <- function(spatial_df, 
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
                 "entropy" = data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types))),
                 "cross_K" = data.frame(matrix(nrow = length(radii), ncol = length(cross_K_df_colnames))),
                 "cross_L" = data.frame(matrix(nrow = length(radii), ncol = length(cross_K_df_colnames))),
                 "cross_G" = list(),
                 "co_occurrence" = data.frame(matrix(nrow = length(radii), ncol = length(co_occurrence_df_colnames))))
  colnames(result[["cells_in_neighbourhood"]]) <- target_cell_types
  colnames(result[["cells_in_neighbourhood_proportion"]]) <- target_cell_types
  colnames(result[["entropy"]]) <- target_cell_types
  colnames(result[["cross_K"]]) <- cross_K_df_colnames
  colnames(result[["cross_L"]]) <- cross_K_df_colnames
  colnames(result[["co_occurrence"]]) <- co_occurrence_df_colnames
  
  # Define indiviudal data frames for mixing_score and cross_G
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
    df <- calculate_all_single_radius_cc_metrics3D(spatial_df,
                                                   reference_cell_type,
                                                   target_cell_types,
                                                   radii[i],
                                                   feature_colname)
    
    if (is.null(df)) return(NULL)
    
    df[["cells_in_neighbourhood"]]$ref_cell_id <- NULL
    
    result[["cells_in_neighbourhood"]][i, ] <- apply(df[["cells_in_neighbourhood"]], 2, mean)
    result[["cells_in_neighbourhood_proportion"]][i, ] <- apply(df[["cells_in_neighbourhood_proportion"]][ , paste(target_cell_types, "_prop", sep = "")], 2, mean, na.rm = T)
    result[["entropy"]][i, ] <- apply(df[["entropy"]][ , paste(target_cell_types, "_entropy", sep = "")], 2, mean, na.rm = T)
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
    
    expected_entropy <- calculate_entropy_background3D(spatial_df, target_cell_types, feature_colname)
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

calculate_all_single_radius_cc_metrics3D <- function(spatial_df, 
                                                     reference_cell_type, 
                                                     target_cell_types, 
                                                     radius, 
                                                     feature_colname = "Cell.Type") {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
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
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  
  ## For reference_cell_type, check it is found in the spatial_df object
  if (!(reference_cell_type %in% spatial_df[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spatial_df object"))
    return(NULL)
  }
  ## For target_cell_types, check they are found in the spatial_df object
  unknown_cell_types <- setdiff(target_cell_types, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spatial_df object:\n   ",
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
  spatial_df_coords <- spatial_df[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
  length <- round(max(spatial_df_coords$Cell.X.Position) - min(spatial_df_coords$Cell.X.Position))
  width  <- round(max(spatial_df_coords$Cell.Y.Position) - min(spatial_df_coords$Cell.Y.Position))
  height <- round(max(spatial_df_coords$Cell.Z.Position) - min(spatial_df_coords$Cell.Z.Position))
  ## Get volume of the window the cells are in
  volume <- length * width * height
  
  # All single radius cc metrics stem from calculate_entropy3D function
  entropy_df <- calculate_entropy3D(spatial_df, 
                                    reference_cell_type, 
                                    target_cell_types, 
                                    radius, 
                                    feature_colname)  
  
  ## Cells in neighbourhood ----------
  result[["cells_in_neighbourhood"]] <- entropy_df[ , c("ref_cell_id", target_cell_types)]
  
  ## Cells in neighbourhood proportion ----------
  result[["cells_in_neighbourhood_proportion"]] <- entropy_df[ , c("ref_cell_id", paste(target_cell_types, "_prop", sep = ""))]
  
  ## Entropy --------------
  result[["entropy"]] <- entropy_df[ , c("ref_cell_id", paste(target_cell_types, "_entropy", sep = ""))]
  
  ## Mixing score -----------------
  for (target_cell_type in target_cell_types) {
    mixing_score_df <- data.frame(matrix(nrow = 1, ncol = length(mixing_score_df_colnames)))
    colnames(mixing_score_df) <- mixing_score_df_colnames
    mixing_score_df$ref_cell_type <- reference_cell_type
    
    # No need to fill in mixing_score_df if the reference and target cell is the same
    if (reference_cell_type != target_cell_type) {
      mixing_score_df$tar_cell_type <- target_cell_type
      mixing_score_df$n_ref_cells <- sum(spatial_df[[feature_colname]] == reference_cell_type)
      mixing_score_df$n_tar_cells <- sum(spatial_df[[feature_colname]] == target_cell_type)
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
    cross_K_df[[target_cell_type]] <- (((volume * sum(entropy_df[[target_cell_type]])) / sum(spatial_df[[feature_colname]] == reference_cell_type)) / sum(spatial_df[[feature_colname]] == target_cell_type)) 
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
    n_target_cells <- sum(spatial_df[[feature_colname]] == target_cell_type)
    target_cell_type_intensity <- n_target_cells / volume
    observed_cross_G <- sum(reference_target_interactions != 0) / length(reference_target_interactions)
    expected_cross_G <- 1 - exp(-1 * target_cell_type_intensity * (4 / 3) * pi * radius^3)
    
    cross_G_df$observed_cross_G <- observed_cross_G
    cross_G_df$expected_cross_G <- expected_cross_G
    result[["cross_G"]][[target_cell_type]] <- cross_G_df
  } 
  
  ## Co_occurrence ---------------
  all_cell_types <- unique(spatial_df[[feature_colname]])
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood3D(spatial_df,
                                                                  reference_cell_type,
                                                                  all_cell_types,
                                                                  radius,
                                                                  feature_colname,
                                                                  F,
                                                                  F)
  
  cells_in_neighbourhood_df$total <- rowSums(cells_in_neighbourhood_df[, -1], na.rm = TRUE)
  
  co_occurrence_df <- data.frame(matrix(nrow = 1, ncol = length(co_occurrence_df_colnames)))
  colnames(co_occurrence_df) <- co_occurrence_df_colnames
  co_occurrence_df$reference <- reference_cell_type
  
  n_cells_in_spe <- length(spatial_df[[feature_colname]])
  n_cells_in_reference_cell_type_radius <- sum(cells_in_neighbourhood_df$total)
  
  for (target_cell_type in target_cell_types) {
    n_target_cells_in_reference_cell_type_radius <- sum(cells_in_neighbourhood_df[[target_cell_type]])
    target_cell_type_proportion_in_reference_cell_type_radius <- n_target_cells_in_reference_cell_type_radius / n_cells_in_reference_cell_type_radius
    n_target_cells_in_spe <- sum(spatial_df[[feature_colname]] == target_cell_type)
    target_cell_type_proportion_in_spe <- n_target_cells_in_spe / n_cells_in_spe
    target_cell_type_co_occurrence <- target_cell_type_proportion_in_reference_cell_type_radius / target_cell_type_proportion_in_spe
    
    co_occurrence_df[[target_cell_type]] <- target_cell_type_co_occurrence
  }
  result[["co_occurrence"]] <- co_occurrence_df
  
  return(result)
}


calculate_border_of_clusters3D <- function(spatial_df, 
                                           radius,
                                           cluster_colname, 
                                           feature_colname = "Cell.Type", 
                                           plot_image = T) {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
  }
  if (!(is.numeric(radius) && length(radius) == 1 && radius > 0)) {
    stop("`radius` is not a positive numeric.")
  }
  if (!is.character(cluster_colname)) {
    stop("`cluster_colname` is not a character. This should be 'alpha_hull_cluster', 'dbscan_cluster', or 'grid_based_cluster', depending on the chosen method.")
  }
  if (is.null(spatial_df[[cluster_colname]])) {
    stop(paste("No column called", cluster_colname, "found in spatial_df object."))
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  ## Get spatial coords of spatial_df
  spatial_df_coords <- spatial_df[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
  
  ## Get coords of non-cluster cells
  non_cluster_coords <- spatial_df_coords[spatial_df[[cluster_colname]] == 0, ]
  
  # New column for spatial_df object: 'cluster_border'. Default is 'outside'
  spatial_df$cluster_border <- "outside"
  
  # Label cells part of a cluster (e.g. 'cluster1')
  spatial_df$cluster_border[spatial_df[[cluster_colname]] != 0] <- paste("inside_C", spatial_df[[cluster_colname]][spatial_df[[cluster_colname]] != 0], sep = "")
  
  ## Iterate for each cluster
  n_clusters <- max(spatial_df[[cluster_colname]])
  
  for (i in seq_len(n_clusters)) {
    
    ## Subset for cells in the current cluster of interest
    cluster_coords <- spatial_df_coords[spatial_df[[cluster_colname]] == i, ]
    
    # For each cell in the current cluster, check how many other cells in the cluster are in its radius
    cluster_to_cluster_interactions <- dbscan::frNN(cluster_coords, radius)
    
    # Determine the median minimum number of cluster cells found in the radius of cluster cell. Use this as the threshold for non-cluster cells.
    non_cluster_threshold <- quantile(unlist(lapply(cluster_to_cluster_interactions$dist, length)), 0.5)
    
    # For each non-cluster cell, check how many cluster cells are in its radius.
    non_cluster_to_cluster_interactions <- dbscan::frNN(cluster_coords, radius, non_cluster_coords)
    
    # If number of cluster cells found in the radius of non-cluster cells is greater than threshold, non-cluster cell has probably infiltrated cluster too
    n_cluster_cells_in_non_cluster_cell_radius <- unlist(lapply(non_cluster_to_cluster_interactions$id, length))
    
    spatial_df$cluster_border[as.numeric(names(non_cluster_to_cluster_interactions$id)[n_cluster_cells_in_non_cluster_cell_radius > non_cluster_threshold])] <- paste("infiltrated_C", i, sep = "")
    
    # If number of cluster cells found in the radius of non-cluster cells is less than threshold, but greater than 0, non-cluster cell is probably on the border
    spatial_df$cluster_border[as.numeric(names(non_cluster_to_cluster_interactions$id)[n_cluster_cells_in_non_cluster_cell_radius > 0 & n_cluster_cells_in_non_cluster_cell_radius < non_cluster_threshold])] <- paste("border_C", i, sep = "")
  }
  
  ## Plot
  if (plot_image) {
    fig <- plot_cells3D(spatial_df, feature_colname = "cluster_border")
    methods::show(fig)
  }
  
  return(spatial_df)
}
calculate_cell_proportion_grid_metrics3D <- function(spatial_df, 
                                                     n_splits,
                                                     reference_cell_types,
                                                     target_cell_types,
                                                     feature_colname = "Cell.Type",
                                                     plot_image = TRUE) {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  ## Check reference_cell_types are found in the spatial_df object
  unknown_cell_types <- setdiff(reference_cell_types, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in reference_cell_types are not found in the spatial_df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  ## Check target_cell_types are found in the spatial_df object
  unknown_cell_types <- setdiff(target_cell_types, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spatial_df object:\n   ",
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
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Add grid metrics to spatial_df
  grid_metrics <- get_spatial_df_grid_metrics3D(spatial_df, n_splits, feature_colname)
  
  # Get grid_prism_cell_matrix from spatial_df
  grid_prism_cell_matrix <- grid_metrics$grid_prism_cell_matrix
  
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
  result <- cbind(result, grid_metrics$grid_prism_coordinates)
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_metrics_continuous3D(result, "proportion")
    methods::show(fig)
  }
  
  return(result)
}
calculate_cell_proportions_of_clusters3D <- function(spatial_df, cluster_colname, feature_colname = "Cell.Type", plot_image = T) {
  
  # Get number of clusters
  n_clusters <- max(spatial_df[[cluster_colname]])
  
  ## Get different cell types found in the clusters (alphabetical for consistency)
  cell_types <- unique(spatial_df[[feature_colname]][spatial_df[[cluster_colname]] != 0])
  cell_types <- cell_types[order(cell_types)]
  
  ## For each cluster, determine the size and cell proportion of each cluster
  result <- data.frame(matrix(nrow = n_clusters, ncol = 2 + length(cell_types)))
  colnames(result) <- c("cluster_number", "n_cells", cell_types)
  result$cluster_number <- as.character(seq(n_clusters))
  
  for (i in seq(n_clusters)) {
    cells_in_cluster <- spatial_df[[feature_colname]][spatial_df[[cluster_colname]] == i]
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
#' @param spatial_df A SpatialExperiment object containing 3D spatial information for the cells.
#' @param cell_types_of_interest A character vector spatial_dfcifying the cell types of interest.
#'    If NULL, all cell types in the `feature_colname` column will be considered.
#' @param feature_colname A string spatial_dfcifying the name of the column in the `colData` slot of the SpatialExperiment
#'    object that contains the cell type information.
#' @param plot_image A logical indicating whether to plot violin plots of the minimum distances 
#'    between cell type pairs. Defaults to TRUE.
#'
#' @return A data frame containing the cell types, their frequencies, proportions, and percentages.
#'
#' @examples
#' cell_proportions <- calculate_cell_proportions3D(
#'     spatial_df = SPIAT3D::simulated_spatial_df, 
#'     cell_types_of_interest = NULL, 
#'     feature_colname = "Cell.Type", 
#'     plot_image = TRUE
#' )
#' 
#' @export


calculate_cell_proportions3D <- function(spatial_df,
                                         cell_types_of_interest = NULL, 
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
  }
  if (ncol(spatial_df) == 0) {
    stop("No cells found for calculating cell proportions.")
  }
  if (!(is.null(cell_types_of_interest) || is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Creates frequency/bar plot of all cell types in the entire image
  cell_proportions <- data.frame(table(spatial_df[[feature_colname]]))
  names(cell_proportions) <- c("cell_type", 'frequency')
  
  # Only include cell types the user has chosen
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spatial_df object
    unknown_cell_types <- setdiff(cell_types_of_interest, cell_proportions$cell_type)
    if (length(unknown_cell_types) != 0) {
      stop(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
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
calculate_cells_in_neighbourhood_gradient3D <- function(spatial_df, 
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
    cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood3D(spatial_df,
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
calculate_cells_in_neighbourhood_proportions_gradient3D <- function(spatial_df, 
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
    cell_proportions_neighbourhood_proportions_df <- calculate_cells_in_neighbourhood_proportions3D(spatial_df,
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
calculate_cells_in_neighbourhood_proportions3D <- function(spatial_df, 
                                                           reference_cell_type, 
                                                           target_cell_types, 
                                                           radius, 
                                                           feature_colname = "Cell.Type") {
  
  ## Get cells in neighbourhood df
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood3D(spatial_df,
                                                                  reference_cell_type,
                                                                  c(reference_cell_type, target_cell_types),
                                                                  radius,
                                                                  feature_colname,
                                                                  FALSE,
                                                                  FALSE)
  
  if (is.null(cells_in_neighbourhood_df)) return(NULL)
  
  cells_in_neighbourhood_df[ , paste(target_cell_types, "_prop", sep = "")] <- 
    cells_in_neighbourhood_df[ , target_cell_types] / (cells_in_neighbourhood_df[ , target_cell_types] + cells_in_neighbourhood_df[ , reference_cell_type])
  
  # If reference cell type is in target cell types, proportion should be 1
  if (reference_cell_type %in% target_cell_types) {
    cells_in_neighbourhood_df[cells_in_neighbourhood_df[[reference_cell_type]] != 0, paste(reference_cell_type, "_prop", sep = "")] <- 1
  }
  
  return(cells_in_neighbourhood_df)
}
calculate_cells_in_neighbourhood3D <- function(spatial_df, 
                                               reference_cell_type, 
                                               target_cell_types, 
                                               radius, 
                                               feature_colname = "Cell.Type",
                                               show_summary = TRUE,
                                               plot_image = TRUE) {
  
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
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
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(show_summary)) {
    stop("`show_summary` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  ## For reference_cell_type, check it is found in the spatial_df object
  if (!(reference_cell_type %in% spatial_df[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spatial_df object"))
    return(NULL)
  }
  ## For target_cell_types, check they are found in the spatial_df object
  unknown_cell_types <- setdiff(target_cell_types, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spatial_df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
  }
  
  if (is.null(spatial_df[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spatial_df")
    spatial_df$Cell.ID <- paste("Cell", seq(nrow(spatial_df)), sep = "_")
  }  
  
  # Get spatial_df coords
  spatial_df_coords <- spatial_df[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
  
  # Get reference_cell_type coords
  reference_cell_type_coords <- spatial_df_coords[spatial_df[[feature_colname]] == reference_cell_type, ]
  
  result <- data.frame(ref_cell_id = spatial_df$Cell.ID[spatial_df[[feature_colname]] == reference_cell_type])
  
  for (target_cell_type in target_cell_types) {
    
    if (sum(spatial_df[[feature_colname]] == target_cell_type) == 0) {
      result[[target_cell_type]] <- NA
      next
    }
    
    ## Get target_cell_type coords
    target_cell_type_coords <- spatial_df_coords[spatial_df[[feature_colname]] == target_cell_type, ]
    
    ## Determine number of target cells spatial_dfcified distance for each reference cell
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

calculate_center_of_clusters3D <- function(spatial_df, cluster_colname) {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
  }
  if (!is.character(cluster_colname)) {
    stop("`cluster_colname` is not a character. This should be 'alpha_hull_cluster', 'dbscan_cluster', or 'grid_based_cluster', depending on the chosen method.")
  }
  if (is.null(spatial_df[[cluster_colname]])) {
    stop(paste("No column called", cluster_colname, "found in spatial_df object."))
  }
  
  # Get number of clusters
  n_clusters <- max(spatial_df[[cluster_colname]])
  
  # Get spatial_df coords
  spatial_df_coords <- spatial_df[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
  
  ## For each cluster, determine the number of cells in each cluster of each cluster
  result <- data.frame(matrix(nrow = n_clusters, ncol = 4))
  colnames(result) <- c("cluster_number", "Centre.X.Position", "Centre.Y.Position", "Centre.Z.Position")
  
  result$cluster_number <- as.character(seq(n_clusters))
  for (i in seq(n_clusters)) {
    spatial_df_cluster_coords <- spatial_df_coords[spatial_df[[cluster_colname]] == i, ]
    result[i, c("Centre.X.Position", "Centre.Y.Position", "Centre.Z.Position")] <- 
      apply(spatial_df_cluster_coords, 2, mean)
  }
  
  return(result)
}
calculate_co_occurrence_gradient3D <- function(spatial_df, 
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
    co_occurrence_df <- calculate_co_occurrence3D(spatial_df,
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
calculate_co_occurrence3D <- function(spatial_df, 
                                      reference_cell_type, 
                                      target_cell_types, 
                                      radius, 
                                      feature_colname = "Cell.Type") {
  
  # Get all cell types in spatial_df
  all_cell_types <- unique(spatial_df[[feature_colname]])
  
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood3D(spatial_df,
                                                                  reference_cell_type,
                                                                  all_cell_types,
                                                                  radius,
                                                                  feature_colname,
                                                                  F,
                                                                  F)
  
  cells_in_neighbourhood_df$total <- rowSums(cells_in_neighbourhood_df[, -1], na.rm = TRUE)
  
  result <- data.frame(reference = reference_cell_type)
  
  # Get total number of cells in spatial_df
  n_cells_in_spatial_df <- length(spatial_df[[feature_colname]])
  
  # Get total number of cells in radius around reference cell type
  n_cells_in_reference_cell_type_radius <- sum(cells_in_neighbourhood_df$total)
  
  for (target_cell_type in target_cell_types) {
    
    # Get total number of target cells in radius around reference cell type
    n_target_cells_in_reference_cell_type_radius <- sum(cells_in_neighbourhood_df[[target_cell_type]])
    
    # Get proportion of target cells in radius around reference cell type
    target_cell_type_proportion_in_reference_cell_type_radius <- n_target_cells_in_reference_cell_type_radius / n_cells_in_reference_cell_type_radius
    
    # Get proportion of target cell type in spatial_df
    n_target_cells_in_spatial_df <- sum(spatial_df[[feature_colname]] == target_cell_type)
    target_cell_type_proportion_in_spatial_df <- n_target_cells_in_spatial_df / n_cells_in_spatial_df
    
    # Get co-occurence value for taget cell type
    target_cell_type_co_occurrence <- target_cell_type_proportion_in_reference_cell_type_radius / target_cell_type_proportion_in_spatial_df
    
    # Add to result data frame
    result[[target_cell_type]] <- target_cell_type_co_occurrence
  }
  
  return(result)
}
calculate_cross_G_gradient3D <- function(spatial_df, 
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
    cross_G_df <- calculate_cross_G3D(spatial_df,
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
calculate_cross_G3D <- function(spatial_df,
                                reference_cell_type,
                                target_cell_type,
                                radius,
                                feature_colname = "Cell.Type") {
  
  ### Calculate the observed cross_G
  # Get the number of target cells in the radius around each reference cell
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood3D(spatial_df,
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
  spatial_df_coords <- spatial_df[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
  
  length <- round(max(spatial_df_coords$Cell.X.Position) - min(spatial_df_coords$Cell.X.Position))
  width  <- round(max(spatial_df_coords$Cell.Y.Position) - min(spatial_df_coords$Cell.Y.Position))
  height <- round(max(spatial_df_coords$Cell.Z.Position) - min(spatial_df_coords$Cell.Z.Position))
  
  # Get volume of the window the cells are in
  volume <- length * width * height
  
  # Get the number of target cells
  n_target_cells <- sum(spatial_df[[feature_colname]] == target_cell_type)
  
  # Get target_cell_type intensity (density)
  target_cell_type_intensity <- n_target_cells / volume
  
  # Apply formula
  expected_cross_G <- 1 - exp(-1 * target_cell_type_intensity * (4 / 3) * pi * radius^3)
  
  result <- data.frame(observed_cross_G = observed_cross_G,
                       expected_cross_G = expected_cross_G)
  
  return(result)
}
calculate_cross_K_gradient3D <- function(spatial_df, 
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
    cross_K_df <- calculate_cross_K3D(spatial_df,
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
calculate_cross_K3D <- function(spatial_df, 
                                reference_cell_type, 
                                target_cell_types, 
                                radius, 
                                feature_colname = "Cell.Type") {
  
  if (is.null(spatial_df[[feature_colname]])) stop(paste("No column called", feature_colname, "found in spatial_df object"))
  
  if (is.null(spatial_df[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spatial_df")
    spatial_df$Cell.ID <- paste("Cell", seq(nrow(spatial_df)), sep = "_")
  }  
  
  
  ## Get expected cross K-function
  expected_cross_K <- (4/3) * pi * radius^3
  
  ## For reference_cell_type, check it is found in the spatial_df object
  if (!(reference_cell_type %in% spatial_df[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spatial_df object"))
    result <- data.frame(observed_cross_K = NA,
                         expected_cross_K = expected_cross_K,
                         cross_K_ratio = NA)
    return(result)
  }
  
  ## Get rough dimensions of the window the points are in
  spatial_df_coords <- spatial_df[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
  
  length <- round(max(spatial_df_coords$Cell.X.Position) - min(spatial_df_coords$Cell.X.Position))
  width  <- round(max(spatial_df_coords$Cell.Y.Position) - min(spatial_df_coords$Cell.Y.Position))
  height <- round(max(spatial_df_coords$Cell.Z.Position) - min(spatial_df_coords$Cell.Z.Position))
  ## Get volume of the window the cells are in
  volume <- length * width * height
  
  # Number of reference cell types is constant
  n_ref_cells <- sum(spatial_df[[feature_colname]] == reference_cell_type)
  
  # Define result data frame
  result <- data.frame(reference = reference_cell_type, expected = expected_cross_K)
  
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood3D(spatial_df,
                                                                  reference_cell_type,
                                                                  target_cell_types,
                                                                  radius,
                                                                  feature_colname,
                                                                  show_summary = FALSE,
                                                                  plot_image = FALSE)
  
  for (target_cell_type in target_cell_types) {
    
    n_ref_tar_interactions <- sum(cells_in_neighbourhood_df[[target_cell_type]])
    
    n_tar_cells <- sum(spatial_df[[feature_colname]] == target_cell_type)
    
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
calculate_cross_L_gradient3D <- function(spatial_df, 
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
    cross_L_df <- calculate_cross_L3D(spatial_df,
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
calculate_cross_L3D <- function(spatial_df, 
                                reference_cell_type, 
                                target_cell_types, 
                                radius, 
                                feature_colname = "Cell.Type") {
  
  result <- calculate_cross_K3D(spatial_df = spatial_df,
                                reference_cell_type = reference_cell_type,
                                target_cell_types = target_cell_types,
                                radius = radius,
                                feature_colname = feature_colname)
  
  result[ , c("expected", target_cell_types)] <- (result[ , c("expected", target_cell_types)] / (4 * pi / 3)) ^ (1/3)
  
  return(result)
}
calculate_entropy_background3D <- function(spatial_df,
                                           cell_types_of_interest, 
                                           feature_colname = "Cell.Type") {
  
  # NULL case: entropy is undefined
  if (is.null(cell_types_of_interest)) return(NA)
  
  # One cell type case: entropy is 0
  if (is.character(cell_types_of_interest) && length(cell_types_of_interest) == 1) return(0)
  
  cell_proportions_data <- calculate_cell_proportions3D(spatial_df, cell_types_of_interest, feature_colname, FALSE)
  
  # Calculate entropy of the entire image
  entropy <- -1 * sum(cell_proportions_data$proportion * log(cell_proportions_data$proportion, length(cell_proportions_data$proportion)))
  
  return(entropy) 
}
calculate_entropy_gradient3D <- function(spatial_df,
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
    entropy_df <- calculate_entropy3D(spatial_df,
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
    expected_entropy <- calculate_entropy_background3D(spatial_df, target_cell_types, feature_colname)
    fig <- plot_entropy_gradient3D(result, expected_entropy, reference_cell_type, target_cell_types)
    methods::show(fig)
  }
  
  return(result)
}
calculate_entropy_grid_metrics3D <- function(spatial_df, 
                                             n_splits,
                                             cell_types_of_interest,
                                             feature_colname = "Cell.Type",
                                             plot_image = TRUE) {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  ## Check cell_types_of_interest are found in the spatial_df object
  unknown_cell_types <- setdiff(cell_types_of_interest, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  ## If cell types have been chosen, check they are found in the spatial_df object
  unknown_cell_types <- setdiff(cell_types_of_interest, unique(spatial_df[[feature_colname]]))
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Add grid metrics to spatial_df
  grid_metrics <- get_spatial_df_grid_metrics3D(spatial_df, n_splits, feature_colname)
  
  # Get grid_prism_cell_matrix from spatial_df
  grid_prism_cell_matrix <- grid_metrics$grid_prism_cell_matrix
  
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
  result <- cbind(result, grid_metrics$grid_prism_coordinates)
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_metrics_continuous3D(result, "entropy")
    methods::show(fig)
  }
  
  return(result)
}
calculate_entropy3D <- function(spatial_df,
                                reference_cell_type,
                                target_cell_types,
                                radius,
                                feature_colname = "Cell.Type") {
  
  # Check target_cell_types
  if (!(is.character(target_cell_types) && length(target_cell_types) >= 2)) {
    stop("`target_cell_types` is not a character vector with at least 2 cell types.")
  }
  
  ## Users should ensure include the reference_cell_type as one of the target_cell_types
  cells_in_neighbourhood_proportion_df <- calculate_cells_in_neighbourhood_proportions3D(spatial_df,
                                                                                         reference_cell_type,
                                                                                         target_cell_types,
                                                                                         radius,
                                                                                         feature_colname)
  
  if (is.null(cells_in_neighbourhood_proportion_df)) return(NULL)
  
  ## Get entropy for target_cell_type
  cells_in_neighbourhood_proportion_df[ , paste(target_cell_types, "_entropy", sep = "")] <- 
    -1 * 
    (cells_in_neighbourhood_proportion_df[ , paste(target_cell_types, "_prop", sep = "")] * log(cells_in_neighbourhood_proportion_df[ , paste(target_cell_types, "_prop", sep = "")], 2) +
       (1 - cells_in_neighbourhood_proportion_df[ , paste(target_cell_types, "_prop", sep = "")]) * log(1 - (cells_in_neighbourhood_proportion_df[ , paste(target_cell_types, "_prop", sep = "")]), 2))
  
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
#'    It allows you to spatial_dfcify a subset of cell types to analyse and provides the option to summarise 
#'    the results and plot violin plots of the minimum distances between cell types.
#'
#' @param spatial_df A SpatialExperiment object containing 3D spatial information for the cells.
#' @param cell_types_of_interest A character vector spatial_dfcifying the cell types of interest.
#'   If NULL, all cell types in the `feature_colname` column will be considered.
#' @param feature_colname A string spatial_dfcifying the name of the column in the `colData` slot of the SpatialExperiment
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
#'     spatial_df = SPIAT3D::simulated_spatial_df,
#'     cell_types_of_interest = NULL,
#'     feature_colname = "Cell.Type",
#'     show_summary = TRUE,
#'     plot_image = TRUE
#' )
#' 
#' @export


calculate_minimum_distances_between_cell_types3D <- function(spatial_df,
                                                             cell_types_of_interest = NULL,
                                                             feature_colname = "Cell.Type",
                                                             show_summary = TRUE,
                                                             plot_image = TRUE) {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
  }
  if (ncol(spatial_df) < 2) {
    stop("There must be at least two cells in spatial_df.")
  }
  if (!(is.null(cell_types_of_interest) || is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(show_summary)) {
    stop("`show_summary` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  if (is.null(spatial_df[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spatial_df")
    spatial_df$Cell.ID <- paste("Cell", seq(nrow(spatial_df)), sep = "_")
  }  
  
  # De-factor feature column in spatial_df object
  spatial_df[[feature_colname]] <- as.character(spatial_df[[feature_colname]])
  
  # Subset spatial_df to only contain the cells of interest
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spatial_df object
    unknown_cell_types <- setdiff(cell_types_of_interest, spatial_df[[feature_colname]])
    if (length(unknown_cell_types) != 0) {
      warning(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
                    paste(unknown_cell_types, collapse = ", ")))
    }
    
    spatial_df <- spatial_df[spatial_df[[feature_colname]] %in% cell_types_of_interest, ]
  }
  # If cell_types_of_interest is NULL, use all cells in spatial_df
  else {
    cell_types_of_interest <- unique(spatial_df[[feature_colname]])
  }
  
  # Create a list containing the cell IDs of each cell type
  cell_type_ids <- list()
  for (cell_type in cell_types_of_interest) {
    cell_type_ids[[cell_type]] <- as.character(spatial_df$Cell.ID[spatial_df[[feature_colname]] == cell_type])
  }
  
  # Get spatial_df coords
  spatial_df_coords <- spatial_df[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
  
  # Get different possible cell type combinations
  # Each row represents a combination
  # If a row is [1 , 2], then we are comparing cell type 1 and cell type 2
  permu <- gtools::permutations(length(cell_types_of_interest), 2, repeats.allowed = TRUE)
  
  result <- data.frame()
  
  for (i in seq(nrow(permu))) {
    cell_type1 <- cell_types_of_interest[permu[i, 1]]
    cell_type2 <- cell_types_of_interest[permu[i, 2]]
    
    # Don't have one of the cells
    if (sum(spatial_df[[feature_colname]] == cell_type1) == 0 || sum(spatial_df[[feature_colname]] == cell_type2) == 0) {
      result <- rbind(result, data.frame(ref_cell_id = NA, ref_cell_type = cell_type1, nearest_cell_id = NA, nearest_cell_type = cell_type2, distance = NA))
      next
    }
    
    # Get x, y, z coords for all cells of cell_type1 and cell_type2
    cell_type1_coords <- spatial_df_coords[spatial_df[[feature_colname]] == cell_type1, ]
    cell_type2_coords <- spatial_df_coords[spatial_df[[feature_colname]] == cell_type2, ]
    
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
calculate_minimum_distances_to_clusters3D <- function(spatial_df, 
                                                      cell_types_inside_cluster, 
                                                      cell_types_outside_cluster, 
                                                      cluster_colname, 
                                                      feature_colname = "Cell.Type", 
                                                      plot_image = T) {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
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
  if (is.null(spatial_df[[cluster_colname]])) {
    stop(paste("No column called", cluster_colname, "found in spatial_df object."))
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  ## Add Cell.ID column
  if (is.null(spatial_df[["Cell.ID"]])) {
    warning("Temporarily adding Cell.Id column to your spatial_df")
    spatial_df$Cell.ID <- paste("Cell", seq(nrow(spatial_df)), sep = "_")
  }
  
  ## For each cell type outside clusters, get their set of coords. These exclude cell types in clusters
  spatial_df_coords <- spatial_df[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
  
  # Cells outside cluster have a cluster number of 0 (i.e. they are not in a cluster)
  spatial_df_outside_cluster <- spatial_df[ , spatial_df[[cluster_colname]] == 0]
  
  cell_types_outside_cluster_coords <- list()
  for (cell_type in cell_types_outside_cluster) {
    cell_types_outside_cluster_coords[[cell_type]] <- spatialCoords(spatial_df_outside_cluster)[spatial_df_outside_cluster[[feature_colname]] == cell_type, ]
  }
  
  ## For each cluster, determine the minimum distance of each outside_cell_type  
  result <- vector()
  
  # Get number of clusters
  n_clusters <- max(spatial_df[[cluster_colname]])
  
  for (i in seq(n_clusters)) {
    cluster_coords <- spatial_df_coords[spatial_df[[cluster_colname]] == i & spatial_df[[feature_colname]] %in% cell_types_inside_cluster, ]
    cluster_cell_types <- spatial_df[["Cell.Type"]][spatial_df[[cluster_colname]] == i & spatial_df[[feature_colname]] %in% cell_types_inside_cluster]
    cluster_cell_ids <- spatial_df[["Cell.ID"]][spatial_df[[cluster_colname]] == i & spatial_df[[feature_colname]] %in% cell_types_inside_cluster]
    
    for (outside_cell_type in cell_types_outside_cluster) {
      curr_cell_type_coords <- cell_types_outside_cluster_coords[[outside_cell_type]]
      
      all_closest <- RANN::nn2(data = cluster_coords, 
                               query = curr_cell_type_coords, 
                               k = 1) 
      
      local_dist_mins <- data.frame(
        cluster_number = i,
        outside_cell_id = as.character(spatial_df_outside_cluster$Cell.ID[spatial_df_outside_cluster[["Cell.Type"]] == outside_cell_type]),
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
calculate_mixing_scores_gradient3D <- function(spatial_df, 
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
    mixing_scores <- calculate_mixing_scores3D(spatial_df,
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
calculate_mixing_scores3D <- function(spatial_df, 
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
      n_ref <- sum(spatial_df[[feature_colname]] == reference_cell_type)
      n_tar <- sum(spatial_df[[feature_colname]] == target_cell_type)
      
      
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
      cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood3D(spatial_df,
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
          methods::show(paste("There are no reference to reference interactions for", target_cell_type, "in the spatial_dfcified radius, cannot calculate mixing score"))
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
calculate_pairwise_distances_between_cell_types3D <- function(spatial_df,
                                                              cell_types_of_interest = NULL,
                                                              feature_colname = "Cell.Type",
                                                              show_summary = TRUE,
                                                              plot_image = TRUE) {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
  }
  if (ncol(spatial_df) < 2) {
    stop("There must be at least two cells in spatial_df.")
  }
  if (!(is.null(cell_types_of_interest) || is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(show_summary)) {
    stop("`show_summary` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  if (is.null(spatial_df[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spatial_df")
    spatial_df$Cell.ID <- paste("Cell", seq(nrow(spatial_df)), sep = "_")
  }
  
  # De-factor feature column in spatial_df object
  spatial_df[[feature_colname]] <- as.character(spatial_df[[feature_colname]])
  
  # Subset spatial_df to only contain the cells of interest
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spatial_df object
    unknown_cell_types <- setdiff(cell_types_of_interest, spatial_df[[feature_colname]])
    if (length(unknown_cell_types) != 0) {
      warning(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
                    paste(unknown_cell_types, collapse = ", ")))
    }
    
    spatial_df <- spatial_df[spatial_df[[feature_colname]] %in% cell_types_of_interest, ]
  }
  # If cell_types_of_interest is NULL, use all cells in spatial_df
  else {
    cell_types_of_interest <- unique(spatial_df[[feature_colname]])
  }
  
  # Create a list containing the cell IDs of each cell type
  cell_type_ids <- list()
  for (cell_type in cell_types_of_interest) {
    cell_type_ids[[cell_type]] <- as.character(spatial_df$Cell.ID[spatial_df[[feature_colname]] == cell_type])
  }
  
  # Calculate cell to cell distances
  distance_matrix <- -1 * apcluster::negDistMat(spatial_df[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")])
  rownames(distance_matrix) <- spatial_df$Cell.ID
  colnames(distance_matrix) <- spatial_df$Cell.ID
  
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
  ## Grid prisms within this spatial_dfcified threshold have a weight of 1, otherwise, weight of 0
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
calculate_volume_of_clusters3D <- function(spatial_df, cluster_colname) {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
  }
  if (!is.character(cluster_colname)) {
    stop("`cluster_colname` is not a character. This should be 'alpha_hull_cluster', 'dbscan_cluster', or 'grid_based_cluster', depending on the chosen method.")
  }
  if (is.null(spatial_df[[cluster_colname]])) {
    stop(paste("No column called", cluster_colname, "found in spatial_df object."))
  }
  
  # Get number of clusters
  n_clusters <- max(spatial_df[[cluster_colname]])
  
  ### 1. Estimate volume of each cluster by density of the window. ------------
  
  ## For each cluster, determine the number of cells in each cluster of each cluster
  result <- data.frame(matrix(nrow = n_clusters, ncol = 2))
  colnames(result) <- c("cluster_number", "n_cells")
  
  for (i in seq(n_clusters)) {
    result[i, "n_cells"] <- sum(spatial_df[[cluster_colname]] == i)
  }
  result$cluster_number <- as.character(seq(n_clusters))
  
  ## Assume window is a rectangular prism
  spatial_df_coords <- spatial_df[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
  
  length <- round(max(spatial_df_coords$Cell.X.Position) - min(spatial_df_coords$Cell.X.Position))
  width  <- round(max(spatial_df_coords$Cell.Y.Position) - min(spatial_df_coords$Cell.Y.Position))
  height <- round(max(spatial_df_coords$Cell.Z.Position) - min(spatial_df_coords$Cell.Z.Position))
  
  window_volume <- length * width * height
  
  result$volume_by_density <- (result$n_cells / ncol(spatial_df)) * window_volume
  
  
  ### 2. If cluster_colname == "alpha_hull_cluster", use the volume method found in the alphashape3d package
  if (cluster_colname == "alpha_hull_cluster") {
    result$volume_by_alpha_hull <- volume_ashape3d(spatial_df@metadata$alpha_hull$ashape3d_object, byComponents = T)
  }
  
  
  ### 3. If cluster_colname == "grid_based_cluster", sum the volume of each grid prism to get volume of each cluster
  if (cluster_colname == "grid_based_cluster") {
    result$volume_by_grid <- 0
    i <- 1
    for (grid_cluster in spatial_df@metadata$grid_prisms) {
      result[i, "volume_by_grid"] <- sum(grid_cluster$l * grid_cluster$w * grid_cluster$h)
      i <- i + 1
    }
  }
  
  return(result)
}
library(dbscan)

dbscan_clustering3D <- function(spatial_df,
                                cell_types_of_interest,
                                radius,
                                minimum_cells_in_radius,
                                minimum_cells_in_cluster,
                                feature_colname = "Cell.Type",
                                plot_image = T) {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
  }
  ## Check cell types of interst are found in the spatial_df object
  unknown_cell_types <- setdiff(cell_types_of_interest, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
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
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  spatial_df_subset <- spatial_df[ , spatial_df[[feature_colname]] %in% cell_types_of_interest]
  spatial_df_subset_coords <- spatialCoords(spatial_df_subset)
  
  db <- dbscan::dbscan(spatial_df_subset_coords, eps = radius, minPts = minimum_cells_in_radius, borderPoints = F)
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
  
  ## Convert spatial_df object to data frame
  df <- data.frame(spatial_df[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")], colData(spatial_df))
  
  df_cell_types_of_interest <- df[df[[feature_colname]] %in% cell_types_of_interest, ]
  df_other_cell_types <- df[!(df[[feature_colname]] %in% cell_types_of_interest), ]
  
  df_cell_types_of_interest$dbscan_cluster <- db$cluster
  df_other_cell_types$dbscan_cluster <- 0
  
  ## Convert data frame to spatial_df object
  df <- rbind(df_cell_types_of_interest, df_other_cell_types)
  
  spatial_df <- SpatialExperiment(
    assay = matrix(data = NA, nrow = nrow(df), ncol = nrow(df)),
    colData = df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = spatial_df@metadata)
  
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
  
  return(spatial_df)
}
get_spatial_df_grid_metrics3D <- function(spatial_df, 
                                          n_splits, 
                                          feature_colname = "Cell.Type") {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  
  spatial_df_coords <- spatial_df[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
  
  ## Get dimensions of the window
  min_x <- min(spatial_df_coords[ , "Cell.X.Position"])
  min_y <- min(spatial_df_coords[ , "Cell.Y.Position"])
  min_z <- min(spatial_df_coords[ , "Cell.Z.Position"])
  
  max_x <- max(spatial_df_coords[ , "Cell.X.Position"])
  max_y <- max(spatial_df_coords[ , "Cell.Y.Position"])
  max_z <- max(spatial_df_coords[ , "Cell.Z.Position"])
  
  length <- round(max_x - min_x)
  width  <- round(max_y - min_y)
  height <- round(max_z - min_z)
  
  ## Get distance of row, col and lay
  d_row <- length / n_splits
  d_col <- width / n_splits
  d_lay <- height / n_splits
  
  # Shift spatial_df_coords so they begin at the origin
  spatial_df_coords[, "Cell.X.Position"] <- spatial_df_coords[, "Cell.X.Position"] - min_x
  spatial_df_coords[, "Cell.Y.Position"] <- spatial_df_coords[, "Cell.Y.Position"] - min_y
  spatial_df_coords[, "Cell.Z.Position"] <- spatial_df_coords[, "Cell.Z.Position"] - min_z
  
  ## Figure out which 'grid prism number' each cell is inside
  spatial_df$grid_prism_num <- floor(spatial_df_coords[ , "Cell.X.Position"] / d_row) +
    floor(spatial_df_coords[ , "Cell.Y.Position"] / d_col) * n_splits + 
    floor(spatial_df_coords[ , "Cell.Z.Position"] / d_lay) * n_splits^2 + 1
  
  ## Determine the cell types found in each grid prism
  n_grid_prisms <- n_splits^3
  grid_prism_cell_matrix <- as.data.frame.matrix(table(spatial_df[[feature_colname]], factor(spatial_df$grid_prism_num, levels = seq(n_grid_prisms))))
  
  grid_prism_cell_matrix <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                       t(grid_prism_cell_matrix), check.names = FALSE)
  
  ## Determine centre coordinates of each grid prism
  grid_prism_coordinates <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                       x_coord = ((seq(n_grid_prisms) - 1) %% n_splits + 0.5) * d_row + round(min_x),
                                       y_coord = (floor(((seq(n_grid_prisms) - 1) %% (n_splits)^2) / n_splits) + 0.5) * d_col + round(min_y),
                                       z_coord = (floor((seq(n_grid_prisms) - 1) / (n_splits^2)) + 0.5) * d_lay + round(min_z))
  
  grid_metrics <- list("grid_prism_cell_matrix" = grid_prism_cell_matrix,
                       "grid_prism_coordinates" = grid_prism_coordinates)
  
  return(grid_metrics)
}
grid_based_cluster_recursion3D <- function(df,  # Using a df is much faster than using a spatial_df
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
  
  # Get cell types from spatial_df grid prism
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
grid_based_clustering3D <- function(spatial_df,
                                    cell_types_of_interest,
                                    n_splits,
                                    minimum_cells_in_cluster,
                                    feature_colname = "Cell.Type",
                                    plot_image = TRUE) {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
  }
  ## Check cell types of interst are found in the spatial_df object
  unknown_cell_types <- setdiff(cell_types_of_interest, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
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
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Add grid metrics to spatial_df
  spatial_df <- get_spatial_df_grid_metrics3D(spatial_df, n_splits, feature_colname)
  
  # Get grid_prism_cell_matrix from spatial_df
  grid_prism_cell_matrix <- spatial_df@metadata$grid_metrics$grid_prism_cell_matrix
  
  ## Calculate proportions for each grid prism
  if (length(cell_types_of_interest) == 1) {
    grid_prism_cell_proportions <- grid_prism_cell_matrix[ , cell_types_of_interest]
  }
  else {
    grid_prism_cell_proportions <- rowSums(grid_prism_cell_matrix[ , cell_types_of_interest])
  }
  grid_prism_cell_proportions <- grid_prism_cell_proportions / rowSums(grid_prism_cell_matrix[ , unique(spatial_df[[feature_colname]])])
  n_grid_prisms <- n_splits^3
  names(grid_prism_cell_proportions) <- seq(n_grid_prisms)
  
  
  ## Create template for final result
  result <- list()
  n_clusters <- 1
  
  ## Get dimensions of the window
  spatial_df_coords <- spatial_df[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
  
  min_x <- min(spatial_df_coords$Cell.X.Position)
  min_y <- min(spatial_df_coords$Cell.Y.Position)
  min_z <- min(spatial_df_coords$Cell.Z.Position)
  
  max_x <- max(spatial_df_coords$Cell.X.Position)
  max_y <- max(spatial_df_coords$Cell.Y.Position)
  max_z <- max(spatial_df_coords$Cell.Z.Position)
  
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
    df <- spatial_df_coords
    df[[feature_colname]] <- spatial_df[[feature_colname]] 
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
  
  ## Add grid_based_cluster column to spatial_df, indicating which cluster each cell belongs to
  spatial_df$grid_based_cluster <- 0
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
      
      spatial_df$grid_based_cluster <- ifelse(spatial_df_coords$Cell.X.Position >= x &
                                                spatial_df_coords$Cell.X.Position < (x + l) &
                                                spatial_df_coords$Cell.Y.Position >= y &
                                                spatial_df_coords$Cell.Y.Position < (y + w) &
                                                spatial_df_coords$Cell.Z.Position >= z &
                                                spatial_df_coords$Cell.Z.Position < (z + h) &
                                                spatial_df[[feature_colname]] %in% cell_types_of_interest, 
                                              cluster_number, 
                                              spatial_df$grid_based_cluster)
      
    }
    # Check if current cluster surpasses the minimum_cells_in_cluster threshold
    if (sum(spatial_df$grid_based_cluster == cluster_number) < minimum_cells_in_cluster) {
      spatial_df$grid_based_cluster[spatial_df$grid_based_cluster == cluster_number] <- 0
      result[[paste("cluster", i, sep = "_")]] <- NULL
      
    }
    else {
      cluster_number <- cluster_number + 1 
    }
  }
  
  n_clusters <- max(spatial_df$grid_based_cluster)
  if (n_clusters == 0) {
    stop("All clusters identified do not meet the `minimum_cells_in_cluster` threshold. Consider lowering the `minimum_cells_in_cluster` parameter.")
  }
  
  # re-name each grid_based cluster
  names(result) <- paste("cluster", seq_len(length(result)), sep = "_")
  
  # Add grid_clustering result to spatial_df metadata
  spatial_df@metadata[["grid_prisms"]] <- result
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_based_clusters3D(spatial_df, feature_colname = feature_colname)
    methods::show(fig)
  }
  
  return(spatial_df)
}
plot_alpha_hull_clusters3D <- function(spatial_df_with_alpha_hull, 
                                       plot_cell_types = NULL,
                                       plot_colours = NULL,
                                       feature_colname = "Cell.Type") {
  
  # Check input parameters
  if (class(spatial_df_with_alpha_hull) != "SpatialExperiment") {
    stop("`spatial_df_with_alpha_hull` is not a SpatialExperiment object.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df_with_alpha_hull[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  
  ## If no cell types chosen, use all cell types found in data frame
  if (is.null(plot_cell_types)) plot_cell_types <- unique(spatial_df_with_alpha_hull[[feature_colname]])
  
  ## If cell types have been chosen, check they are found in the spatial_df object
  unknown_cell_types <- setdiff(plot_cell_types, spatial_df_with_alpha_hull[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following plot_cell_types are not found in the spatial_df object:\n   ",
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
  
  ## Convert spatial_df object to data frame
  df <- data.frame(spatialCoords(spatial_df_with_alpha_hull), "Cell.Type" = spatial_df_with_alpha_hull[[feature_colname]])
  
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
  alpha_hull_clusters <- spatial_df_with_alpha_hull$alpha_hull_cluster[spatial_df_with_alpha_hull$alpha_hull_cluster != 0]
  
  # Get number of alpha hulls
  n_alpha_hulls <- length(unique(alpha_hull_clusters))
  
  vertices <- spatial_df_with_alpha_hull@metadata$alpha_hull$vertices
  faces <- data.frame(spatial_df_with_alpha_hull@metadata$alpha_hull$faces)
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
plot_cells3D <- function(spatial_df,
                         plot_cell_types = NULL,
                         plot_colours = NULL,
                         feature_colname = "Cell.Type",
                         aspectmode = "data") {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
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
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste(feature_colname, "is not a valid column in your spatial_df object."))
  }
  
  ## Convert spatial_df object to data frame
  df <- data.frame(spatial_df[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")], "Cell.Type" = spatial_df[[feature_colname]])
  
  ## If no cell types chosen, use all cell types found in data frame
  if (is.null(plot_cell_types)) {
    warning("plot_cell_types not spatial_dfcified, all cell types found in the spatial_df object will be used.")
    plot_cell_types <- unique(df[["Cell.Type"]])
  }
  ## If no colours inputted, use rainbow palette
  if (is.null(plot_colours)) {
    warning("plot_colours not spatial_dfcified, rainbow palette will be used.")
    plot_colours <- rainbow(length(plot_cell_types))
  }
  ## User inputs mismatching cell types and colours
  if (length(plot_cell_types) != length(plot_colours)) {
    stop("Length of plot_cell_types is not equal to length of plot_colours")
  }
  
  ## If cell types have been chosen, check they are found in the spatial_df object
  spatial_df_cell_types <- unique(spatial_df[[feature_colname]])
  unknown_cell_types <- setdiff(plot_cell_types, spatial_df_cell_types)
  
  if (length(unknown_cell_types) == length(plot_cell_types)) {
    stop("None of the plot_cell_types are found in the spatial_df object")
  }
  
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following plot_cell_types are not found in the spatial_df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    plot_colours <- plot_colours[which(plot_cell_types %in% spatial_df_cell_types)]
    plot_cell_types <- intersect(plot_cell_types, spatial_df_cell_types)
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
                                                  titlefont = list(size = 20), tickfont = list(size = 15)),
                                     aspectmode = aspectmode))
  
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
    scale_colour_manual(values = c("expected" = "black", setNames(RColorBrewer::brewer.pal(length(target_cell_types), "Set3"), target_cell_types)), 
                        name = "") +
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
    scale_colour_manual(values = c("expected" = "black", setNames(RColorBrewer::brewer.pal(length(target_cell_types), "Set3"), target_cell_types)), 
                        name = "") +
    theme_bw()
  
  return(fig) 
}
plot_cross_L_gradient3D <- function(cross_L_gradient_df) {
  
  target_cell_types <- colnames( cross_L_gradient_df)[!colnames(cross_L_gradient_df) %in% c("reference", "expected", "radius")]
  
  plot_result <- reshape2::melt(cross_L_gradient_df, "radius", c(target_cell_types, "expected"))
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Cross L-function gradient", x = "Radius", y = "Cross L-function value") +
    scale_colour_manual(values = c("expected" = "black", setNames(RColorBrewer::brewer.pal(length(target_cell_types), "Set3"), target_cell_types)), 
                        name = "") +
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
    fig <- fig + labs(subtitle = paste("Reference: ", reference_cell_type, ", Target: ", paste(target_cell_types, collapse = ", "), sep = ""))
  }
  
  return(fig)
}
plot_grid_based_clusters3D <- function(spatial_df_with_grid, 
                                       plot_cell_types = NULL,
                                       plot_colours = NULL,
                                       feature_colname = "Cell.Type") {
  
  # Check input parameters
  if (class(spatial_df_with_grid) != "SpatialExperiment") {
    stop("`spatial_df_with_grid` is not a SpatialExperiment object.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df_with_grid[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  
  ## If no cell types chosen, use all cell types found in data frame
  if (is.null(plot_cell_types)) plot_cell_types <- unique(spatial_df_with_grid[[feature_colname]])
  
  ## If cell types have been chosen, check they are found in the spatial_df object
  unknown_cell_types <- setdiff(plot_cell_types, spatial_df_with_grid[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following plot_cell_types are not found in the spatial_df object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  
  ## If no colours inputted, use rainbow palette
  if (is.null(plot_colours)) plot_colours <- rainbow(length(plot_cell_types))
  
  ## User inputs mismatching cell types and colours
  if (length(plot_cell_types) != length(plot_colours)) stop("Length of plot_cell_types is not equal to length of plot_colours")
  
  ## Convert spatial_df object to data frame
  df <- data.frame(spatialCoords(spatial_df_with_grid), colData(spatial_df_with_grid))
  
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
  n_grid_based_clusters <- length(spatial_df_with_grid@metadata[["grid_prisms"]])
  
  faces <- data.frame(edge1 = c(1, 1, 1, 1, 1, 1, 8, 8, 8, 8, 8, 8),
                      edge2 = c(2, 5, 2, 3, 3, 5, 6, 4 ,7, 6, 7, 4),
                      edge3 = c(6, 6, 4, 4, 7, 7, 2, 2, 5, 5, 3, 3))
  grid_based_colours <- rainbow(n_grid_based_clusters)
  
  ## Add grid-based clusters to fig, one by one  
  for (i in seq(n_grid_based_clusters)) {
    
    grid_based_cluster <- spatial_df_with_grid@metadata[["grid_prisms"]][[i]]
    
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




# SPIAT-2D functions -----
library(SpatialExperiment)
library(dbscan)

library(apcluster)
library(plotly)
library(dplyr)
library(reshape2)
library(gtools)
library(cowplot)
library(Hmisc)


calculate_all_gradient_cc_metrics2D <- function(spatial_df, 
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
                 "entropy" = data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types))),
                 "cross_K" = data.frame(matrix(nrow = length(radii), ncol = length(cross_K_df_colnames))),
                 "cross_L" = data.frame(matrix(nrow = length(radii), ncol = length(cross_K_df_colnames))),
                 "cross_G" = list(),
                 "co_occurrence" = data.frame(matrix(nrow = length(radii), ncol = length(co_occurrence_df_colnames))))
  colnames(result[["cells_in_neighbourhood"]]) <- target_cell_types
  colnames(result[["cells_in_neighbourhood_proportion"]]) <- target_cell_types
  colnames(result[["entropy"]]) <- target_cell_types
  colnames(result[["cross_K"]]) <- cross_K_df_colnames
  colnames(result[["cross_L"]]) <- cross_K_df_colnames
  colnames(result[["co_occurrence"]]) <- co_occurrence_df_colnames
  
  # Define indiviudal data frames for mixing_score and cross_G
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
    df <- calculate_all_single_radius_cc_metrics2D(spatial_df,
                                                   reference_cell_type,
                                                   target_cell_types,
                                                   radii[i],
                                                   feature_colname)
    
    if (is.null(df)) return(NULL)
    
    df[["cells_in_neighbourhood"]]$ref_cell_id <- NULL
    
    result[["cells_in_neighbourhood"]][i, ] <- apply(df[["cells_in_neighbourhood"]], 2, mean)
    result[["cells_in_neighbourhood_proportion"]][i, ] <- apply(df[["cells_in_neighbourhood_proportion"]][ , paste(target_cell_types, "_prop", sep = "")], 2, mean, na.rm = T)
    result[["entropy"]][i, ] <- apply(df[["entropy"]][ , paste(target_cell_types, "_entropy", sep = "")], 2, mean, na.rm = T)
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
    
    expected_entropy <- calculate_entropy_background2D(spatial_df, target_cell_types, feature_colname)
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

calculate_all_single_radius_cc_metrics2D <- function(spatial_df, 
                                                     reference_cell_type, 
                                                     target_cell_types, 
                                                     radius, 
                                                     feature_colname = "Cell.Type") {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
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
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  
  ## For reference_cell_type, check it is found in the spatial_df object
  if (!(reference_cell_type %in% spatial_df[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spatial_df object"))
    return(NULL)
  }
  ## For target_cell_types, check they are found in the spatial_df object
  unknown_cell_types <- setdiff(target_cell_types, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spatial_df object:\n   ",
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
  spatial_df_coords <- spatial_df[ , c("Cell.X.Position", "Cell.Y.Position")]
  length <- round(max(spatial_df_coords$Cell.X.Position) - min(spatial_df_coords$Cell.X.Position))
  width  <- round(max(spatial_df_coords$Cell.Y.Position) - min(spatial_df_coords$Cell.Y.Position))
  
  ## Get volume of the window the cells are in
  volume <- length * width
  
  # All single radius cc metrics stem from calculate_entropy2D function
  entropy_df <- calculate_entropy2D(spatial_df, 
                                    reference_cell_type, 
                                    target_cell_types, 
                                    radius, 
                                    feature_colname)  
  
  ## Cells in neighbourhood ----------
  result[["cells_in_neighbourhood"]] <- entropy_df[ , c("ref_cell_id", target_cell_types)]
  
  ## Cells in neighbourhood proportion ----------
  result[["cells_in_neighbourhood_proportion"]] <- entropy_df[ , c("ref_cell_id", paste(target_cell_types, "_prop", sep = ""))]
  
  ## Entropy --------------
  result[["entropy"]] <- entropy_df[ , c("ref_cell_id", paste(target_cell_types, "_entropy", sep = ""))]
  
  ## Mixing score -----------------
  for (target_cell_type in target_cell_types) {
    mixing_score_df <- data.frame(matrix(nrow = 1, ncol = length(mixing_score_df_colnames)))
    colnames(mixing_score_df) <- mixing_score_df_colnames
    mixing_score_df$ref_cell_type <- reference_cell_type
    
    # No need to fill in mixing_score_df if the reference and target cell is the same
    if (reference_cell_type != target_cell_type) {
      mixing_score_df$tar_cell_type <- target_cell_type
      mixing_score_df$n_ref_cells <- sum(spatial_df[[feature_colname]] == reference_cell_type)
      mixing_score_df$n_tar_cells <- sum(spatial_df[[feature_colname]] == target_cell_type)
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
  cross_K_df$expected <- pi * radius^2
  
  for (target_cell_type in target_cell_types) {
    cross_K_df[[target_cell_type]] <- (((volume * sum(entropy_df[[target_cell_type]])) / sum(spatial_df[[feature_colname]] == reference_cell_type)) / sum(spatial_df[[feature_colname]] == target_cell_type)) 
  }
  result[["cross_K"]] <- cross_K_df
  
  ## Cross_L ---------------------
  cross_L_df <- cross_K_df
  cross_L_df[ , c("expected", target_cell_types)] <- (cross_L_df[ , c("expected", target_cell_types)] / (pi)) ^ (1/2)
  result[["cross_L"]] <- cross_L_df
  
  ## Cross_G ---------------------
  for (target_cell_type in target_cell_types) {
    cross_G_df <- data.frame(matrix(nrow = 1, ncol = length(cross_G_df_colnames)))
    colnames(cross_G_df) <- cross_G_df_colnames
    
    reference_target_interactions <- entropy_df[[target_cell_type]]
    n_target_cells <- sum(spatial_df[[feature_colname]] == target_cell_type)
    target_cell_type_intensity <- n_target_cells / volume
    observed_cross_G <- sum(reference_target_interactions != 0) / length(reference_target_interactions)
    expected_cross_G <- 1 - exp(-1 * target_cell_type_intensity * pi * radius^2)
    
    cross_G_df$observed_cross_G <- observed_cross_G
    cross_G_df$expected_cross_G <- expected_cross_G
    result[["cross_G"]][[target_cell_type]] <- cross_G_df
  } 
  
  ## Co_occurrence ---------------
  all_cell_types <- unique(spatial_df[[feature_colname]])
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood2D(spatial_df,
                                                                  reference_cell_type,
                                                                  all_cell_types,
                                                                  radius,
                                                                  feature_colname,
                                                                  F,
                                                                  F)
  
  cells_in_neighbourhood_df$total <- rowSums(cells_in_neighbourhood_df[, -1], na.rm = TRUE)
  
  
  co_occurrence_df <- data.frame(matrix(nrow = 1, ncol = length(co_occurrence_df_colnames)))
  colnames(co_occurrence_df) <- co_occurrence_df_colnames
  co_occurrence_df$reference <- reference_cell_type
  
  n_cells_in_spe <- length(spatial_df[[feature_colname]])
  n_cells_in_reference_cell_type_radius <- sum(cells_in_neighbourhood_df$total)
  
  for (target_cell_type in target_cell_types) {
    n_target_cells_in_reference_cell_type_radius <- sum(cells_in_neighbourhood_df[[target_cell_type]])
    target_cell_type_proportion_in_reference_cell_type_radius <- n_target_cells_in_reference_cell_type_radius / n_cells_in_reference_cell_type_radius
    n_target_cells_in_spe <- sum(spatial_df[[feature_colname]] == target_cell_type)
    target_cell_type_proportion_in_spe <- n_target_cells_in_spe / n_cells_in_spe
    target_cell_type_co_occurrence <- target_cell_type_proportion_in_reference_cell_type_radius / target_cell_type_proportion_in_spe
    
    co_occurrence_df[[target_cell_type]] <- target_cell_type_co_occurrence
  }
  result[["co_occurrence"]] <- co_occurrence_df
  
  return(result)
}


calculate_cell_proportion_grid_metrics2D <- function(spatial_df, 
                                                     n_splits,
                                                     reference_cell_types,
                                                     target_cell_types,
                                                     feature_colname = "Cell.Type",
                                                     plot_image = TRUE) {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  ## Check reference_cell_types are found in the spatial_df object
  unknown_cell_types <- setdiff(reference_cell_types, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in reference_cell_types are not found in the spatial_df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  ## Check target_cell_types are found in the spatial_df object
  unknown_cell_types <- setdiff(target_cell_types, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spatial_df object:\n   ",
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
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Add grid metrics to spatial_df
  grid_metrics <- get_spatial_df_grid_metrics2D(spatial_df, n_splits, feature_colname)
  
  # Get grid_prism_cell_matrix from spatial_df
  grid_prism_cell_matrix <- grid_metrics$grid_prism_cell_matrix
  
  ## Define data frame which contains all results
  n_grid_prisms <- n_splits^2
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
  result <- cbind(result, grid_metrics$grid_prism_coordinates)
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_metrics_continuous2D(result, "proportion")
    methods::show(fig)
  }
  
  return(result)
}


calculate_cell_proportions2D <- function(spatial_df,
                                         cell_types_of_interest = NULL, 
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
  }
  if (ncol(spatial_df) == 0) {
    stop("No cells found for calculating cell proportions.")
  }
  if (!(is.null(cell_types_of_interest) || is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Creates frequency/bar plot of all cell types in the entire image
  cell_proportions <- data.frame(table(spatial_df[[feature_colname]]))
  names(cell_proportions) <- c("cell_type", 'frequency')
  
  # Only include cell types the user has chosen
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spatial_df object
    unknown_cell_types <- setdiff(cell_types_of_interest, cell_proportions$cell_type)
    if (length(unknown_cell_types) != 0) {
      stop(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
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
calculate_cells_in_neighbourhood_gradient2D <- function(spatial_df, 
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
    cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood2D(spatial_df,
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
calculate_cells_in_neighbourhood_proportions_gradient2D <- function(spatial_df, 
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
    cell_proportions_neighbourhood_proportions_df <- calculate_cells_in_neighbourhood_proportions2D(spatial_df,
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
calculate_cells_in_neighbourhood_proportions2D <- function(spatial_df, 
                                                           reference_cell_type, 
                                                           target_cell_types, 
                                                           radius, 
                                                           feature_colname = "Cell.Type") {
  
  ## Get cells in neighbourhood df
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood2D(spatial_df,
                                                                  reference_cell_type,
                                                                  c(reference_cell_type, target_cell_types),
                                                                  radius,
                                                                  feature_colname,
                                                                  FALSE,
                                                                  FALSE)
  
  if (is.null(cells_in_neighbourhood_df)) return(NULL)
  
  cells_in_neighbourhood_df[ , paste(target_cell_types, "_prop", sep = "")] <- 
    cells_in_neighbourhood_df[ , target_cell_types] / (cells_in_neighbourhood_df[ , target_cell_types] + cells_in_neighbourhood_df[ , reference_cell_type])
  
  # If reference cell type is in target cell types, proportion should be 1
  if (reference_cell_type %in% target_cell_types) {
    cells_in_neighbourhood_df[cells_in_neighbourhood_df[[reference_cell_type]] != 0, paste(reference_cell_type, "_prop", sep = "")] <- 1
  }
  
  return(cells_in_neighbourhood_df)
}
calculate_cells_in_neighbourhood2D <- function(spatial_df, 
                                               reference_cell_type, 
                                               target_cell_types, 
                                               radius, 
                                               feature_colname = "Cell.Type",
                                               show_summary = TRUE,
                                               plot_image = TRUE) {
  
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
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
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(show_summary)) {
    stop("`show_summary` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  ## For reference_cell_type, check it is found in the spatial_df object
  if (!(reference_cell_type %in% spatial_df[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spatial_df object"))
    return(NULL)
  }
  ## For target_cell_types, check they are found in the spatial_df object
  unknown_cell_types <- setdiff(target_cell_types, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spatial_df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
  }
  
  if (is.null(spatial_df[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spatial_df")
    spatial_df$Cell.ID <- paste("Cell", seq(nrow(spatial_df)), sep = "_")
  }  
  
  # Get spatial_df coords
  spatial_df_coords <- spatial_df[ , c("Cell.X.Position", "Cell.Y.Position")]
  
  # Get reference_cell_type coords
  reference_cell_type_coords <- spatial_df_coords[spatial_df[[feature_colname]] == reference_cell_type, ]
  
  result <- data.frame(ref_cell_id = spatial_df$Cell.ID[spatial_df[[feature_colname]] == reference_cell_type])
  
  for (target_cell_type in target_cell_types) {
    
    if (sum(spatial_df[[feature_colname]] == target_cell_type) == 0) {
      result[[target_cell_type]] <- NA
      next
    }
    
    ## Get target_cell_type coords
    target_cell_type_coords <- spatial_df_coords[spatial_df[[feature_colname]] == target_cell_type, ]
    
    ## Determine number of target cells spatial_dfcified distance for each reference cell
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

calculate_co_occurrence_gradient2D <- function(spatial_df, 
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
    co_occurrence_df <- calculate_co_occurrence2D(spatial_df,
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
calculate_co_occurrence2D <- function(spatial_df, 
                                      reference_cell_type, 
                                      target_cell_types, 
                                      radius, 
                                      feature_colname = "Cell.Type") {
  
  # Get all cell types in spatial_df
  all_cell_types <- unique(spatial_df[[feature_colname]])
  
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood2D(spatial_df,
                                                                  reference_cell_type,
                                                                  all_cell_types,
                                                                  radius,
                                                                  feature_colname,
                                                                  F,
                                                                  F)
  
  cells_in_neighbourhood_df$total <- rowSums(cells_in_neighbourhood_df[, -1], na.rm = TRUE)
  
  result <- data.frame(reference = reference_cell_type)
  
  # Get total number of cells in spatial_df
  n_cells_in_spatial_df <- length(spatial_df[[feature_colname]])
  
  # Get total number of cells in radius around reference cell type
  n_cells_in_reference_cell_type_radius <- sum(cells_in_neighbourhood_df$total)
  
  for (target_cell_type in target_cell_types) {
    
    # Get total number of target cells in radius around reference cell type
    n_target_cells_in_reference_cell_type_radius <- sum(cells_in_neighbourhood_df[[target_cell_type]])
    
    # Get proportion of target cells in radius around reference cell type
    target_cell_type_proportion_in_reference_cell_type_radius <- n_target_cells_in_reference_cell_type_radius / n_cells_in_reference_cell_type_radius
    
    # Get proportion of target cell type in spatial_df
    n_target_cells_in_spatial_df <- sum(spatial_df[[feature_colname]] == target_cell_type)
    target_cell_type_proportion_in_spatial_df <- n_target_cells_in_spatial_df / n_cells_in_spatial_df
    
    # Get co-occurence value for taget cell type
    target_cell_type_co_occurrence <- target_cell_type_proportion_in_reference_cell_type_radius / target_cell_type_proportion_in_spatial_df
    
    # Add to result data frame
    result[[target_cell_type]] <- target_cell_type_co_occurrence
  }
  
  return(result)
}
calculate_cross_G_gradient2D <- function(spatial_df, 
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
    cross_G_df <- calculate_cross_G2D(spatial_df,
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
calculate_cross_G2D <- function(spatial_df,
                                reference_cell_type,
                                target_cell_type,
                                radius,
                                feature_colname = "Cell.Type") {
  
  ### Calculate the observed cross_G
  # Get the number of target cells in the radius around each reference cell
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood2D(spatial_df,
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
  spatial_df_coords <- spatial_df[ , c("Cell.X.Position", "Cell.Y.Position")]
  
  length <- round(max(spatial_df_coords$Cell.X.Position) - min(spatial_df_coords$Cell.X.Position))
  width  <- round(max(spatial_df_coords$Cell.Y.Position) - min(spatial_df_coords$Cell.Y.Position))
  
  # Get volume of the window the cells are in
  volume <- length * width
  
  # Get the number of target cells
  n_target_cells <- sum(spatial_df[[feature_colname]] == target_cell_type)
  
  # Get target_cell_type intensity (density)
  target_cell_type_intensity <- n_target_cells / volume
  
  # Apply formula
  expected_cross_G <- 1 - exp(-1 * target_cell_type_intensity * pi * radius^2)
  
  result <- data.frame(observed_cross_G = observed_cross_G,
                       expected_cross_G = expected_cross_G)
  
  return(result)
}
calculate_cross_K_gradient2D <- function(spatial_df, 
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
    cross_K_df <- calculate_cross_K2D(spatial_df,
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
calculate_cross_K2D <- function(spatial_df, 
                                reference_cell_type, 
                                target_cell_types, 
                                radius, 
                                feature_colname = "Cell.Type") {
  
  if (is.null(spatial_df[[feature_colname]])) stop(paste("No column called", feature_colname, "found in spatial_df object"))
  
  if (is.null(spatial_df[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spatial_df")
    spatial_df$Cell.ID <- paste("Cell", seq(nrow(spatial_df)), sep = "_")
  }  
  
  
  ## Get expected cross K-function
  expected_cross_K <- pi * radius^2
  
  ## For reference_cell_type, check it is found in the spatial_df object
  if (!(reference_cell_type %in% spatial_df[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spatial_df object"))
    result <- data.frame(observed_cross_K = NA,
                         expected_cross_K = expected_cross_K,
                         cross_K_ratio = NA)
    return(result)
  }
  
  ## Get rough dimensions of the window the points are in
  spatial_df_coords <- spatial_df[ , c("Cell.X.Position", "Cell.Y.Position")]
  
  length <- round(max(spatial_df_coords$Cell.X.Position) - min(spatial_df_coords$Cell.X.Position))
  width  <- round(max(spatial_df_coords$Cell.Y.Position) - min(spatial_df_coords$Cell.Y.Position))
  
  ## Get volume of the window the cells are in
  volume <- length * width
  
  # Number of reference cell types is constant
  n_ref_cells <- sum(spatial_df[[feature_colname]] == reference_cell_type)
  
  # Define result data frame
  result <- data.frame(reference = reference_cell_type, expected = expected_cross_K)
  
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood2D(spatial_df,
                                                                  reference_cell_type,
                                                                  target_cell_types,
                                                                  radius,
                                                                  feature_colname,
                                                                  show_summary = FALSE,
                                                                  plot_image = FALSE)
  
  for (target_cell_type in target_cell_types) {
    
    n_ref_tar_interactions <- sum(cells_in_neighbourhood_df[[target_cell_type]])
    
    n_tar_cells <- sum(spatial_df[[feature_colname]] == target_cell_type)
    
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
calculate_cross_L_gradient2D <- function(spatial_df, 
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
    cross_L_df <- calculate_cross_L2D(spatial_df,
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
calculate_cross_L2D <- function(spatial_df, 
                                reference_cell_type, 
                                target_cell_types, 
                                radius, 
                                feature_colname = "Cell.Type") {
  
  result <- calculate_cross_K2D(spatial_df = spatial_df,
                                reference_cell_type = reference_cell_type,
                                target_cell_types = target_cell_types,
                                radius = radius,
                                feature_colname = feature_colname)
  
  result[ , c("expected", target_cell_types)] <- (result[ , c("expected", target_cell_types)] / (pi)) ^ (1/2)
  
  return(result)
}
calculate_entropy_background2D <- function(spatial_df,
                                           cell_types_of_interest, 
                                           feature_colname = "Cell.Type") {
  
  # NULL case: entropy is undefined
  if (is.null(cell_types_of_interest)) return(NA)
  
  # One cell type case: entropy is 0
  if (is.character(cell_types_of_interest) && length(cell_types_of_interest) == 1) return(0)
  
  cell_proportions_data <- calculate_cell_proportions2D(spatial_df, cell_types_of_interest, feature_colname, FALSE)
  
  # Calculate entropy of the entire image
  entropy <- -1 * sum(cell_proportions_data$proportion * log(cell_proportions_data$proportion, length(cell_proportions_data$proportion)))
  
  return(entropy) 
}
calculate_entropy_gradient2D <- function(spatial_df,
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
    entropy_df <- calculate_entropy2D(spatial_df,
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
    expected_entropy <- calculate_entropy_background2D(spatial_df, target_cell_types, feature_colname)
    fig <- plot_entropy_gradient2D(result, expected_entropy, reference_cell_type, target_cell_types)
    methods::show(fig)
  }
  
  return(result)
}
calculate_entropy_grid_metrics2D <- function(spatial_df, 
                                             n_splits,
                                             cell_types_of_interest,
                                             feature_colname = "Cell.Type",
                                             plot_image = TRUE) {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  ## Check cell_types_of_interest are found in the spatial_df object
  unknown_cell_types <- setdiff(cell_types_of_interest, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  ## If cell types have been chosen, check they are found in the spatial_df object
  unknown_cell_types <- setdiff(cell_types_of_interest, unique(spatial_df[[feature_colname]]))
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Add grid metrics to spatial_df
  grid_metrics <- get_spatial_df_grid_metrics2D(spatial_df, n_splits, feature_colname)
  
  # Get grid_prism_cell_matrix from spatial_df
  grid_prism_cell_matrix <- grid_metrics$grid_prism_cell_matrix
  
  ## Define data frame which contains all results
  n_grid_prisms <- n_splits^2
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
  result <- cbind(result, grid_metrics$grid_prism_coordinates)
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_metrics_continuous2D(result, "entropy")
    methods::show(fig)
  }
  
  return(result)
}
calculate_entropy2D <- function(spatial_df,
                                reference_cell_type,
                                target_cell_types,
                                radius,
                                feature_colname = "Cell.Type") {
  
  # Check target_cell_types
  if (!(is.character(target_cell_types) && length(target_cell_types) >= 2)) {
    stop("`target_cell_types` is not a character vector with at least 2 cell types.")
  }
  
  ## Users should ensure include the reference_cell_type as one of the target_cell_types
  cells_in_neighbourhood_proportion_df <- calculate_cells_in_neighbourhood_proportions2D(spatial_df,
                                                                                         reference_cell_type,
                                                                                         target_cell_types,
                                                                                         radius,
                                                                                         feature_colname)
  
  if (is.null(cells_in_neighbourhood_proportion_df)) return(NULL)
  
  ## Get entropy for target_cell_type
  cells_in_neighbourhood_proportion_df[ , paste(target_cell_types, "_entropy", sep = "")] <- 
    -1 * 
    (cells_in_neighbourhood_proportion_df[ , paste(target_cell_types, "_prop", sep = "")] * log(cells_in_neighbourhood_proportion_df[ , paste(target_cell_types, "_prop", sep = "")], 2) +
       (1 - cells_in_neighbourhood_proportion_df[ , paste(target_cell_types, "_prop", sep = "")]) * log(1 - (cells_in_neighbourhood_proportion_df[ , paste(target_cell_types, "_prop", sep = "")]), 2))
  
  return(cells_in_neighbourhood_proportion_df)
}


calculate_minimum_distances_between_cell_types2D <- function(spatial_df,
                                                             cell_types_of_interest = NULL,
                                                             feature_colname = "Cell.Type",
                                                             show_summary = TRUE,
                                                             plot_image = TRUE) {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
  }
  if (ncol(spatial_df) < 2) {
    stop("There must be at least two cells in spatial_df.")
  }
  if (!(is.null(cell_types_of_interest) || is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(show_summary)) {
    stop("`show_summary` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  if (is.null(spatial_df[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spatial_df")
    spatial_df$Cell.ID <- paste("Cell", seq(nrow(spatial_df)), sep = "_")
  }  
  
  # De-factor feature column in spatial_df object
  spatial_df[[feature_colname]] <- as.character(spatial_df[[feature_colname]])
  
  # Subset spatial_df to only contain the cells of interest
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spatial_df object
    unknown_cell_types <- setdiff(cell_types_of_interest, spatial_df[[feature_colname]])
    if (length(unknown_cell_types) != 0) {
      warning(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
                    paste(unknown_cell_types, collapse = ", ")))
    }
    
    spatial_df <- spatial_df[spatial_df[[feature_colname]] %in% cell_types_of_interest, ]
  }
  # If cell_types_of_interest is NULL, use all cells in spatial_df
  else {
    cell_types_of_interest <- unique(spatial_df[[feature_colname]])
  }
  
  # Create a list containing the cell IDs of each cell type
  cell_type_ids <- list()
  for (cell_type in cell_types_of_interest) {
    cell_type_ids[[cell_type]] <- as.character(spatial_df$Cell.ID[spatial_df[[feature_colname]] == cell_type])
  }
  
  # Get spatial_df coords
  spatial_df_coords <- spatial_df[ , c("Cell.X.Position", "Cell.Y.Position")]
  
  # Get different possible cell type combinations
  # Each row represents a combination
  # If a row is [1 , 2], then we are comparing cell type 1 and cell type 2
  permu <- gtools::permutations(length(cell_types_of_interest), 2, repeats.allowed = TRUE)
  
  result <- data.frame()
  
  for (i in seq(nrow(permu))) {
    cell_type1 <- cell_types_of_interest[permu[i, 1]]
    cell_type2 <- cell_types_of_interest[permu[i, 2]]
    
    # Don't have one of the cells
    if (sum(spatial_df[[feature_colname]] == cell_type1) == 0 || sum(spatial_df[[feature_colname]] == cell_type2) == 0) {
      result <- rbind(result, data.frame(ref_cell_id = NA, ref_cell_type = cell_type1, nearest_cell_id = NA, nearest_cell_type = cell_type2, distance = NA))
      next
    }
    
    # Get x, y, z coords for all cells of cell_type1 and cell_type2
    cell_type1_coords <- spatial_df_coords[spatial_df[[feature_colname]] == cell_type1, ]
    cell_type2_coords <- spatial_df_coords[spatial_df[[feature_colname]] == cell_type2, ]
    
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

calculate_mixing_scores_gradient2D <- function(spatial_df, 
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
    mixing_scores <- calculate_mixing_scores2D(spatial_df,
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
calculate_mixing_scores2D <- function(spatial_df, 
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
      n_ref <- sum(spatial_df[[feature_colname]] == reference_cell_type)
      n_tar <- sum(spatial_df[[feature_colname]] == target_cell_type)
      
      
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
      cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood2D(spatial_df,
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
          methods::show(paste("There are no reference to reference interactions for", target_cell_type, "in the spatial_dfcified radius, cannot calculate mixing score"))
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
calculate_pairwise_distances_between_cell_types2D <- function(spatial_df,
                                                              cell_types_of_interest = NULL,
                                                              feature_colname = "Cell.Type",
                                                              show_summary = TRUE,
                                                              plot_image = TRUE) {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
  }
  if (ncol(spatial_df) < 2) {
    stop("There must be at least two cells in spatial_df.")
  }
  if (!(is.null(cell_types_of_interest) || is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(show_summary)) {
    stop("`show_summary` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  if (is.null(spatial_df[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spatial_df")
    spatial_df$Cell.ID <- paste("Cell", seq(nrow(spatial_df)), sep = "_")
  }
  
  # De-factor feature column in spatial_df object
  spatial_df[[feature_colname]] <- as.character(spatial_df[[feature_colname]])
  
  # Subset spatial_df to only contain the cells of interest
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spatial_df object
    unknown_cell_types <- setdiff(cell_types_of_interest, spatial_df[[feature_colname]])
    if (length(unknown_cell_types) != 0) {
      warning(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
                    paste(unknown_cell_types, collapse = ", ")))
    }
    
    spatial_df <- spatial_df[spatial_df[[feature_colname]] %in% cell_types_of_interest, ]
  }
  # If cell_types_of_interest is NULL, use all cells in spatial_df
  else {
    cell_types_of_interest <- unique(spatial_df[[feature_colname]])
  }
  
  # Create a list containing the cell IDs of each cell type
  cell_type_ids <- list()
  for (cell_type in cell_types_of_interest) {
    cell_type_ids[[cell_type]] <- as.character(spatial_df$Cell.ID[spatial_df[[feature_colname]] == cell_type])
  }
  
  # Calculate cell to cell distances
  distance_matrix <- -1 * apcluster::negDistMat(spatial_df[ , c("Cell.X.Position", "Cell.Y.Position")])
  rownames(distance_matrix) <- spatial_df$Cell.ID
  colnames(distance_matrix) <- spatial_df$Cell.ID
  
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
  
  ## Get splitting number (should be the cube root of n_grid_prisms)
  n_splits <- (n_grid_prisms)^(1/3)
  
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
  ## Adjacent points are within sqrt(2) unit apart. e.g. (0, 0) vs (0, 1)
  else if (weight_method == "queen") {
    weight_matrix <- ifelse(weight_matrix > sqrt(2), 0, 1)  
  }
  ## If a number (x) between 0 and 1 is supplied, set a threshold to be x quantile value of c(weight_matrix)
  ## Grid prisms within this spatial_dfcified threshold have a weight of 1, otherwise, weight of 0
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

get_spatial_df_grid_metrics2D <- function(spatial_df, 
                                          n_splits, 
                                          feature_colname = "Cell.Type") {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  
  spatial_df_coords <- spatial_df[ , c("Cell.X.Position", "Cell.Y.Position")]
  
  ## Get dimensions of the window
  min_x <- min(spatial_df_coords[ , "Cell.X.Position"])
  min_y <- min(spatial_df_coords[ , "Cell.Y.Position"])
  
  max_x <- max(spatial_df_coords[ , "Cell.X.Position"])
  max_y <- max(spatial_df_coords[ , "Cell.Y.Position"])
  
  length <- round(max_x - min_x)
  width  <- round(max_y - min_y)
  
  ## Get distance of row, col and lay
  d_row <- length / n_splits
  d_col <- width / n_splits
  
  # Shift spatial_df_coords so they begin at the origin
  spatial_df_coords[, "Cell.X.Position"] <- spatial_df_coords[, "Cell.X.Position"] - min_x
  spatial_df_coords[, "Cell.Y.Position"] <- spatial_df_coords[, "Cell.Y.Position"] - min_y
  
  ## Figure out which 'grid prism number' each cell is inside
  spatial_df$grid_prism_num <- floor(spatial_df_coords[ , "Cell.X.Position"] / d_row) +
    floor(spatial_df_coords[ , "Cell.Y.Position"] / d_col) * n_splits
  
  ## Determine the cell types found in each grid prism
  n_grid_prisms <- n_splits^2
  grid_prism_cell_matrix <- as.data.frame.matrix(table(spatial_df[[feature_colname]], factor(spatial_df$grid_prism_num, levels = seq(n_grid_prisms))))
  grid_prism_cell_matrix <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                       t(grid_prism_cell_matrix), check.names = FALSE)
  
  ## Determine centre coordinates of each grid prism
  grid_prism_coordinates <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                       x_coord = ((seq(n_grid_prisms) - 1) %% n_splits + 0.5) * d_row + round(min_x),
                                       y_coord = (floor(((seq(n_grid_prisms) - 1) %% (n_splits)^2) / n_splits) + 0.5) * d_col + round(min_y))
  
  grid_metrics <- list("grid_prism_cell_matrix" = grid_prism_cell_matrix,
                       "grid_prism_coordinates" = grid_prism_coordinates)
  
  return(grid_metrics)
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


# Analysis pipeline function -----
# ******** READ THIS *********
# This function assumes the following:
# data3D is a data frame with column names: "Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position", ""Cell.Type"
# "Cell.Z.Position" contains discrete values.
# cell_types consists of only cells in the "Cell.Type" column, and at least two cells.

analyse_3D_data <- function(
    data3D,
    cell_types,
    radii = seq(20, 100, 10),
    n_splits = 10,
    thresholds = seq(0.01, 1, 0.01)
) {
  
  slice_z_coords <- unique(data3D$Cell.Z.Position)
  n_slices <- length(unique(data3D$Cell.Z.Position))
  n_cell_type_combinations <- length(cell_types)^2
  radii_colnames <- paste("r", radii, sep = "")
  thresholds_colnames <- paste("t", thresholds, sep = "")
  
  create_empty_metric_df_list <- function(
    cell_types,
    n_slices,
    radii_colnames,
    thresholds_colnames
  ) {
    n_cell_type_combinations <- length(cell_types)^2
    
    # Define AMD data frames as well as constants
    AMD_df_colnames <- c("slice", "reference", "target", "AMD")
    AMD_df <- data.frame(matrix(nrow = (n_slices + 1) * n_cell_type_combinations, ncol = length(AMD_df_colnames)))
    colnames(AMD_df) <- AMD_df_colnames
    
    # Define MS, NMS, ACIN, ACINP, CKR, CLR, CGR, COO, AE data frames as well as constants
    radii_colnames <- paste("r", radii, sep = "")
    
    MS_df_colnames <- c("slice", "reference", "target", radii_colnames)
    MS_df <- data.frame(matrix(nrow = (n_slices + 1) * n_cell_type_combinations, ncol = length(MS_df_colnames)))
    colnames(MS_df) <- MS_df_colnames
    
    NMS_df <- ACIN_df <- AE_df <- ACINP_df <- CKR_df <- CLR_df <- COO_df <- CGR_df <- MS_df
    
    # Define SAC and prevalence data frames as well as constants
    thresholds_colnames <- paste("t", thresholds, sep = "")
    
    PBSAC_df_colnames <- c("slice", "reference", "target", "PBSAC")
    PBSAC_df <- data.frame(matrix(nrow = (n_slices + 1) * n_cell_type_combinations, ncol = length(PBSAC_df_colnames)))
    colnames(PBSAC_df) <- PBSAC_df_colnames
    
    PBP_df_colnames <- c("slice", "reference", "target", thresholds_colnames)
    PBP_df <- data.frame(matrix(nrow = (n_slices + 1) * n_cell_type_combinations, ncol = length(PBP_df_colnames)))
    colnames(PBP_df) <- PBP_df_colnames
    
    EBSAC_df_colnames <- c("slice", "cell_types", "EBSAC")
    EBSAC_df <- data.frame(matrix(nrow = (n_slices + 1) * n_cell_type_combinations, ncol = length(EBSAC_df_colnames)))
    colnames(EBSAC_df) <- EBSAC_df_colnames
    
    EBP_df_colnames <- c("slice", "cell_types", thresholds_colnames)
    EBP_df <- data.frame(matrix(nrow = (n_slices + 1) * n_cell_type_combinations, ncol = length(EBP_df_colnames)))
    colnames(EBP_df) <- EBP_df_colnames
    
    
    # Add all to list:
    metric_df_list <- list(AMD = AMD_df,
                           MS = MS_df,
                           NMS = NMS_df,
                           ACINP = ACINP_df,
                           AE = AE_df,
                           ACIN = ACIN_df,
                           CKR = CKR_df,
                           CLR = CLR_df,
                           CGR = CGR_df,
                           COO = COO_df,
                           PBSAC = PBSAC_df,
                           PBP = PBP_df,
                           EBSAC = EBSAC_df,
                           EBP = EBP_df)
    
    return(metric_df_list)
  }
  
  metric_df_list <- create_empty_metric_df_list(
    cell_types,
    n_slices,
    radii_colnames,
    thresholds_colnames
  )
  
  for (i in seq(n_slices + 1, 1)) {
    print(i)
    
    # i represents the current slice index
    if (i == n_slices + 1) {
      df <- data3D
    }
    # if i is the last index, analyse in 3D instead
    else {
      df <- data3D[data3D$Cell.Z.Position == slice_z_coords[i], ]
    }
    
    ### 3D analysis -----------------------------
    if (i == n_slices + 1) {
      
      index <- n_cell_type_combinations * (i - 1) + 1 
      
      minimum_distance_data <- calculate_minimum_distances_between_cell_types3D(df,
                                                                                cell_types,
                                                                                show_summary = F,
                                                                                plot_image = F)
      
      minimum_distance_data_summary <- summarise_distances_between_cell_types3D(minimum_distance_data)
      
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "slice"] <- i
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "reference"] <- minimum_distance_data_summary$reference
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "target"] <- minimum_distance_data_summary$target
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "AMD"] <- minimum_distance_data_summary$mean
      
      # Need a new index for gradient-based data which increments after each target cell type
      gradient_index <- n_cell_type_combinations * (i - 1) + 1 
      
      for (reference_cell_type in cell_types) {
        gradient_data <- calculate_all_gradient_cc_metrics3D(df,
                                                             reference_cell_type,
                                                             cell_types,
                                                             radii,
                                                             plot_image = F)
        
        for (target_cell_type in cell_types) {
          print(paste(reference_cell_type, target_cell_type, sep = "/"))
          metric_df_list[["ACIN"]][gradient_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["ACINP"]][gradient_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["CKR"]][gradient_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["CLR"]][gradient_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["COO"]][gradient_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["CGR"]][gradient_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["MS"]][gradient_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["NMS"]][gradient_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["AE"]][gradient_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, 
                                                                                         paste(reference_cell_type, target_cell_type, sep = ","))
          
          if (reference_cell_type == target_cell_type) {
            metric_df_list[["ACINP"]][gradient_index, radii_colnames] <- Inf
            metric_df_list[["MS"]][gradient_index, radii_colnames] <- Inf
            metric_df_list[["NMS"]][gradient_index, radii_colnames] <- Inf
            metric_df_list[["AE"]][gradient_index, radii_colnames] <- Inf
          }
          
          if (is.null(gradient_data)) {
            metric_df_list[["ACIN"]][gradient_index, radii_colnames] <- NA
            metric_df_list[["CKR"]][gradient_index, radii_colnames] <- NA
            metric_df_list[["CLR"]][gradient_index, radii_colnames] <- NA
            metric_df_list[["COO"]][gradient_index, radii_colnames] <- NA
            metric_df_list[["CGR"]][gradient_index, radii_colnames] <- NA
          }
          else {
            metric_df_list[["ACIN"]][gradient_index, radii_colnames] <- gradient_data[["cells_in_neighbourhood"]][[target_cell_type]]
            metric_df_list[["CKR"]][gradient_index, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]] / gradient_data[["cross_K"]][["expected"]]
            metric_df_list[["CLR"]][gradient_index, radii_colnames] <- gradient_data[["cross_L"]][[target_cell_type]] / gradient_data[["cross_L"]][["expected"]]
            metric_df_list[["COO"]][gradient_index, radii_colnames] <- gradient_data[["co_occurrence"]][[target_cell_type]]
            metric_df_list[["CGR"]][gradient_index, radii_colnames] <- gradient_data[["cross_G"]][[target_cell_type]][["observed_cross_G"]] / gradient_data[["cross_G"]][[target_cell_type]][["expected_cross_G"]]
            
            if (reference_cell_type != target_cell_type) {
              metric_df_list[["ACINP"]][gradient_index, radii_colnames] <- gradient_data[["cells_in_neighbourhood_proportion"]][[target_cell_type]]
              metric_df_list[["MS"]][gradient_index, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$mixing_score
              metric_df_list[["NMS"]][gradient_index, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$normalised_mixing_score
              metric_df_list[["AE"]][gradient_index, radii_colnames] <- gradient_data[["entropy"]][[target_cell_type]]
            }
          }
          
          # Spatial heterogeneity metrics
          metric_df_list[["PBSAC"]][gradient_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["PBP"]][gradient_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          
          metric_df_list[["EBSAC"]][gradient_index, c("slice", "cell_types")] <- c(i, paste(reference_cell_type, target_cell_type, sep = ","))
          metric_df_list[["EBP"]][gradient_index, c("slice", "cell_types")] <- c(i, paste(reference_cell_type, target_cell_type, sep = ","))
          
          if (reference_cell_type != target_cell_type) {
            proportion_grid_metrics <- calculate_cell_proportion_grid_metrics3D(df, 
                                                                                n_splits,
                                                                                reference_cell_type, 
                                                                                target_cell_type,
                                                                                plot_image = F)
            
            if (is.null(proportion_grid_metrics)) {
              metric_df_list[["PBSAC"]][gradient_index, "PBSAC"] <- NA
              metric_df_list[["PBP"]][gradient_index, thresholds_colnames] <- NA
            }
            else {
              PBSAC <- calculate_spatial_autocorrelation3D(proportion_grid_metrics, 
                                                           "proportion",
                                                           weight_method = 0.1)
              
              PBP_df <- calculate_prevalence_gradient3D(proportion_grid_metrics,
                                                        "proportion",
                                                        show_AUC = F,
                                                        plot_image = F)
              
              
              metric_df_list[["PBSAC"]][gradient_index, "PBSAC"] <- PBSAC
              metric_df_list[["PBP"]][gradient_index, thresholds_colnames] <- PBP_df$prevalence
            } 
            
            entropy_grid_metrics <- calculate_entropy_grid_metrics3D(df, 
                                                                     n_splits,
                                                                     c(reference_cell_type, target_cell_type), 
                                                                     plot_image = F)
            
            if (is.null(entropy_grid_metrics)) {
              metric_df_list[["EBSAC"]][gradient_index, "EBSAC"] <- NA
              metric_df_list[["EBP"]][gradient_index, thresholds_colnames] <- NA
            }
            else {
              EBSAC <- calculate_spatial_autocorrelation3D(entropy_grid_metrics, 
                                                           "entropy",
                                                           weight_method = 0.1)
              
              EBP_df <- calculate_prevalence_gradient3D(entropy_grid_metrics,
                                                        "entropy",
                                                        show_AUC = F,
                                                        plot_image = F)
              
              metric_df_list[["EBSAC"]][gradient_index, "EBSAC"] <- EBSAC
              metric_df_list[["EBP"]][gradient_index, thresholds_colnames] <- EBP_df$prevalence
            }    
          }
          else {
            metric_df_list[["PBSAC"]][gradient_index, "PBSAC"] <- Inf
            metric_df_list[["PBP"]][gradient_index, thresholds_colnames] <- Inf
            metric_df_list[["EBSAC"]][gradient_index, "EBSAC"] <- Inf
            metric_df_list[["EBP"]][gradient_index, thresholds_colnames] <- Inf
          }
          
          gradient_index <- gradient_index + 1
        }
      }
    }
    ### 2D analysis -----------------------------
    else {
      
      index <- n_cell_type_combinations * (i - 1) + 1 
      
      minimum_distance_data <- calculate_minimum_distances_between_cell_types2D(df,
                                                                                cell_types,
                                                                                show_summary = F,
                                                                                plot_image = F)
      
      minimum_distance_data_summary <- summarise_distances_between_cell_types2D(minimum_distance_data)
      
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "slice"] <- i
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "reference"] <- minimum_distance_data_summary$reference
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "target"] <- minimum_distance_data_summary$target
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "AMD"] <- minimum_distance_data_summary$mean
      
      # Need a new index for gradient-based data which increments after each target cell type
      gradient_index <- n_cell_type_combinations * (i - 1) + 1 
      
      for (reference_cell_type in cell_types) {
        gradient_data <- calculate_all_gradient_cc_metrics2D(df,
                                                             reference_cell_type,
                                                             cell_types,
                                                             radii,
                                                             plot_image = F)
        
        for (target_cell_type in cell_types) {
          print(paste(reference_cell_type, target_cell_type, sep = "/"))
          metric_df_list[["ACIN"]][gradient_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["ACINP"]][gradient_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["CKR"]][gradient_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["CLR"]][gradient_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["COO"]][gradient_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["CGR"]][gradient_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["MS"]][gradient_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["NMS"]][gradient_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["AE"]][gradient_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, 
                                                                                         paste(reference_cell_type, target_cell_type, sep = ","))
          
          if (reference_cell_type == target_cell_type) {
            metric_df_list[["ACINP"]][gradient_index, radii_colnames] <- Inf
            metric_df_list[["MS"]][gradient_index, radii_colnames] <- Inf
            metric_df_list[["NMS"]][gradient_index, radii_colnames] <- Inf
            metric_df_list[["AE"]][gradient_index, radii_colnames] <- Inf
          }
          
          if (is.null(gradient_data)) {
            metric_df_list[["ACIN"]][gradient_index, radii_colnames] <- NA
            metric_df_list[["CKR"]][gradient_index, radii_colnames] <- NA
            metric_df_list[["CLR"]][gradient_index, radii_colnames] <- NA
            metric_df_list[["COO"]][gradient_index, radii_colnames] <- NA
            metric_df_list[["CGR"]][gradient_index, radii_colnames] <- NA
          }
          else {
            metric_df_list[["ACIN"]][gradient_index, radii_colnames] <- gradient_data[["cells_in_neighbourhood"]][[target_cell_type]]
            metric_df_list[["CKR"]][gradient_index, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]] / gradient_data[["cross_K"]][["expected"]]
            metric_df_list[["CLR"]][gradient_index, radii_colnames] <- gradient_data[["cross_L"]][[target_cell_type]] / gradient_data[["cross_L"]][["expected"]]
            metric_df_list[["COO"]][gradient_index, radii_colnames] <- gradient_data[["co_occurrence"]][[target_cell_type]]
            metric_df_list[["CGR"]][gradient_index, radii_colnames] <- gradient_data[["cross_G"]][[target_cell_type]][["observed_cross_G"]] / gradient_data[["cross_G"]][[target_cell_type]][["expected_cross_G"]]
            
            if (reference_cell_type != target_cell_type) {
              metric_df_list[["ACINP"]][gradient_index, radii_colnames] <- gradient_data[["cells_in_neighbourhood_proportion"]][[target_cell_type]]
              metric_df_list[["MS"]][gradient_index, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$mixing_score
              metric_df_list[["NMS"]][gradient_index, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$normalised_mixing_score
              metric_df_list[["AE"]][gradient_index, radii_colnames] <- gradient_data[["entropy"]][[target_cell_type]]
            }
          }
          
          # Spatial heterogeneity metrics
          metric_df_list[["PBSAC"]][gradient_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["PBP"]][gradient_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          
          metric_df_list[["EBSAC"]][gradient_index, c("slice", "cell_types")] <- c(i, paste(reference_cell_type, target_cell_type, sep = ","))
          metric_df_list[["EBP"]][gradient_index, c("slice", "cell_types")] <- c(i, paste(reference_cell_type, target_cell_type, sep = ","))
          
          if (reference_cell_type != target_cell_type) {
            proportion_grid_metrics <- calculate_cell_proportion_grid_metrics2D(df, 
                                                                                n_splits,
                                                                                reference_cell_type, 
                                                                                target_cell_type,
                                                                                plot_image = F)
            
            if (is.null(proportion_grid_metrics)) {
              metric_df_list[["PBSAC"]][gradient_index, "PBSAC"] <- NA
              metric_df_list[["PBP"]][gradient_index, thresholds_colnames] <- NA
            }
            else {
              PBSAC <- calculate_spatial_autocorrelation2D(proportion_grid_metrics, 
                                                           "proportion",
                                                           weight_method = 0.1)
              
              PBP_df <- calculate_prevalence_gradient2D(proportion_grid_metrics,
                                                        "proportion",
                                                        show_AUC = F,
                                                        plot_image = F)
              
              
              metric_df_list[["PBSAC"]][gradient_index, "PBSAC"] <- PBSAC
              metric_df_list[["PBP"]][gradient_index, thresholds_colnames] <- PBP_df$prevalence
            } 
            
            entropy_grid_metrics <- calculate_entropy_grid_metrics2D(df, 
                                                                     n_splits,
                                                                     c(reference_cell_type, target_cell_type), 
                                                                     plot_image = F)
            
            if (is.null(entropy_grid_metrics)) {
              metric_df_list[["EBSAC"]][gradient_index, "EBSAC"] <- NA
              metric_df_list[["EBP"]][gradient_index, thresholds_colnames] <- NA
            }
            else {
              EBSAC <- calculate_spatial_autocorrelation2D(entropy_grid_metrics, 
                                                           "entropy",
                                                           weight_method = 0.1)
              
              EBP_df <- calculate_prevalence_gradient2D(entropy_grid_metrics,
                                                        "entropy",
                                                        show_AUC = F,
                                                        plot_image = F)
              
              metric_df_list[["EBSAC"]][gradient_index, "EBSAC"] <- EBSAC
              metric_df_list[["EBP"]][gradient_index, thresholds_colnames] <- EBP_df$prevalence
            }    
          }
          else {
            metric_df_list[["PBSAC"]][gradient_index, "PBSAC"] <- Inf
            metric_df_list[["PBP"]][gradient_index, thresholds_colnames] <- Inf
            metric_df_list[["EBSAC"]][gradient_index, "EBSAC"] <- Inf
            metric_df_list[["EBP"]][gradient_index, thresholds_colnames] <- Inf
          }
          
          gradient_index <- gradient_index + 1
        }
      }
    }
  }
  return(metric_df_list)
}



# Analysis and upload ----
metric_df_list <- analyse_3D_data(data3D, cell_types)
setwd(save_directory)
saveRDS(metric_df_list, file_name)