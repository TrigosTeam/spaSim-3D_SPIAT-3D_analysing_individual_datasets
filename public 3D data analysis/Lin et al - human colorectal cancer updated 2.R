### Read data -----

# In raw data, each cell is defined by a number corresponding to its cell type
# Use the cell type dictionary to determine the cell type of a cell from its number
setwd("~/R/Lin et al - human colorectal cancer/other_data")
cell_type_dictionary <- read.csv("cell_type_dictionary.csv")
specific_cell_types <- unique(cell_type_dictionary$Type_Name) # E.g. Tumour/Epi, Ki67+ Tumour/Epi...
generic_cell_types <- unique(cell_type_dictionary$Category) # E.g. Tumour, Stroma, Immune

setwd("~/R/Lin et al - human colorectal cancer/CRC1_data_updated/raw_slice_data")

slice_nums <- c("002", "007", "014", "020", "025", "029",
                "034", "039", "044", "049", "050", "051",
                "052", "054", "059", "064", "069", "074",
                "078", "084", "086", "091", "097", "102", "106")
file_names <- paste("CRC01-", slice_nums, ".csv", sep = "")

# Get z coord of each slice in micrometers using the slice number
get_slice_z_coord <- function(slice_num) {
  
  slice_num <- as.integer(slice_num)
  
  if (slice_num <= 85) {
    return(slice_num * 5)
  }
  
  else {
    return(85 * 5 + (slice_num - 85) * 4)
  }
}

# Get 2D data
# Store each slice as a data.frame, all in a big list
list_dfs2D <- list()

# Loop through each slice
for (i in seq(length(slice_nums))) {
  
  file_name <- file_names[i]
  
  curr_slice_df <- read.csv(file_name, sep = ",")
  curr_slice_df <- curr_slice_df[ , c("NewType", "Xtt", "Ytt")]
  curr_slice_df$Ztt <- get_slice_z_coord(slice_nums[i])
  
  curr_slice_df$Cell.Type.Specific <- specific_cell_types[curr_slice_df$NewType]
  curr_slice_df$Cell.Type.Generic <- generic_cell_types[ifelse(curr_slice_df$NewType <= 3, 1, 
                                                               ifelse(curr_slice_df$NewType <= 6, 2, 3))]
  
  curr_slice_df$NewType <- NULL
  
  list_dfs2D[[i]] <- curr_slice_df
}

# Get 3D data
# Store as one data frame
df3D <- do.call(rbind, list_dfs2D)

setwd("~/R/Lin et al - human colorectal cancer/CRC1_data_final")
saveRDS(list_dfs2D, "list_dfs2D.RDS")
saveRDS(df3D, "df3D.RDS")


### Set up dfs to contain analysis results -----
setwd("~/R/Lin et al - human colorectal cancer/CRC1_data_final")
df3D <- readRDS("df3D.RDS")
colnames(df3D)[1:3] <- c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")
n_slices <- length(unique(df3D$Cell.Z.Position))

# Define AMD data frames as well as constants
cell_types <- c("Tumour", "Immune")
AMD_pairs <- c("Tumour/Tumour", "Tumour/Immune", "Immune/Tumour", "Immune/Immune")

AMD_df_colnames <- c("slice", "reference", "target", "AMD")
AMD_df <- data.frame(matrix(nrow = length(AMD_pairs) * (n_slices + 1), ncol = length(AMD_df_colnames)))
colnames(AMD_df) <- AMD_df_colnames

# Define MS, NMS, ACIN, ACINP, CKR, AE data frames as well as constants
radii <- seq(20, 100, 10)
radii_colnames <- paste("r", radii, sep = "")

MS_df_colnames <- c("slice", "reference", "target", radii_colnames)
MS_df <- data.frame(matrix(nrow = length(cell_types) * (n_slices + 1), ncol = length(MS_df_colnames)))
colnames(MS_df) <- MS_df_colnames

NMS_df_colnames <- c("slice", "reference", "target", radii_colnames)
NMS_df <- data.frame(matrix(nrow = length(cell_types) * (n_slices + 1), ncol = length(MS_df_colnames)))
colnames(NMS_df) <- NMS_df_colnames

# Target is always Tumour and Immune together
# Only choose prop(Immune) as prop(Tumour) = 1 - prop(Immune) always
ACINP_df_colnames <- c("slice", "reference", "target", radii_colnames)
ACINP_df <- data.frame(matrix(nrow = length(cell_types) * (n_slices + 1), ncol = length(ACINP_df_colnames)))
colnames(ACINP_df) <- ACINP_df_colnames

# Target is always Tumour and Immune together
AE_df_colnames <- c("slice", "reference", "target", radii_colnames)
AE_df <- data.frame(matrix(nrow = length(cell_types) * (n_slices + 1), ncol = length(AE_df_colnames)))
colnames(AE_df) <- AE_df_colnames

## ACIN and CKR are twice as large
# (ref Tumour and tar Tumour or Immune) OR (ref Immune and tar Immune or Tumour)
ACIN_df_colnames <- c("slice", "reference", "target", radii_colnames)
ACIN_df <- data.frame(matrix(nrow = length(cell_types)^2 * (n_slices + 1), ncol = length(ACIN_df_colnames)))
colnames(ACIN_df) <- ACIN_df_colnames

# (ref Tumour and tar Tumour or Immune) OR (ref Immune and tar Immune or Tumour)
CKR_df_colnames <- c("slice", "reference", "target", radii_colnames)
CKR_df <- data.frame(matrix(nrow = length(cell_types)^2 * (n_slices + 1), ncol = length(CKR_df_colnames)))
colnames(CKR_df) <- CKR_df_colnames

# Define SAC and prevalence data frames as well as constants
n_splits <- 10
thresholds <- seq(0.01, 1, 0.01)
thresholds_colnames <- paste("t", thresholds, sep = "")

prop_cell_types <- data.frame(ref = c("Tumour"), tar = c("Immune"))

prop_SAC_df_colnames <- c("slice", "reference", "target", "prop_SAC")
prop_SAC_df <- data.frame(matrix(nrow = nrow(prop_cell_types) * (n_slices + 1), ncol = length(prop_SAC_df_colnames)))
colnames(prop_SAC_df) <- prop_SAC_df_colnames

prop_prevalence_df_colnames <- c("slice", "reference", "target", thresholds_colnames)
prop_prevalence_df <- data.frame(matrix(nrow = nrow(prop_cell_types) * (n_slices + 1), ncol = length(prop_prevalence_df_colnames)))
colnames(prop_prevalence_df) <- prop_prevalence_df_colnames


entropy_cell_types <- data.frame(cell_types = c("Tumour,Immune"))

entropy_SAC_df_colnames <- c("slice", "cell_types", "entropy_SAC")
entropy_SAC_df <- data.frame(matrix(nrow = nrow(entropy_cell_types) * (n_slices + 1), ncol = length(entropy_SAC_df_colnames)))
colnames(entropy_SAC_df) <- entropy_SAC_df_colnames

entropy_prevalence_df_colnames <- c("slice", "cell_types", thresholds_colnames)
entropy_prevalence_df <- data.frame(matrix(nrow = nrow(entropy_cell_types) * (n_slices + 1), ncol = length(entropy_prevalence_df_colnames)))
colnames(entropy_prevalence_df) <- entropy_prevalence_df_colnames


# Add all to list:
metric_df_list <- list(AMD = AMD_df,
                       MS = MS_df,
                       NMS = NMS_df,
                       ACINP = ACINP_df,
                       AE = AE_df,
                       ACIN = ACIN_df,
                       CKR = CKR_df,
                       prop_SAC = prop_SAC_df,
                       prop_prevalence = prop_prevalence_df,
                       entropy_SAC = entropy_SAC_df,
                       entropy_prevalence = entropy_prevalence_df)


### SPIAT-3D / SPIAT functions (using df as input) ------

calculate_minimum_distances_between_cell_types <- function(df,
                                                           cell_types_of_interest = NULL,
                                                           feature_colname = "Cell.Type",
                                                           dimension) {
  
  if (is.null(df[[feature_colname]])) stop(paste("No column called", feature_colname, "found in df object"))
  
  if (is.null(df[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your df")
    df$Cell.ID <- paste("Cell", seq(nrow(df)), sep = "_")
  }  
  
  # If there are less than two cells, give error
  if (nrow(df) < 2) stop("There must be at least two cells in df")
  
  # Subset df to only contain the cells of interest
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the df object
    unknown_cell_types <- setdiff(cell_types_of_interest, df[[feature_colname]])
    if (length(unknown_cell_types) != 0) {
      warning(paste("The following cell types in cell_types_of_interest are not found in the df object:\n   ",
                    paste(unknown_cell_types, collapse = ", ")))
    }
    
    df <- df[df[[feature_colname]] %in% cell_types_of_interest , ]
  }
  # If cell_types_of_interest is NULL, use all cells in df
  else {
    cell_types_of_interest <- unique(df[[feature_colname]])
  }
  
  # Create a list containing the cell IDs of each cell type
  cell_type_ids <- list()
  for (cell_type in cell_types_of_interest) {
    cell_type_ids[[cell_type]] <- as.character(df$Cell.ID[df[[feature_colname]] == cell_type])
  }
  
  # Get df coords
  if (dimension == "3D") {
    df_coords <- df[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]  
  }
  else if (dimension == "2D") {
    df_coords <- df[ , c("Cell.X.Position", "Cell.Y.Position")]
  }
  else {
    stop("Invalid dimension. Choose either 3D or 2D")
  }
  
  
  # Get different possible cell type combinations
  # Each row represents a combination
  # If a row is [1 , 2], then we are comparing cell type 1 and cell type 2
  permu <- gtools::permutations(length(cell_types_of_interest), 2, repeats.allowed = TRUE)
  result <- data.frame()
  
  for (i in seq(nrow(permu))) {
    cell_type1 <- cell_types_of_interest[permu[i, 1]]
    cell_type2 <- cell_types_of_interest[permu[i, 2]]
    
    # Don't have one of the cells
    if (sum(df[[feature_colname]] == cell_type1) == 0 || sum(df[[feature_colname]] == cell_type2) == 0) {
      result <- rbind(result, data.frame(ref_cell_id = NA, ref_cell_type = cell_type1, nearest_cell_id = NA, nearest_cell_type = cell_type2, distance = NA))
      next
    }
    
    # Get x, y, z coords for all cells of cell_type1 and cell_type2
    cell_type1_coords <- df_coords[df[[feature_colname]] == cell_type1, ]
    cell_type2_coords <- df_coords[df[[feature_colname]] == cell_type2, ]
    
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
    
    curr_pair_df <- data.frame(
      ref_cell_id = cell_type_ids[[cell_type1]],
      ref_cell_type = cell_type1,
      nearest_cell_id = cell_type_ids[[cell_type2]][c(nearest_neighbours$nn.idx)],
      nearest_cell_type = cell_type2,
      distance = nearest_neighbours$nn.dists
    )
    result <- rbind(result, curr_pair_df)
  }
  
  result$pair <- paste(result$ref_cell_type, result$nearest_cell_type,sep = "/")
  
  return(result)
}




summarise_distances_between_cell_types <- function(distances_df) {
  
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



calculate_all_gradient_cc_metrics <- function(df, 
                                              reference_cell_type, 
                                              target_cell_types, 
                                              radii, 
                                              feature_colname = "Cell.Type", 
                                              dimension) {
  
  
  ## Define result
  result <- list("mixing_score" = list(),
                 "cells_in_neighbourhood" = data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types))),
                 "cells_in_neighbourhood_proportion" = data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types))),
                 "entropy" = data.frame(matrix(nrow = length(radii), ncol = 1)),
                 "cross_K" = list())
  colnames(result[["cells_in_neighbourhood"]]) <- target_cell_types
  colnames(result[["cells_in_neighbourhood_proportion"]]) <- target_cell_types
  colnames(result[["entropy"]]) <- "entropy"
  
  # Define other constants
  mixing_score_df_colnames <- c("ref_cell_type", 
                                "tar_cell_type", 
                                "n_ref_cells",
                                "n_tar_cells", 
                                "n_ref_tar_interactions",
                                "n_ref_ref_interactions", 
                                "mixing_score", 
                                "normalised_mixing_score")
  cross_K_df_colnames <- c("ref_cell_type",
                           "tar_cell_type",
                           "observed_cross_K",
                           "expected_cross_K",
                           "cross_K_ratio")
  
  # Define indiviudal data frames for mixing_score and cross_K
  for (target_cell_type in target_cell_types) {
    if (reference_cell_type != target_cell_type) {
      result[["mixing_score"]][[target_cell_type]] <- data.frame(matrix(nrow = length(radii), ncol = length(mixing_score_df_colnames)))
      colnames(result[["mixing_score"]][[target_cell_type]]) <- mixing_score_df_colnames
    }
    result[["cross_K"]][[target_cell_type]] <- data.frame(matrix(nrow = length(radii), ncol = length(cross_K_df_colnames)))
    colnames(result[["cross_K"]][[target_cell_type]]) <- cross_K_df_colnames
  }
  
  # Get gradient results for each metric
  for (i in seq(length(radii))) {
    df_gradient_results <- calculate_all_single_radius_cc_metrics(df,
                                                                  reference_cell_type,
                                                                  target_cell_types,
                                                                  radii[i],
                                                                  feature_colname,
                                                                  dimension)
    
    if (is.null(df_gradient_results)) return(NULL)
    
    df_gradient_results[["cells_in_neighbourhood"]]$ref_cell_id <- NULL
    
    result[["cells_in_neighbourhood"]][i, ] <- apply(df_gradient_results[["cells_in_neighbourhood"]], 2, mean)
    result[["cells_in_neighbourhood_proportion"]][i, ] <- apply(df_gradient_results[["cells_in_neighbourhood_proportion"]][ , paste(target_cell_types, "_prop", sep = "")], 2, mean, na.rm = T)
    result[["entropy"]][i, "entropy"] <- mean(df_gradient_results[["entropy"]]$entropy, na.rm = T)
    
    for (target_cell_type in names(df_gradient_results[["mixing_score"]])) {
      result[["mixing_score"]][[target_cell_type]][i, ] <- df_gradient_results[["mixing_score"]][[target_cell_type]]
    }
    
    for (target_cell_type in names(df_gradient_results[["cross_K"]])) {
      result[["cross_K"]][[target_cell_type]][i, ] <- df_gradient_results[["cross_K"]][[target_cell_type]]
    }
  }
  
  # Add radius column to each data frame
  result[["cells_in_neighbourhood"]]$radius <- radii
  result[["cells_in_neighbourhood_proportion"]]$radius <- radii
  result[["entropy"]]$radius <- radii
  for (target_cell_type in names(df_gradient_results[["mixing_score"]])) {
    result[["mixing_score"]][[target_cell_type]]$radius <- radii
  }
  
  for (target_cell_type in names(df_gradient_results[["cross_K"]])) {
    result[["cross_K"]][[target_cell_type]]$radius <- radii
  }
  
  return(result)
}


### Calculate all single radius cell-colocalisation metrics
# If a function only requires one target cell type, iterate through each cell type in target_cell_types, else use all target_cell_types

calculate_all_single_radius_cc_metrics <- function(df, 
                                                   reference_cell_type, 
                                                   target_cell_types, 
                                                   radius, 
                                                   feature_colname = "Cell.Type",
                                                   dimension) {
  
  if (is.null(df[[feature_colname]])) stop(paste("No column called", feature_colname, "found in df object"))
  
  ## For reference_cell_type, check it is found in the df object
  if (!(reference_cell_type %in% df[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the df object"))
    return(NULL)
  }
  
  ## For target_cell_types, check they are found in the df object
  unknown_cell_types <- setdiff(target_cell_types, df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
  }
  
  # Check if radius is numeric
  if (!is.numeric(radius)) stop(paste(radius, " is not of type 'numeric'"))
  
  
  # Define result
  result <- list("cells_in_neighbourhood" = list(),
                 "cells_in_neighbourhood_proportion" = list(),
                 "entropy" = list(),
                 "mixing_score" = list(),
                 "cross_K" = list())
  
  # Define other constants
  mixing_score_df_colnames <- c("ref_cell_type", 
                                "tar_cell_type", 
                                "n_ref_cells",
                                "n_tar_cells", 
                                "n_ref_tar_interactions",
                                "n_ref_ref_interactions", 
                                "mixing_score", 
                                "normalised_mixing_score")
  cross_K_df_colnames <- c("ref_cell_type",
                           "tar_cell_type",
                           "observed_cross_K",
                           "expected_cross_K",
                           "cross_K_ratio")
  
  # Get rough dimensions of window for cross_K
  if (dimension == "3D") {
    df_coords <- df[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]  
    length <- round(max(df_coords$Cell.X.Position) - min(df_coords$Cell.X.Position))
    width  <- round(max(df_coords$Cell.Y.Position) - min(df_coords$Cell.Y.Position))
    height <- round(max(df_coords$Cell.Z.Position) - min(df_coords$Cell.Z.Position))
    
    ## Get volume of the window the cells are in
    volume <- length * width * height
  }
  else if (dimension == "2D") {
    df_coords <- df[ , c("Cell.X.Position", "Cell.Y.Position")]
    length <- round(max(df_coords$Cell.X.Position) - min(df_coords$Cell.X.Position))
    width  <- round(max(df_coords$Cell.Y.Position) - min(df_coords$Cell.Y.Position))
    
    
    ## Get area of the window the cells are in
    area <- length * width
  }
  else {
    stop("Invalid dimension. Choose either 3D or 2D")
  }
  
  
  
  
  # All single radius cc metrics stem from calculate_entropy3D function
  entropy_df <- calculate_entropy(df, 
                                  reference_cell_type, 
                                  target_cell_types, 
                                  radius, 
                                  feature_colname,
                                  dimension)  
  
  ## Cells in neighbourhood ----------
  result[["cells_in_neighbourhood"]] <- entropy_df[ , c("ref_cell_id", target_cell_types)]
  
  ## Cells in neighbourhood proportion ----------
  result[["cells_in_neighbourhood_proportion"]] <- entropy_df[ , c("ref_cell_id", target_cell_types, paste(target_cell_types, "_prop", sep = ""))]
  
  ## Entropy --------------
  result[["entropy"]] <- entropy_df
  
  
  ## These metrics focus on a particular cell type 
  for (target_cell_type in target_cell_types) {
    mixing_score_df <- data.frame(matrix(nrow = 1, ncol = length(mixing_score_df_colnames)))
    colnames(mixing_score_df) <- mixing_score_df_colnames
    mixing_score_df$ref_cell_type <- reference_cell_type
    
    cross_K_df <- data.frame(matrix(nrow = 1, ncol = length(cross_K_df_colnames)))
    colnames(cross_K_df) <- cross_K_df_colnames
    cross_K_df$ref_cell_type <- reference_cell_type
    
    ## Mixing score -----------------
    # No need to fill in mixing_score_df if the reference and target cell is the same
    if (reference_cell_type != target_cell_type) {
      mixing_score_df$tar_cell_type <- target_cell_type
      mixing_score_df$n_ref_cells <- sum(df[[feature_colname]] == reference_cell_type)
      mixing_score_df$n_tar_cells <- sum(df[[feature_colname]] == target_cell_type)
      mixing_score_df$n_ref_tar_interactions <- sum(entropy_df[[target_cell_type]])
      mixing_score_df$n_ref_ref_interactions <- sum(entropy_df[[reference_cell_type]])
      mixing_score_df$mixing_score <- mixing_score_df$n_ref_tar_interactions / (0.5 * mixing_score_df$n_ref_ref_interactions)
      mixing_score_df$normalised_mixing_score <- 0.5 * mixing_score_df$mixing_score * mixing_score_df$n_ref_cells / mixing_score_df$n_tar_cell
      if (is.infinite(mixing_score_df$mixing_score)) mixing_score_df$mixing_score <- NA
      if (is.infinite(mixing_score_df$normalised_mixing_score)) mixing_score_df$normalised_mixing_score <- NA
      result[["mixing_score"]][[target_cell_type]] <- mixing_score_df
    }
    
    ## Cross_K ---------------------
    cross_K_df$tar_cell_type <- target_cell_type
    if (dimension == "3D") {
      cross_K_df$observed_cross_K <- (((volume * sum(entropy_df[[target_cell_type]])) / sum(df[[feature_colname]] == reference_cell_type)) / sum(df[[feature_colname]] == target_cell_type))
      cross_K_df$expected_cross_K <- (4/3) * pi * radius^3
    }
    else if (dimension == "2D") {
      cross_K_df$observed_cross_K <- (((area * sum(entropy_df[[target_cell_type]])) / sum(df[[feature_colname]] == reference_cell_type)) / sum(df[[feature_colname]] == target_cell_type))
      cross_K_df$expected_cross_K <- pi * radius^2
    }
    
    cross_K_df$cross_K_ratio <- cross_K_df$observed_cross_K / cross_K_df$expected_cross_K
    result[["cross_K"]][[target_cell_type]] <- cross_K_df
  }
  
  return(result)
}




calculate_entropy <- function(df,
                              reference_cell_type,
                              target_cell_types,
                              radius,
                              feature_colname = "Cell.Type",
                              dimension) {
  
  # Check
  if (length(target_cell_types) < 2) stop("Need at least two target cell types")
  
  ## Users should ensure include the reference_cell_type as one of the target_cell_types
  cells_in_neighbourhood_proportion_df <- calculate_cells_in_neighbourhood_proportions(df,
                                                                                       reference_cell_type,
                                                                                       target_cell_types,
                                                                                       radius,
                                                                                       feature_colname,
                                                                                       dimension)
  
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





calculate_cells_in_neighbourhood_proportions <- function(df, 
                                                         reference_cell_type, 
                                                         target_cell_types, 
                                                         radius, 
                                                         feature_colname = "Cell.Type",
                                                         dimension) {
  
  ## Get cells in neighbourhood df
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood(df,
                                                                reference_cell_type,
                                                                target_cell_types,
                                                                radius,
                                                                feature_colname,
                                                                dimension)
  
  if (is.null(cells_in_neighbourhood_df)) return(NULL)
  
  ## Get total number of target cells for each row (first column is the reference cell id column, so we exclude it)
  cells_in_neighbourhood_df$total <- apply(cells_in_neighbourhood_df[ , c(-1)], 1, sum)
  
  cells_in_neighbourhood_df[ , paste(target_cell_types, "_prop", sep = "")] <- cells_in_neighbourhood_df[ , target_cell_types] / cells_in_neighbourhood_df$total
  
  return(cells_in_neighbourhood_df)
}



calculate_cells_in_neighbourhood <- function(df, 
                                             reference_cell_type, 
                                             target_cell_types, 
                                             radius, 
                                             feature_colname = "Cell.Type",
                                             dimension) {
  
  if (is.null(df[[feature_colname]])) stop(paste("No column called", feature_colname, "found in df object"))
  
  if (is.null(df[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your df")
    df$Cell.ID <- paste("Cell", seq(nrow(df)), sep = "_")
  }  
  
  ## For reference_cell_type, check it is found in the df object
  if (!(reference_cell_type %in% df[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the df object"))
    return(NULL)
  }
  
  ## For target_cell_types, check they are found in the df object
  unknown_cell_types <- setdiff(target_cell_types, df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
  }
  
  # Check if radius is numeric
  if (!is.numeric(radius)) {
    stop(paste(radius, " is not of type 'numeric'"))
  }
  
  # Get df coords
  if (dimension == "3D") {
    df_coords <- df[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]  
  }
  else if (dimension == "2D") {
    df_coords <- df[ , c("Cell.X.Position", "Cell.Y.Position")]
  }
  else {
    stop("Invalid dimension. Choose either 3D or 2D")
  }
  
  
  # Get reference_cell_type coords
  reference_cell_type_coords <- df_coords[df[[feature_colname]] == reference_cell_type, ]
  
  result <- data.frame(matrix(nrow = nrow(reference_cell_type_coords), ncol = 0))
  
  for (target_cell_type in target_cell_types) {
    
    if (sum(df[[feature_colname]] == target_cell_type) == 0) {
      result[[target_cell_type]] <- NA
      next
    }
    
    ## Get target_cell_type coords
    target_cell_type_coords <- df_coords[df[[feature_colname]] == target_cell_type, ]
    
    ## Determine number of target cells dfcified distance for each reference cell
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
  
  result <- data.frame(ref_cell_id = df$Cell.ID[df[[feature_colname]] == reference_cell_type], result)
  
  return(result)
}




get_df_grid_metrics <- function(df, 
                                n_splits, 
                                feature_colname = "Cell.Type",
                                dimension) {
  
  if (is.null(df[[feature_colname]])) stop(paste("No column called", feature_colname, "found in df object"))
  
  # Check if n_splits is numeric
  if (!is.numeric(n_splits)) {
    stop(paste(n_splits, " n_splits is not of type 'numeric'"))
  }
  
  
  if (dimension == "3D") {
    df_coords <- df[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]  
    
    
    ## Get dimensions of the window
    min_x <- min(df_coords[ , "Cell.X.Position"])
    min_y <- min(df_coords[ , "Cell.Y.Position"])
    min_z <- min(df_coords[ , "Cell.Z.Position"])
    
    max_x <- max(df_coords[ , "Cell.X.Position"])
    max_y <- max(df_coords[ , "Cell.Y.Position"])
    max_z <- max(df_coords[ , "Cell.Z.Position"])
    
    length <- round(max_x - min_x)
    width  <- round(max_y - min_y)
    height <- round(max_z - min_z)
    
    ## Get distance of row, col and lay
    d_row <- length / n_splits
    d_col <- width / n_splits
    d_lay <- height / n_splits
    
    # Shift df_coords so they begin at the origin
    df_coords[, "Cell.X.Position"] <- df_coords[, "Cell.X.Position"] - min_x
    df_coords[, "Cell.Y.Position"] <- df_coords[, "Cell.Y.Position"] - min_y
    df_coords[, "Cell.Z.Position"] <- df_coords[, "Cell.Z.Position"] - min_z
    
    ## Figure out which 'grid prism number' each cell is inside
    df$grid_prism_num <- floor(df_coords[ , "Cell.X.Position"] / d_row) +
      floor(df_coords[ , "Cell.Y.Position"] / d_col) * n_splits + 
      floor(df_coords[ , "Cell.Z.Position"] / d_lay) * n_splits^2 + 1
    
    ## Determine the cell types found in each grid prism
    n_grid_prisms <- n_splits^3
    grid_prism_cell_matrix <- as.data.frame.matrix(table(df[[feature_colname]], factor(df$grid_prism_num, levels = seq(n_grid_prisms))))
    grid_prism_cell_matrix <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                         t(grid_prism_cell_matrix))
    
    ## Determine centre coordinates of each grid prism
    grid_prism_coordinates <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                         x_coord = ((seq(n_grid_prisms) - 1) %% n_splits + 0.5) * d_row + round(min_x),
                                         y_coord = (floor(((seq(n_grid_prisms) - 1) %% (n_splits)^2) / n_splits) + 0.5) * d_col + round(min_y),
                                         z_coord = (floor((seq(n_grid_prisms) - 1) / (n_splits^2)) + 0.5) * d_lay + round(min_z))
    
    grid_prism_data <- list("grid_prism_cell_matrix" = grid_prism_cell_matrix,
                            "grid_prism_coordinates" = grid_prism_coordinates)
  }
  
  else if (dimension == "2D") {
    df_coords <- df[ , c("Cell.X.Position", "Cell.Y.Position")]  
    
    
    ## Get dimensions of the window
    min_x <- min(df_coords[ , "Cell.X.Position"])
    min_y <- min(df_coords[ , "Cell.Y.Position"])
    
    max_x <- max(df_coords[ , "Cell.X.Position"])
    max_y <- max(df_coords[ , "Cell.Y.Position"])
    
    length <- round(max_x - min_x)
    width  <- round(max_y - min_y)
    
    ## Get distance of row, col
    d_row <- length / n_splits
    d_col <- width / n_splits
    
    # Shift df_coords so they begin at the origin
    df_coords[, "Cell.X.Position"] <- df_coords[, "Cell.X.Position"] - min_x
    df_coords[, "Cell.Y.Position"] <- df_coords[, "Cell.Y.Position"] - min_y
    
    ## Figure out which 'grid prism number' each cell is inside
    df$grid_prism_num <- floor(df_coords[ , "Cell.X.Position"] / d_row) +
      floor(df_coords[ , "Cell.Y.Position"] / d_col) * n_splits
    
    ## Determine the cell types found in each grid prism
    n_grid_prisms <- n_splits^2
    grid_prism_cell_matrix <- as.data.frame.matrix(table(df[[feature_colname]], factor(df$grid_prism_num, levels = seq(n_grid_prisms))))
    grid_prism_cell_matrix <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                         t(grid_prism_cell_matrix))
    
    ## Determine centre coordinates of each grid prism
    grid_prism_coordinates <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                         x_coord = ((seq(n_grid_prisms) - 1) %% n_splits + 0.5) * d_row + round(min_x),
                                         y_coord = (floor(((seq(n_grid_prisms) - 1) %% (n_splits)^2) / n_splits) + 0.5) * d_col + round(min_y))
    
    grid_prism_data <- list("grid_prism_cell_matrix" = grid_prism_cell_matrix,
                            "grid_prism_coordinates" = grid_prism_coordinates)
  }
  else {
    stop("Invalid dimension. Choose either 3D or 2D")
  }
  return(grid_prism_data)
}



calculate_cell_proportion_grid_metrics <- function(df, 
                                                   n_splits,
                                                   reference_cell_types,
                                                   target_cell_types,
                                                   feature_colname = "Cell.Type",
                                                   dimension) {
  
  if (is.null(df[[feature_colname]])) stop(paste("No column called", feature_colname, "found in df object"))
  
  ## Check reference_cell_types are found in the df object
  unknown_cell_types <- setdiff(reference_cell_types, df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in reference_cell_types are not found in the df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  ## Check target_cell_types are found in the df object
  unknown_cell_types <- setdiff(target_cell_types, df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  # Check if there is intersection between reference_cell_types and target_cell_types
  if (length(intersect(reference_cell_types, target_cell_types)) > 0) {
    stop("Cannot have same cells in both reference_cell_types and target_cell_types")
  }
  
  # Add grid metrics to df
  grid_prism_data <- get_df_grid_metrics(df, n_splits, feature_colname, dimension)
  
  # Get grid_prism_cell_matrix from df
  grid_prism_cell_matrix <- grid_prism_data$grid_prism_cell_matrix
  
  ## Define data frame which contains all results
  if (dimension == "3D") {
    n_grid_prisms <- n_splits^3  
  }
  else if (dimension == "2D") {
    n_grid_prisms <- n_splits^2
  }
  else {
    stop("Invalid dimension. Choose either 3D or 2D")
  }
  
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
  result <- cbind(result, grid_prism_data$grid_prism_coordinates)
  
  return(result)
}


calculate_entropy_grid_metrics <- function(df, 
                                           n_splits,
                                           cell_types_of_interest,
                                           feature_colname = "Cell.Type",
                                           dimension) {
  
  if (is.null(df[[feature_colname]])) stop(paste("No column called", feature_colname, "found in df object"))
  
  ## If cell types have been chosen, check they are found in the df object
  unknown_cell_types <- setdiff(cell_types_of_interest, unique(df[[feature_colname]]))
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in cell_types_of_interest are not found in the df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  
  # Add grid metrics to df
  grid_prism_data <- get_df_grid_metrics(df, n_splits, feature_colname, dimension)
  
  # Get grid_prism_cell_matrix from df
  grid_prism_cell_matrix <- grid_prism_data$grid_prism_cell_matrix
  
  ## Define data frame which contains all results
  if (dimension == "3D") {
    n_grid_prisms <- n_splits^3  
  }
  else if (dimension == "2D") {
    n_grid_prisms <- n_splits^2
  }
  else {
    stop("Invalid dimension. Choose either 3D or 2D")
  }
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
  result <- cbind(result, grid_prism_data$grid_prism_coordinates)
  
  return(result)
}



calculate_spatial_autocorrelation <- function(grid_metrics,
                                              metric_colname,
                                              weight_method = 0.1,
                                              dimension) {
  
  
  ## Get number of grid prisms
  n_grid_prisms <- nrow(grid_metrics)
  
  if (dimension == "3D") {
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
    else {
      stop(paste(weight_method, " weight_method is not an appropriate method"))
    }
  }
  else if (dimension == "2D") {
    ## Get splitting number (should be the square root of n_grid_prisms)
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
    ## Use Rook method: adjacent points get a weight of 1, otherwise, weight of 0
    ## Adjacent points are within 1 unit apart.
    else if (weight_method == "rook") {
      weight_matrix <- ifelse(weight_matrix > 1, 0, 1)  
    }
    else if (weight_method == "queen") {
      weight_matrix <- ifelse(weight_matrix > sqrt(2), 0, 1)  
    }
    ## If a number (x) between 0 and 1 is supplied, set a threshold to be x quantile value of c(weight_matrix)
    ## Grid prisms within this specified threshold have a weight of 1, otherwise, weight of 0
    else if (as.numeric(weight_method) && 0 < weight_method && weight_method < 1) {
      threshold <- quantile(c(weight_matrix), weight_method)
      weight_matrix <- ifelse(weight_matrix > threshold, 0, 1)
    }
    else {
      stop(paste(weight_method, " weight_method is not an appropriate method"))
    }
  }
  else {
    stop("Invalid dimension. Choose either 3D or 2D")
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


calculate_prevalence_gradient <- function(grid_metrics,
                                          metric_colname) {
  
  # Thresholds range from 0 to 1
  thresholds <- seq(0.01, 1, 0.01)
  
  # Define result
  result <- data.frame(threshold = thresholds)
  
  # Get prevalences for each threshold
  result$prevalence <- sapply(thresholds, function(threshold) { 
    calculate_prevalence(grid_metrics, metric_colname, threshold) 
  })
  
  return(result)
}


calculate_prevalence <- function(grid_metrics,
                                 metric_colname,
                                 threshold,
                                 above = TRUE) {
  
  ## Exclude rows with NA values
  grid_metrics <- grid_metrics[!is.na(grid_metrics[[metric_colname]]), ]
  
  if (above) {
    p <- sum(grid_metrics[[metric_colname]] >= threshold) / nrow(grid_metrics) * 100
  }
  else {
    p <- sum(grid_metrics[[metric_colname]] < threshold) / nrow(grid_metrics) * 100    
  }
  
  return(p)
}

### Analyse 2D and 3D data -------
slice_z_coords <- unique(df3D$Cell.Z.Position)
feature_colname <- "Cell.Type.Generic"
df3D[df3D$Cell.Type.Generic == "Tumor", "Cell.Type.Generic"] <- "Tumour"


for (i in seq(n_slices + 1)) {
  
  # i represents the current slice index
  if (i != n_slices + 1) {
    df <- df3D[df3D$Cell.Z.Position == slice_z_coords[i], ]
    dimension <- "2D"
    slice <- i
  }
  # if i == n_slices + 1, analyse in 3D instead
  else {
    df <- df3D
    dimension <- "3D"
    slice <- 0
  }
  
  minimum_distance_data <- calculate_minimum_distances_between_cell_types(df,
                                                                          cell_types,
                                                                          feature_colname,
                                                                          dimension = dimension)
  
  minimum_distance_data_summary <- summarise_distances_between_cell_types(minimum_distance_data)
  ## Fill in 4 rows at a time for AMD df (as we have Tumour/Tumour, Tumour/Immune, Immune/Tumour, Immune/IMmune)
  index <- 4 * (i - 1) + 1 # index is 1, 5, 9, 13
  
  metric_df_list[["AMD"]][index:(index + 3), "slice"] <- slice
  metric_df_list[["AMD"]][index:(index + 3), "reference"] <- minimum_distance_data_summary$reference
  metric_df_list[["AMD"]][index:(index + 3), "target"] <- minimum_distance_data_summary$target
  metric_df_list[["AMD"]][index:(index + 3), "AMD"] <- minimum_distance_data_summary$mean
  
  
  
  index1 <- 2 * (i - 1) + 1 # index1 is 1, 3, 5, ...
  index2 <- 4 * (i - 1) + 1 # index2 is 1, 5, 9, 13...
  for (reference_cell_type in cell_types) {
    gradient_data <- calculate_all_gradient_cc_metrics(df,
                                                       reference_cell_type,
                                                       cell_types,
                                                       radii,
                                                       feature_colname,
                                                       dimension = dimension)
    
    target_cell_type <- setdiff(cell_types, reference_cell_type)
    
    if (!is.null(gradient_data)) {
      metric_df_list[["MS"]][index1, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$mixing_score
      metric_df_list[["NMS"]][index1, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$normalised_mixing_score
      metric_df_list[["ACINP"]][index1, radii_colnames] <- gradient_data[["cells_in_neighbourhood_proportion"]][["Immune"]]
      metric_df_list[["AE"]][index1, radii_colnames] <- gradient_data[["entropy"]]$entropy
    }
    else {
      metric_df_list[["MS"]][index1, radii_colnames] <- NA
      metric_df_list[["NMS"]][index1, radii_colnames] <- NA
      metric_df_list[["ACINP"]][index1, radii_colnames] <- NA
      metric_df_list[["AE"]][index1, radii_colnames] <- NA
    }
    
    metric_df_list[["MS"]][index1, "slice"] <- slice
    metric_df_list[["MS"]][index1, c("reference", "target")] <- c(reference_cell_type, target_cell_type)
    
    metric_df_list[["NMS"]][index1, "slice"] <- slice
    metric_df_list[["NMS"]][index1, c("reference", "target")] <- c(reference_cell_type, target_cell_type)
    
    metric_df_list[["ACINP"]][index1, "slice"] <- slice
    metric_df_list[["ACINP"]][index1, c("reference", "target")] <- c(reference_cell_type, "Immune")
    
    metric_df_list[["AE"]][index1, "slice"] <- slice
    metric_df_list[["AE"]][index1, c("reference", "target")] <- c(reference_cell_type, "Tumour,Immune")

    index1 <- index1 + 1
    
    for (target_cell_type in cell_types) {
      ## Calculate ACIN and CKR as target cell type can also be the reference cell type
      
      # ACIN
      metric_df_list[["ACIN"]][index2, "slice"] <- slice
      metric_df_list[["ACIN"]][index2, c("reference", "target")] <- c(reference_cell_type, target_cell_type)

      # CKR
      metric_df_list[["CKR"]][index2, "slice"] <- slice
      metric_df_list[["CKR"]][index2, c("reference", "target")] <- c(reference_cell_type, target_cell_type)

      
      if (!is.null(gradient_data)) {
        metric_df_list[["ACIN"]][index2, radii_colnames] <- gradient_data[["cells_in_neighbourhood"]][[target_cell_type]]
        metric_df_list[["CKR"]][index2, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]]$cross_K_ratio
      }
      else {
        metric_df_list[["ACIN"]][index2, radii_colnames] <- NA
        metric_df_list[["CKR"]][index2, radii_colnames] <- NA
      }
      
      index2 <- index2 + 1
    }
  }
  
  
  # Get proportion grid metrics
  for (j in seq_len(nrow(prop_cell_types))) {
    proportion_grid_metrics <- calculate_cell_proportion_grid_metrics(df, 
                                                                      n_splits,
                                                                      strsplit(prop_cell_types$ref[j], ",")[[1]], 
                                                                      strsplit(prop_cell_types$tar[j], ",")[[1]],
                                                                      feature_colname,
                                                                      dimension = dimension)
    
    
    if (!is.null(proportion_grid_metrics)) {
      proportion_SAC <- calculate_spatial_autocorrelation(proportion_grid_metrics, 
                                                          "proportion",
                                                          weight_method = 0.1,
                                                          dimension = dimension)
      
      proportion_prevalence_df <- calculate_prevalence_gradient(proportion_grid_metrics,
                                                                "proportion")
    }
    else {
      proportion_SAC <- NA
      proportion_prevalence_df <- data.frame(threshold = seq(0.01, 1, 0.01), prevalence = NA)
    }

    index <- nrow(prop_cell_types) * (i - 1) + j
    metric_df_list[["prop_SAC"]][index, "slice"] <- slice
    metric_df_list[["prop_SAC"]][index, c("reference", "target")] <- c(prop_cell_types$ref[j], prop_cell_types$tar[j])
    metric_df_list[["prop_SAC"]][index, "prop_SAC"] <- proportion_SAC
    
    metric_df_list[["prop_prevalence"]][index, "slice"] <- slice
    metric_df_list[["prop_prevalence"]][index, c("reference", "target")] <- c(prop_cell_types$ref[j], prop_cell_types$tar[j])
    metric_df_list[["prop_prevalence"]][index, thresholds_colnames] <- proportion_prevalence_df$prevalence
  }
  
  # Get entropy grid metrics
  for (j in seq_len(nrow(entropy_cell_types))) {
    entropy_grid_metrics <- calculate_entropy_grid_metrics(df, 
                                                           n_splits,
                                                           strsplit(entropy_cell_types$cell_types[j], ",")[[1]],
                                                           feature_colname,
                                                           dimension = dimension)
    
    if (!is.null(entropy_grid_metrics)) {
      entropy_SAC <- calculate_spatial_autocorrelation(entropy_grid_metrics, 
                                                       "entropy",
                                                       weight_method = 0.1,
                                                       dimension = dimension)
      
      entropy_prevalence_df <- calculate_prevalence_gradient(entropy_grid_metrics,
                                                             "entropy")
      
    }
    else {
      entropy_SAC <- NA
      entropy_prevalence_df <- data.frame(threshold = seq(0.01, 1, 0.01), prevalence = NA)
    }
    
    index <- nrow(entropy_cell_types) * (i - 1) + j
    metric_df_list[["entropy_SAC"]][index, "slice"] <- slice
    metric_df_list[["entropy_SAC"]][index, c("cell_types")] <- c(entropy_cell_types$cell_types[j])
    metric_df_list[["entropy_SAC"]][index, "entropy_SAC"] <- entropy_SAC
    
    metric_df_list[["entropy_prevalence"]][index, "slice"] <- slice
    metric_df_list[["entropy_prevalence"]][index, c("cell_types")] <- c(entropy_cell_types$cell_types[j])
    metric_df_list[["entropy_prevalence"]][index, thresholds_colnames] <- entropy_prevalence_df$prevalence
  }  
}


setwd("~/R/Lin et al - human colorectal cancer/CRC1_data_final")
saveRDS(metric_df_list, "metric_df_list.RDS")

### Plot analysis of 2D and 3D data -----
setwd("~/R/Lin et al - human colorectal cancer/CRC1_data_final")
metric_df_list <- readRDS("metric_df_list.RDS")

get_gradient <- function(metric) {
  if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR")) {
    return("radius")
  }
  else if (metric %in% c("prop_prevalence", "entropy_prevalence")) {
    return("threshold")  
  }
  else {
    stop("Invalid metric. Must be gradient-based")
  }
}


## Turn gradient radii metrics into AUC and add to metric_df list
get_AUC_for_radii_gradient_metrics <- function(y) {
  x <- radii
  h <- diff(x)[1]
  n <- length(x)
  
  AUC <- (h / 2) * (y[1] + 2 * sum(y[2:(n - 1)]) + y[n])
  
  return(AUC)
}


radii <- seq(20, 100, 10)
radii_colnames <- paste("r", radii, sep = "")

gradient_radii_metrics <- c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR")


for (metric in gradient_radii_metrics) {
  metric_AUC_name <- paste(metric, "AUC", sep = "_")
  
  if (metric %in% c("MS", "NMS", "ACIN", "CKR")) {
    subset_colnames <- c("slice", "reference", "target", metric_AUC_name)
  }
  else {
    subset_colnames <- c("slice", "reference", metric_AUC_name)
  }
  
  df <- metric_df_list[[metric]]
  df[[metric_AUC_name]] <- apply(df[ , radii_colnames], 1, get_AUC_for_radii_gradient_metrics)
  metric_df_list[[metric_AUC_name]] <- df
  
}

## Turn threshold radii metrics into AUC and add to metric_df list
thresholds <- seq(0.01, 1, 0.01)
threshold_colnames <- paste("t", thresholds, sep = "")

# prop_AUC 3D
prop_prevalence_df <- metric_df_list[["prop_prevalence"]]
prop_prevalence_df$prop_AUC <- apply(prop_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
prop_AUC_df <- prop_prevalence_df[ , c("slice", "reference", "target", "prop_AUC")]
metric_df_list[["prop_AUC"]] <- prop_AUC_df

# entropy_AUC 3D
entropy_prevalence_df <- metric_df_list[["entropy_prevalence"]]
entropy_prevalence_df$entropy_AUC <- apply(entropy_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
entropy_AUC_df <- entropy_prevalence_df[ , c("slice", "cell_types", "entropy_AUC")]
metric_df_list[["entropy_AUC"]] <- entropy_AUC_df


## Functions to plot
# Utility function to get metric cell types
get_metric_cell_types <- function(metric) {
  # Get metric_cell_types
  if (metric %in% c("AMD", "ACIN", "CKR", "ACIN_AUC", "CKR_AUC")) {
    metric_cell_types <- data.frame(ref = c("Tumour"), tar = c("Immune"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("MS", "NMS", "MS_AUC", "NMS_AUC")) {
    metric_cell_types <- data.frame(ref = c("Tumour"), tar = c("Immune"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("ACINP", "ACINP_AUC")) {
    metric_cell_types <- data.frame(ref = c("Tumour"), tar = c("Tumour"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("AE", "AE_AUC")) {
    metric_cell_types <- data.frame(ref = c("Tumour"), tar = c("Tumour,Immune"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("prop_SAC", "prop_prevalence", "prop_AUC")) {
    metric_cell_types <- data.frame(ref = c("Tumour"), tar = c("Immune"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("entropy_SAC", "entropy_prevalence", "entropy_AUC")) {
    metric_cell_types <- data.frame(cell_types = c("Tumour,Immune"))
  }
  else {
    stop("metric not found")
  }
  return(metric_cell_types)
}

# Utility function to subset metric_df
subset_metric_df <- function(metric,
                             metric_df,
                             metric_cell_types,
                             index) {
  
  if (metric %in% c("AMD", "ACIN", "CKR", "MS", "NMS", "ACIN_AUC", "CKR_AUC", "MS_AUC", "NMS_AUC", "prop_SAC", "prop_prevalence", "prop_AUC")) {
    metric_df_subset <- metric_df[metric_df$reference == metric_cell_types[index, "ref"] & metric_df$target == metric_cell_types[index, "tar"], ] 
  }
  else if (metric %in% c("ACINP", "AE", "ACINP_AUC", "AE_AUC")) {
    metric_df_subset <- metric_df[metric_df$reference == metric_cell_types[index, "ref"], ] 
  }
  else if (metric %in% c("entropy_SAC", "entropy_prevalence", "entropy_AUC")) {
    metric_df_subset <- metric_df[metric_df$cell_types == metric_cell_types[index, "cell_types"], ]
  }
  else {
    stop("metric not found")
  }
  
  return(metric_df_subset)
}



plot_3D_vs_2D <- function(metric_df_list,
                          metric) {
  
  # Get metric_df for current metric
  metric_df <- metric_df_list[[metric]]
  
  # Get metric cell types for current metric (should only be one set/ one row)
  metric_cell_types <- get_metric_cell_types(metric)
  
  # Subset metric_df
  metric_df_subset <- subset_metric_df(metric,
                                       metric_df,
                                       metric_cell_types,
                                       1) # Always first row
  
  # Change and further subset columns of metric_df_subset
  colnames(metric_df_subset)[colnames(metric_df_subset) == metric] <- "value"
  metric_df_subset$metric <- metric
  metric_df_subset <- metric_df_subset[ , c("slice", "value", "metric")]
  
  metric_df_subset$dummy <- "dummy"
  metric_df_subset$metric <- factor(metric_df_subset$metric, metrics)
  
  # Create the dot plot, highlighting the maximum slice points with a star shape and using facets
  fig <- ggplot(metric_df_subset[metric_df_subset$slice != 0, ], aes(x = dummy, y = value)) +
    geom_jitter(width = 0.2, height = 0, size = 1.5) +
    geom_point(data = metric_df_subset[metric_df_subset$slice == 0, ], color = "red", shape = 8, size = 6) +
    labs(x = "", y = metric) +
    theme_bw()  +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  
  return(fig)
  
}


plot_3D_vs_error <- function(metric_df_list,
                             metrics) {
  
  # Get metric_df for current metric
  metric_df <- metric_df_list[[metric]]
  
  # Get metric cell types for current metric (should only be one set/ one row)
  metric_cell_types <- get_metric_cell_types(metric)
  
  # Subset metric_df
  metric_df_subset <- subset_metric_df(metric,
                                       metric_df,
                                       metric_cell_types,
                                       1) # Always first row
  
  # Change and further subset columns of metric_df_subset
  colnames(metric_df_subset)[colnames(metric_df_subset) == metric] <- "value"
  metric_df_subset$metric <- metric
  metric_df_subset <- metric_df_subset[ , c("slice", "value", "metric")]
  
  metric_df_subset$dummy <- "dummy"
  metric_df_subset$metric <- factor(metric_df_subset$metric, metrics)
  
  # Calculate error for each slice, and remove 3D row
  value_3D <- metric_df_subset[["value"]][metric_df_subset[["slice"]] == 0]
  metric_df_subset[["value"]] <- ((metric_df_subset[["value"]] - value_3D) / value_3D) * 100
  metric_df_subset <- metric_df_subset[metric_df_subset[["slice"]] != 0, ]
  
  # Create the dot plot
  fig <- ggplot(metric_df_subset, aes(x = dummy, y = value)) +
    geom_point(data = data.frame(x = "dummy", y = 0), aes(x, y), size = 0) + # Ensures plot shows y = 0
    geom_jitter(width = 0.2, height = 0, size = 1.5) +
    geom_abline(intercept = 0, slope = 0, color = "red", linetype = "longdash") +
    labs(x = "", y = paste(metric, "error (%)")) +
    theme_bw()  +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  
  return(fig)
  
}




## Get plot
metrics <- c("AMD", "ACIN_AUC", "ACINP_AUC", "AE_AUC", "MS_AUC", "NMS_AUC", "CKR_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")

plot3D_vs_2D_metric_list <- list()
plot3D_vs_error_metric_list <- list()


for (metric in metrics) {
  plot3D_vs_2D_metric_list[[metric]] <- plot_3D_vs_2D(metric_df_list,
                                                      metric)
  plot3D_vs_error_metric_list[[metric]] <- plot_3D_vs_error(metric_df_list,
                                                         metric)
}

plots3D_vs_2D <- plot_grid(plotlist = plot3D_vs_2D_metric_list,
                           nrow = 3,
                           ncol = 4,
                           labels = LETTERS[1:13])

plots3D_vs_error <- plot_grid(plotlist = plot3D_vs_error_metric_list,
                              nrow = 3,
                              ncol = 4,
                              labels = LETTERS[1:13])

methods::show(plots3D_vs_2D)
methods::show(plots3D_vs_error)


setwd("~/R/Lin et al - human colorectal cancer/CRC1_data_final")
pdf("lin_et_al_plots.pdf", width = 10, height = 8)

print(plots3D_vs_2D)
print(plots3D_vs_error)

dev.off()





### Get error_values ----
plot_3D_vs_error_values <- function(metric_df_list,
                                    metrics) {
  
  # Get metric_df for current metric
  metric_df <- metric_df_list[[metric]]
  
  # Get metric cell types for current metric (should only be one set/ one row)
  metric_cell_types <- get_metric_cell_types(metric)
  
  # Subset metric_df
  metric_df_subset <- subset_metric_df(metric,
                                       metric_df,
                                       metric_cell_types,
                                       1) # Always first row
  
  # Change and further subset columns of metric_df_subset
  colnames(metric_df_subset)[colnames(metric_df_subset) == metric] <- "value"
  metric_df_subset$metric <- metric
  metric_df_subset <- metric_df_subset[ , c("slice", "value", "metric")]
  
  metric_df_subset$dummy <- "dummy"
  metric_df_subset$metric <- factor(metric_df_subset$metric, metrics)
  
  # Calculate error for each slice, and remove 3D row
  value_3D <- metric_df_subset[["value"]][metric_df_subset[["slice"]] == 0]
  metric_df_subset[["value"]] <- ((metric_df_subset[["value"]] - value_3D) / value_3D) * 100
  metric_df_subset <- metric_df_subset[metric_df_subset[["slice"]] != 0, ]
  
  text_mean <- paste("mean:", round(mean(metric_df_subset[["value"]], na.rm = T), 2))
  
  fig <- ggplot() +
    theme(axis.ticks = element_blank(), axis.text = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"),
          panel.grid = element_blank()) +
    labs(x= "", y = "") +
    annotate("text", x = 0.5, y = 0.5, label = text_mean, size = 4, hjust = 0.5, vjust = -2)
  
  
  # Outliers (underestimate heavily)
  if (metric %in% c("AMD", "entropy_AUC")) {
    text_outlier <- paste("outlier:", round(min(metric_df_subset[["value"]], na.rm = T), 2))
    fig <- fig +
      annotate("text", x = 0.5, y = 0.5, label = text_outlier, size = 4, hjust = 0.5, vjust = 0)
    
    text_min <- paste("min:", round(sort(metric_df_subset[["value"]])[2], 2))
    text_max <- paste("max:", round(max(metric_df_subset[["value"]], na.rm = T), 2))
    fig <- fig +
      annotate("text", x = 0.5, y = 0.5, label = text_min, size = 4, hjust = 0.5, vjust = 2) +
      annotate("text", x = 0.5, y = 0.5, label = text_max, size = 4, hjust = 0.5, vjust = 4)
  }
  # Outliers (overestimate heavily)
  if (metric %in% c("ACINP_AUC", "AE_AUC", "MS_AUC", "prop_AUC")) {
    text_outlier <- paste("outlier:", round(max(metric_df_subset[["value"]], na.rm = T), 2))
    fig <- fig +
      annotate("text", x = 0.5, y = 0.5, label = text_outlier, size = 4, hjust = 0.5, vjust = 0)
    
    text_min <- paste("min:", round(min(metric_df_subset[["value"]], na.rm = T), 2))
    text_max <- paste("max:", round(sort(metric_df_subset[["value"]][2], decreasing = T), 2))
    fig <- fig +
      annotate("text", x = 0.5, y = 0.5, label = text_min, size = 4, hjust = 0.5, vjust = 2) +
      annotate("text", x = 0.5, y = 0.5, label = text_max, size = 4, hjust = 0.5, vjust = 4)
  }
  
  
  return(fig)
  
}

metrics <- c("AMD", "ACIN_AUC", "ACINP_AUC", "AE_AUC", "MS_AUC", "NMS_AUC", "CKR_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")
plot3D_vs_error_metric_list <- list()

for (metric in metrics) {
  plot3D_vs_error_metric_list[[metric]] <- plot_3D_vs_error_values(metric_df_list,
                                                                   metric)
}

plots3D_vs_error <- plot_grid(plotlist = plot3D_vs_error_metric_list,
                              nrow = 3,
                              ncol = 4,
                              labels = LETTERS[1:13])

setwd("~/R/Lin et al - human colorectal cancer/CRC1_data_final")
pdf("lin_et_al_plots_error_values.pdf", width = 10, height = 8)

print(plots3D_vs_error)

dev.off()