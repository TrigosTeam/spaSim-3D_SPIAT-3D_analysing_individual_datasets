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

setwd("~/R/spaSim-3D/scripts/simulations and analysis S2/S2 data")
spes_metadata <- readRDS("spes_metadata.RDS")

# Get number of spe groups
n_spe_groups <- length(spes_metadata)

# Get number of spes in each group (should be the same)
n_spes <- length(spes_metadata[[1]])

# Define APD and AMD data frames as well as constants
cell_types <- c("A", "B")
APD_pairs <- c("A/A", "A/B", "B/B")
AMD_pairs <- c("A/A", "A/B", "B/A", "B/B")

APD_colnames <- c("spe", "pair", "APD")
APD_df <- data.frame(matrix(nrow = n_spes * length(APD_pairs), ncol = length(APD_colnames)))
colnames(APD_df) <- APD_colnames

AMD_colnames <- c("spe", "reference", "target", "AMD")
AMD_df <- data.frame(matrix(nrow = n_spes * length(AMD_pairs), ncol = 4))
colnames(AMD_df) <- AMD_colnames


# Define MS, NMS, ACIN, ACINP, CKR, AE data frames as well as constants
radii <- seq(20, 100, 10)
radii_colnames <- paste("r", radii, sep = "")

MS_colnames <- c("spe", "reference", "target", radii_colnames)
MS_df <- data.frame(matrix(nrow = n_spes * length(cell_types), ncol = length(MS_colnames)))
colnames(MS_df) <- MS_colnames

NMS_colnames <- c("spe", "reference", "target", radii_colnames)
NMS_df <- data.frame(matrix(nrow = n_spes * length(cell_types), ncol = length(NMS_colnames)))
colnames(NMS_df) <- NMS_colnames

# Only choose prop(A) as prop(A) = 1 - prop(B) always
ACINP_colnames <- c("spe", "reference", "target", radii_colnames)
ACINP_df <- data.frame(matrix(nrow = n_spes * length(cell_types), ncol = length(ACINP_colnames)))
colnames(ACINP_df) <- ACINP_colnames

# Target is always A and B together
AE_colnames <- c("spe", "reference", "target", radii_colnames)
AE_df <- data.frame(matrix(nrow = n_spes * length(cell_types), ncol = length(AE_colnames)))
colnames(AE_df) <- AE_colnames

## ACIN and CKR are twice as large
# (ref A and tar A or B) OR (ref B and tar B or A)
ACIN_colnames <- c("spe", "reference", "target", radii_colnames)
ACIN_df <- data.frame(matrix(nrow = n_spes * length(cell_types)^2, ncol = length(ACIN_colnames)))
colnames(ACIN_df) <- ACIN_colnames

# (ref A and tar A or B) OR (ref B and tar B or A)
CKR_colnames <- c("spe", "reference", "target", radii_colnames)
CKR_df <- data.frame(matrix(nrow = n_spes * length(cell_types)^2, ncol = length(CKR_colnames)))
colnames(CKR_df) <- CKR_colnames


# Define SAC and prevalence data frames as well as constants
n_splits <- 10
thresholds <- seq(0.01, 1, 0.01)
thresholds_colnames <- paste("t", thresholds, sep = "")

prop_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))

prop_SAC_df_colnames <- c("spe", "reference", "target", "prop_SAC")
prop_SAC_df <- data.frame(matrix(nrow = n_spes * nrow(prop_cell_types), ncol = length(prop_SAC_df_colnames)))
colnames(prop_SAC_df) <- prop_SAC_df_colnames

prop_prevalence_df_colnames <- c("spe", "reference", "target", thresholds_colnames)
prop_prevalence_df <- data.frame(matrix(nrow = n_spes * nrow(prop_cell_types), ncol = length(prop_prevalence_df_colnames)))
colnames(prop_prevalence_df) <- prop_prevalence_df_colnames


entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))

entropy_SAC_df_colnames <- c("spe", "cell_types", "entropy_SAC")
entropy_SAC_df <- data.frame(matrix(nrow = n_spes * nrow(entropy_cell_types), ncol = length(entropy_SAC_df_colnames)))
colnames(entropy_SAC_df) <- entropy_SAC_df_colnames

entropy_prevalence_df_colnames <- c("spe", "cell_types", thresholds_colnames)
entropy_prevalence_df <- data.frame(matrix(nrow = n_spes * nrow(entropy_cell_types), ncol = length(entropy_prevalence_df_colnames)))
colnames(entropy_prevalence_df) <- entropy_prevalence_df_colnames


# Add all to list:
metric_df_list3D <- list(APD = APD_df,
                         AMD = AMD_df,
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


# Define APD and AMD data frames as well as constants
cell_types <- c("A", "B")
APD_pairs <- c("A/A", "A/B", "B/B")
AMD_pairs <- c("A/A", "A/B", "B/A", "B/B")

slices_APD_df_colnames <- c("spe", "slice", "pair", "APD")
slices_APD_df <- data.frame(matrix(nrow = n_spes * length(APD_pairs) * n_slices, ncol = length(slices_APD_df_colnames)))
colnames(slices_APD_df) <- slices_APD_df_colnames

slices_AMD_df_colnames <- c("spe", "slice", "reference", "target", "AMD")
slices_AMD_df <- data.frame(matrix(nrow = n_spes * length(AMD_pairs) * n_slices, ncol = length(slices_AMD_df_colnames)))
colnames(slices_AMD_df) <- slices_AMD_df_colnames

# Define MS, NMS, ACIN, ACINP, CKR, AE data frames as well as constants
radii <- seq(20, 100, 10)
radii_colnames <- paste("r", radii, sep = "")

slices_MS_df_colnames <- c("spe", "slice", "reference", "target", radii_colnames)
slices_MS_df <- data.frame(matrix(nrow = n_spes * length(cell_types) * n_slices, ncol = length(slices_MS_df_colnames)))
colnames(slices_MS_df) <- slices_MS_df_colnames

slices_NMS_df_colnames <- c("spe", "slice", "reference", "target", radii_colnames)
slices_NMS_df <- data.frame(matrix(nrow = n_spes * length(cell_types) * n_slices, ncol = length(slices_MS_df_colnames)))
colnames(slices_NMS_df) <- slices_NMS_df_colnames

# Target is always A and B together
# Only choose prop(A) as prop(A) = 1 - prop(B) always
slices_ACINP_df_colnames <- c("spe", "slice", "reference", "target", radii_colnames)
slices_ACINP_df <- data.frame(matrix(nrow = n_spes * length(cell_types) * n_slices, ncol = length(slices_ACINP_df_colnames)))
colnames(slices_ACINP_df) <- slices_ACINP_df_colnames

# Target is always A and B together
slices_AE_df_colnames <- c("spe", "slice", "reference", "target", radii_colnames)
slices_AE_df <- data.frame(matrix(nrow = n_spes * length(cell_types) * n_slices, ncol = length(slices_AE_df_colnames)))
colnames(slices_AE_df) <- slices_AE_df_colnames

## ACIN and CKR are twice as large
# (ref A and tar A or B) OR (ref B and tar B or A)
slices_ACIN_df_colnames <- c("spe", "slice", "reference", "target", radii_colnames)
slices_ACIN_df <- data.frame(matrix(nrow = n_spes * length(cell_types)^2 * n_slices, ncol = length(slices_ACIN_df_colnames)))
colnames(slices_ACIN_df) <- slices_ACIN_df_colnames

# (ref A and tar A or B) OR (ref B and tar B or A)
slices_CKR_df_colnames <- c("spe", "slice", "reference", "target", radii_colnames)
slices_CKR_df <- data.frame(matrix(nrow = n_spes * length(cell_types)^2 * n_slices, ncol = length(slices_CKR_df_colnames)))
colnames(slices_CKR_df) <- slices_CKR_df_colnames

# Define SAC and prevalence data frames as well as constants
n_splits <- 10
thresholds <- seq(0.01, 1, 0.01)
thresholds_colnames <- paste("t", thresholds, sep = "")

prop_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))

slices_prop_SAC_df_colnames <- c("spe", "slice", "reference", "target", "prop_SAC")
slices_prop_SAC_df <- data.frame(matrix(nrow = n_spes * nrow(prop_cell_types) * n_slices, ncol = length(slices_prop_SAC_df_colnames)))
colnames(slices_prop_SAC_df) <- slices_prop_SAC_df_colnames

slices_prop_prevalence_df_colnames <- c("spe", "slice", "reference", "target", thresholds_colnames)
slices_prop_prevalence_df <- data.frame(matrix(nrow = n_spes * nrow(prop_cell_types) * n_slices, ncol = length(slices_prop_prevalence_df_colnames)))
colnames(slices_prop_prevalence_df) <- slices_prop_prevalence_df_colnames


entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))

slices_entropy_SAC_df_colnames <- c("spe", "slice", "cell_types", "entropy_SAC")
slices_entropy_SAC_df <- data.frame(matrix(nrow = n_spes * nrow(entropy_cell_types) * n_slices, ncol = length(slices_entropy_SAC_df_colnames)))
colnames(slices_entropy_SAC_df) <- slices_entropy_SAC_df_colnames

slices_entropy_prevalence_df_colnames <- c("spe", "slice", "cell_types", thresholds_colnames)
slices_entropy_prevalence_df <- data.frame(matrix(nrow = n_spes * nrow(entropy_cell_types) * n_slices, ncol = length(slices_entropy_prevalence_df_colnames)))
colnames(slices_entropy_prevalence_df) <- slices_entropy_prevalence_df_colnames


# Add all to list:
metric_df_list2D <- list(APD = slices_APD_df,
                         AMD = slices_AMD_df,
                         MS = slices_MS_df,
                         NMS = slices_NMS_df,
                         ACINP = slices_ACINP_df,
                         AE = slices_AE_df,
                         ACIN = slices_ACIN_df,
                         CKR = slices_CKR_df,
                         prop_SAC = slices_prop_SAC_df,
                         prop_prevalence = slices_prop_prevalence_df,
                         entropy_SAC = slices_entropy_SAC_df,
                         entropy_prevalence = slices_entropy_prevalence_df)

metric_df_lists2D <- list(mixed_ellipsoid = metric_df_list2D,
                          mixed_network = metric_df_list2D,
                          ringed_ellipsoid = metric_df_list2D,
                          ringed_network = metric_df_list2D,
                          separated_ellipsoid = metric_df_list2D,
                          separated_network = metric_df_list2D)



### spe 3D and 2D analysis ------------------------------------------
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")

for (arrangement in arrangements) {
  
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    for (i in seq_len(n_spes)) {
      print(i)
      ### 3D analysis -----------------------------
      spe <- simulate_spe_metadata3D(spes_metadata[[spes_metadata_index]][[i]], plot_image = F)
      spe_name <- paste("spe_", i, sep = "")  
      
      
      # pairwise_distance_data <- calculate_pairwise_distances_between_cell_types3D(spe,
      #                                                                             cell_types,
      #                                                                             show_summary = F,
      #                                                                             plot_image = F)
      # pairwise_distance_data_summary <- summarise_distances_between_cell_types3D(pairwise_distance_data)
      # 
      # ## Fill in 3 rows at a time for APD df (as we have A/A, A/B, B/B)
      # index <- 3 * (i - 1) + 1 # index is 1, 4, 7, 10 ...
      # metric_df_lists3D[[spes_metadata_index]][["APD"]][index:(index + 2), "spe"] <- spe_name
      # metric_df_lists3D[[spes_metadata_index]][["APD"]][index:(index + 2), "pair"] <- pairwise_distance_data_summary$pair
      # metric_df_lists3D[[spes_metadata_index]][["APD"]][index:(index + 2), "APD"] <- pairwise_distance_data_summary$mean
      
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
          ## Calculate ACIN and CKR as target cell type can also be the reference cell type
          
          # ACIN
          metric_df_lists3D[[spes_metadata_index]][["ACIN"]][index2, c("spe", "reference", "target")] <- c(spe_name, reference_cell_type, target_cell_type)
          metric_df_lists3D[[spes_metadata_index]][["ACIN"]][index2, radii_colnames] <- gradient_data[["cells_in_neighbourhood"]][[target_cell_type]]
          
          
          # CKR
          metric_df_lists3D[[spes_metadata_index]][["CKR"]][index2, c("spe", "reference", "target")] <- c(spe_name, reference_cell_type, target_cell_type)
          metric_df_lists3D[[spes_metadata_index]][["CKR"]][index2, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]]$cross_K_ratio
          
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
        
        proportion_SAC <- calculate_spatial_autocorrelation3D(proportion_grid_metrics, 
                                                              "proportion",
                                                              weight_method = 0.1)
        
        proportion_prevalence_df <- calculate_prevalence_gradient3D(proportion_grid_metrics,
                                                                    "proportion",
                                                                    show_AUC = F,
                                                                    plot_image = F)
        
        index <- nrow(prop_cell_types) * (i - 1) + j
        metric_df_lists3D[[spes_metadata_index]][["prop_SAC"]][index, c("spe", "reference", "target")] <- c(spe_name, prop_cell_types$ref[j], prop_cell_types$tar[j])
        metric_df_lists3D[[spes_metadata_index]][["prop_SAC"]][index, "prop_SAC"] <- proportion_SAC
        
        metric_df_lists3D[[spes_metadata_index]][["prop_prevalence"]][index, c("spe", "reference", "target")] <- c(spe_name, prop_cell_types$ref[j], prop_cell_types$tar[j])
        metric_df_lists3D[[spes_metadata_index]][["prop_prevalence"]][index, thresholds_colnames] <- proportion_prevalence_df$prevalence
      }
      
      # Get entropy grid metrics
      for (j in seq_len(nrow(entropy_cell_types))) {
        entropy_grid_metrics <- calculate_entropy_grid_metrics3D(spe, 
                                                                 n_splits,
                                                                 strsplit(entropy_cell_types$cell_types[j], ",")[[1]], 
                                                                 plot_image = F)
        
        entropy_SAC <- calculate_spatial_autocorrelation3D(entropy_grid_metrics, 
                                                           "entropy",
                                                           weight_method = 0.1)
        
        entropy_prevalence_df <- calculate_prevalence_gradient3D(entropy_grid_metrics,
                                                                 "entropy",
                                                                 show_AUC = F,
                                                                 plot_image = F)
        
        index <- nrow(entropy_cell_types) * (i - 1) + j
        metric_df_lists3D[[spes_metadata_index]][["entropy_SAC"]][index, c("spe", "cell_types")] <- c(spe_name, entropy_cell_types$cell_types[j])
        metric_df_lists3D[[spes_metadata_index]][["entropy_SAC"]][index, "entropy_SAC"] <- entropy_SAC
        
        metric_df_lists3D[[spes_metadata_index]][["entropy_prevalence"]][index, c("spe", "cell_types")] <- c(spe_name, entropy_cell_types$cell_types[j])
        metric_df_lists3D[[spes_metadata_index]][["entropy_prevalence"]][index, thresholds_colnames] <- entropy_prevalence_df$prevalence
      }  
      
      
      ### 2D slicing analysis -------------------------
      spe_slices <- get_spe_slices_list(spe, slice_bottom_z_coords, slice_top_z_coords)
      for (slice_index in seq(n_slices)) {
        spe_slice <- spe_slices[[slice_index]]
        
        # pairwise_distance_data <- calculate_pairwise_distances_between_cell_types2D(spe_slice,
        #                                                                             cell_types,
        #                                                                             show_summary = F,
        #                                                                             plot_image = F)
        # pairwise_distance_data_summary <- summarise_distances_between_cell_types2D(pairwise_distance_data)
        # 
        # ## Fill in 3 rows at a time for APD df (as we have A/A, A/B, B/B)
        # index <- n_slices * 3 * (i - 1) + 3 * (slice_index - 1) + 1
        #  metric_df_lists2D[[spes_metadata_index]][["APD"]][index:(index + 2), "spe"] <- spe_name
        # metric_df_lists2D[[spes_metadata_index]][["APD"]][index:(index + 2), "slice"] <- slice_index
        # metric_df_lists2D[[spes_metadata_index]][["APD"]][index:(index + 2), "pair"] <- pairwise_distance_data_summary$pair
        # metric_df_lists2D[[spes_metadata_index]][["APD"]][index:(index + 2), "APD"] <- pairwise_distance_data_summary$mean
        
        
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
            ## Calculate ACIN and CKR as target cell type can also be the reference cell type
            
            # ACIN & CKR
            metric_df_lists2D[[spes_metadata_index]][["ACIN"]][index2, c("spe", "slice", "reference", "target")] <- c(spe_name, slice_index, reference_cell_type, target_cell_type)
            metric_df_lists2D[[spes_metadata_index]][["CKR"]][index2, c("spe", "slice", "reference", "target")] <- c(spe_name, slice_index, reference_cell_type, target_cell_type)
            
            if (!is.null(gradient_data)) {
              metric_df_lists2D[[spes_metadata_index]][["ACIN"]][index2, radii_colnames] <- gradient_data[["cells_in_neighbourhood"]][[target_cell_type]]
              metric_df_lists2D[[spes_metadata_index]][["CKR"]][index2, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]]$cross_K_ratio
            }
            else {
              metric_df_lists2D[[spes_metadata_index]][["ACIN"]][index2, radii_colnames] <- NA
              metric_df_lists2D[[spes_metadata_index]][["CKR"]][index2, radii_colnames] <- NA
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
            proportion_SAC <- calculate_spatial_autocorrelation2D(proportion_grid_metrics, 
                                                                  "proportion",
                                                                  weight_method = 0.1)
            
            proportion_prevalence_df <- calculate_prevalence_gradient2D(proportion_grid_metrics,
                                                                        "proportion",
                                                                        show_AUC = F,
                                                                        plot_image = F)
          }
          else {
            proportion_SAC <- NA
            proportion_prevalence_df <- data.frame(threshold = seq(0.01, 1, 0.01), prevalence = NA)
          }
          
          
          
          index <- n_slices * nrow(prop_cell_types) * (i - 1) + nrow(prop_cell_types) * (slice_index - 1) + j
          
          metric_df_lists2D[[spes_metadata_index]][["prop_SAC"]][index, c("spe", "slice", "reference", "target")] <- c(spe_name, slice_index, prop_cell_types$ref[j], prop_cell_types$tar[j])
          metric_df_lists2D[[spes_metadata_index]][["prop_SAC"]][index, "prop_SAC"] <- proportion_SAC
          
          metric_df_lists2D[[spes_metadata_index]][["prop_prevalence"]][index, c("spe", "slice", "reference", "target")] <- c(spe_name, slice_index, prop_cell_types$ref[j], prop_cell_types$tar[j])
          metric_df_lists2D[[spes_metadata_index]][["prop_prevalence"]][index, thresholds_colnames] <- proportion_prevalence_df$prevalence
        }
        
        # Get entropy grid metrics
        for (j in seq_len(nrow(entropy_cell_types))) {
          entropy_grid_metrics <- calculate_entropy_grid_metrics2D(spe_slice, 
                                                                   n_splits,
                                                                   strsplit(entropy_cell_types$cell_types[j], ",")[[1]], 
                                                                   plot_image = F)
          
          if (!is.null(entropy_grid_metrics)) {
            entropy_SAC <- calculate_spatial_autocorrelation2D(entropy_grid_metrics, 
                                                               "entropy",
                                                               weight_method = 0.1)
            
            entropy_prevalence_df <- calculate_prevalence_gradient2D(entropy_grid_metrics,
                                                                     "entropy",
                                                                     show_AUC = F,
                                                                     plot_image = F)
          }
          else {
            entropy_SAC <- NA
            entropy_prevalence_df <- data.frame(threshold = seq(0.01, 1, 0.01), prevalence = NA)
          }
          
          index <- n_slices * nrow(entropy_cell_types) * (i - 1) + nrow(entropy_cell_types) * (slice_index - 1) + j
          
          metric_df_lists2D[[spes_metadata_index]][["entropy_SAC"]][index, c("spe", "slice", "cell_types")] <- c(spe_name, slice_index, entropy_cell_types$cell_types[j])
          metric_df_lists2D[[spes_metadata_index]][["entropy_SAC"]][index, "entropy_SAC"] <- entropy_SAC
          
          metric_df_lists2D[[spes_metadata_index]][["entropy_prevalence"]][index, c("spe", "slice", "cell_types")] <- c(spe_name, slice_index, entropy_cell_types$cell_types[j])
          metric_df_lists2D[[spes_metadata_index]][["entropy_prevalence"]][index, thresholds_colnames] <- entropy_prevalence_df$prevalence
        }  
      }
    }
  }
}

setwd("~/R/spaSim-3D/scripts/simulations and analysis S2/S2 data")
saveRDS(metric_df_lists3D, "metric_df_lists3D.RDS")
saveRDS(metric_df_lists2D, "metric_df_lists2D.RDS")



