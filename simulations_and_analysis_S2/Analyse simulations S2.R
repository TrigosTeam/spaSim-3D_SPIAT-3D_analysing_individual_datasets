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

setwd("~/R/SPIAT-3D_benchmarking/simulations_and_analysis_S2/S2_data")
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


# Define MS, NMS, ACIN, ACINP, CKR, CLR, CGR, COO, AE data frames as well as constants
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
                         CGR = CGR_df,
                         COO = COO_df,
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
          metric_df_lists3D[[spes_metadata_index]][["CKR"]][index2, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]]
          
          # CLR
          metric_df_lists3D[[spes_metadata_index]][["CLR"]][index2, c("spe", "reference", "target")] <- c(spe_name, reference_cell_type, target_cell_type)
          metric_df_lists3D[[spes_metadata_index]][["CLR"]][index2, radii_colnames] <- gradient_data[["cross_L"]][[target_cell_type]]
          
          # COO
          metric_df_lists3D[[spes_metadata_index]][["COO"]][index2, c("spe", "reference", "target")] <- c(spe_name, reference_cell_type, target_cell_type)
          metric_df_lists3D[[spes_metadata_index]][["COO"]][index2, radii_colnames] <- gradient_data[["co_occurrence"]][[target_cell_type]]
          
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
        metric_df_lists3D[[spes_metadata_index]][["PBSAC"]][index, c("spe", "reference", "target")] <- c(spe_name, prop_cell_types$ref[j], prop_cell_types$tar[j])
        metric_df_lists3D[[spes_metadata_index]][["PBSAC"]][index, "PBSAC"] <- proportion_SAC
        
        metric_df_lists3D[[spes_metadata_index]][["PBP"]][index, c("spe", "reference", "target")] <- c(spe_name, prop_cell_types$ref[j], prop_cell_types$tar[j])
        metric_df_lists3D[[spes_metadata_index]][["PBP"]][index, thresholds_colnames] <- proportion_prevalence_df$prevalence
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
            
            # ACIN & CKR
            metric_df_lists2D[[spes_metadata_index]][["ACIN"]][index2, c("spe", "slice", "reference", "target")] <- c(spe_name, slice_index, reference_cell_type, target_cell_type)
            metric_df_lists2D[[spes_metadata_index]][["CKR"]][index2, c("spe", "slice", "reference", "target")] <- c(spe_name, slice_index, reference_cell_type, target_cell_type)
            metric_df_lists2D[[spes_metadata_index]][["CLR"]][index2, c("spe", "slice", "reference", "target")] <- c(spe_name, slice_index, reference_cell_type, target_cell_type)
            metric_df_lists2D[[spes_metadata_index]][["COO"]][index2, c("spe", "slice", "reference", "target")] <- c(spe_name, slice_index, reference_cell_type, target_cell_type)
            
            if (!is.null(gradient_data)) {
              metric_df_lists2D[[spes_metadata_index]][["ACIN"]][index2, radii_colnames] <- gradient_data[["cells_in_neighbourhood"]][[target_cell_type]]
              metric_df_lists2D[[spes_metadata_index]][["CKR"]][index2, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]]
              metric_df_lists2D[[spes_metadata_index]][["CLR"]][index2, radii_colnames] <- gradient_data[["cross_L"]][[target_cell_type]]
              metric_df_lists2D[[spes_metadata_index]][["COO"]][index2, radii_colnames] <- gradient_data[["co_occurrence"]][[target_cell_type]]
            }
            else {
              metric_df_lists2D[[spes_metadata_index]][["ACIN"]][index2, radii_colnames] <- NA
              metric_df_lists2D[[spes_metadata_index]][["CKR"]][index2, radii_colnames] <- NA
              metric_df_lists2D[[spes_metadata_index]][["CLR"]][index2, radii_colnames] <- NA
              metric_df_lists2D[[spes_metadata_index]][["COO"]][index2, radii_colnames] <- NA
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
          
          metric_df_lists2D[[spes_metadata_index]][["PBSAC"]][index, c("spe", "slice", "reference", "target")] <- c(spe_name, slice_index, prop_cell_types$ref[j], prop_cell_types$tar[j])
          metric_df_lists2D[[spes_metadata_index]][["PBSAC"]][index, "PBSAC"] <- proportion_SAC
          
          metric_df_lists2D[[spes_metadata_index]][["PBP"]][index, c("spe", "slice", "reference", "target")] <- c(spe_name, slice_index, prop_cell_types$ref[j], prop_cell_types$tar[j])
          metric_df_lists2D[[spes_metadata_index]][["PBP"]][index, thresholds_colnames] <- proportion_prevalence_df$prevalence
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

setwd("~/R/SPIAT-3D_benchmarking/simulations_and_analysis_S2/S2_data")
saveRDS(metric_df_lists3D, "metric_df_lists3D.RDS")
saveRDS(metric_df_lists2D, "metric_df_lists2D.RDS")



