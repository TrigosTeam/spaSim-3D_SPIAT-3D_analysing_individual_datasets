# Read data -----
setwd("~/R/data3D/CyCIF_colorectal_cancer")
data3D <- read.csv("colorectal_cancer_df.csv")
data3D <- data3D[, -1]
colnames(data3D)[1:3] <- c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")

# Set up data frames to contain results -----
n_slices <- length(unique(data3D$Cell.Z.Position))
cell_types <- c("Tumour", "Immune")

# Define AMD data frames as well as constants
AMD_pairs <- c("Tumour/Tumour", "Tumour/Immune", "Immune/Tumour", "Immune/Immune")

AMD_df_colnames <- c("slice", "reference", "target", "AMD")
AMD_df <- data.frame(matrix(nrow = (n_slices + 1) * length(AMD_pairs), ncol = length(AMD_df_colnames)))
colnames(AMD_df) <- AMD_df_colnames


# Define MS, NMS, ACIN, ACINP, CKR, CLR, CGR, COO, AE data frames as well as constants
radii <- seq(20, 100, 10)
radii_colnames <- paste("r", radii, sep = "")

MS_df_colnames <- c("slice", "reference", "target", radii_colnames)
MS_df <- data.frame(matrix(nrow = (n_slices + 1) * length(cell_types), ncol = length(MS_df_colnames)))
colnames(MS_df) <- MS_df_colnames

# NMS has same data frame output as MS
NMS_df <- MS_df

# Only choose prop(A) as prop(A) = 1 - prop(B) always
ACINP_df_colnames <- c("slice", "reference", "target", radii_colnames)
ACINP_df <- data.frame(matrix(nrow = (n_slices + 1) * length(cell_types), ncol = length(ACINP_df_colnames)))
colnames(ACINP_df) <- ACINP_df_colnames

# AE has same data frame output as ACINP
AE_df <- ACINP_df

## ACIN and CKR are twice as large
# (ref A and tar A or B) OR (ref B and tar B or A)
ACIN_df_colnames <- c("slice", "reference", "target", radii_colnames)
ACIN_df <- data.frame(matrix(nrow = (n_slices + 1) * length(cell_types)^2, ncol = length(ACIN_df_colnames)))
colnames(ACIN_df) <- ACIN_df_colnames

# CKR, CLR, COO, CGR have same data frame ouptut as ACIN
CKR_df <- ACIN_df
CLR_df <- ACIN_df
COO_df <- ACIN_df
CGR_df <- ACIN_df

# Define SAC and prevalence data frames as well as constants
n_splits <- 10
thresholds <- seq(0.01, 1, 0.01)
thresholds_colnames <- paste("t", thresholds, sep = "")

prop_cell_types <- data.frame(ref = c("Tumour"), tar = c("Immune"))

PBSAC_df_colnames <- c("slice", "reference", "target", "PBSAC")
PBSAC_df <- data.frame(matrix(nrow = (n_slices + 1) * nrow(prop_cell_types), ncol = length(PBSAC_df_colnames)))
colnames(PBSAC_df) <- PBSAC_df_colnames

PBP_df_colnames <- c("slice", "reference", "target", thresholds_colnames)
PBP_df <- data.frame(matrix(nrow = (n_slices + 1) * nrow(prop_cell_types), ncol = length(PBP_df_colnames)))
colnames(PBP_df) <- PBP_df_colnames


entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))

EBSAC_df_colnames <- c("slice", "cell_types", "EBSAC")
EBSAC_df <- data.frame(matrix(nrow = (n_slices + 1) * nrow(entropy_cell_types), ncol = length(EBSAC_df_colnames)))
colnames(EBSAC_df) <- EBSAC_df_colnames

EBP_df_colnames <- c("slice", "cell_types", thresholds_colnames)
EBP_df <- data.frame(matrix(nrow = (n_slices + 1) * nrow(entropy_cell_types), ncol = length(EBP_df_colnames)))
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


# Analysis -----
n_slices <- length(unique(data3D$Cell.Z.Position))
slice_z_coords <- unique(data3D$Cell.Z.Position)
colnames(data3D)[4:5] <- c("Cell.Type.Specific", "Cell.Type")
data3D[data3D$Cell.Type == "Tumor/Epi", "Cell.Type"] <- "Tumour"


for (i in seq(n_slices + 1)) {
  print(i)
  
  # i represents the current slice index
  if (i == n_slices + 1) {
    df <- data3D[data3D$Cell.Z.Position == slice_z_coords[i], ]
  }
  # if i == 0, analyse in 3D instead
  else {
    df <- data3D
  }
  
  ### 3D analysis -----------------------------
  if (i == n_slices + 1) {
    
    minimum_distance_data <- calculate_minimum_distances_between_cell_types3D(df,
                                                                              cell_types,
                                                                              show_summary = F,
                                                                              plot_image = F)
    
    minimum_distance_data_summary <- summarise_distances_between_cell_types3D(minimum_distance_data)
    ## Fill in 4 rows at a time for AMD df (as we have A/A, A/B, B/A, B/B)
    index <- 4 * (i - 1) + 1 # index is 1, 5, 9, 13
    
    metric_df_list[["AMD"]][index:(index + 3), "slice"] <- i
    metric_df_list[["AMD"]][index:(index + 3), "reference"] <- minimum_distance_data_summary$reference
    metric_df_list[["AMD"]][index:(index + 3), "target"] <- minimum_distance_data_summary$target
    metric_df_list[["AMD"]][index:(index + 3), "AMD"] <- minimum_distance_data_summary$mean
    
    
    
    
    
    index1 <- 2 * (i - 1) + 1 # index1 is 1, 3, 5, ...
    index2 <- 4 * (i - 1) + 1 # index2 is 1, 5, 9, 13...
    for (reference_cell_type in cell_types) {
      gradient_data <- calculate_all_gradient_cc_metrics3D(df,
                                                           reference_cell_type,
                                                           cell_types,
                                                           radii,
                                                           plot_image = F)
      
      target_cell_type <- setdiff(cell_types, reference_cell_type)
      
      metric_df_list[["MS"]][index1, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
      metric_df_list[["MS"]][index1, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$mixing_score
      
      metric_df_list[["NMS"]][index1, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
      metric_df_list[["NMS"]][index1, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$normalised_mixing_score
      
      metric_df_list[["ACINP"]][index1, c("slice", "reference", "target")] <- c(i, reference_cell_type, "B")
      metric_df_list[["ACINP"]][index1, radii_colnames] <- gradient_data[["cells_in_neighbourhood_proportion"]][["B"]]
      
      metric_df_list[["AE"]][index1, c("slice", "reference", "target")] <- c(i, reference_cell_type, "A,B")
      metric_df_list[["AE"]][index1, radii_colnames] <- gradient_data[["entropy"]]$entropy
      
      index1 <- index1 + 1
      
      for (target_cell_type in cell_types) {
        ## Calculate ACIN, CKR, CKL, COO as target cell type can also be the reference cell type
        
        # ACIN
        metric_df_list[["ACIN"]][index2, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["ACIN"]][index2, radii_colnames] <- gradient_data[["cells_in_neighbourhood"]][[target_cell_type]]
        
        # CKR
        metric_df_list[["CKR"]][index2, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["CKR"]][index2, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]] / gradient_data[["cross_K"]][["expected"]]
        
        # CLR
        metric_df_list[["CLR"]][index2, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["CLR"]][index2, radii_colnames] <- gradient_data[["cross_L"]][[target_cell_type]] / gradient_data[["cross_L"]][["expected"]]
        
        # COO
        metric_df_list[["COO"]][index2, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["COO"]][index2, radii_colnames] <- gradient_data[["co_occurrence"]][[target_cell_type]]
        
        # CGR
        metric_df_list[["CGR"]][index2, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["CGR"]][index2, radii_colnames] <- gradient_data[["cross_G"]][[target_cell_type]][["observed_cross_G"]] / gradient_data[["cross_G"]][[target_cell_type]][["expected_cross_G"]]
        
        index2 <- index2 + 1
      }
    }
    
    
    # Get proportion grid metrics
    for (j in seq_len(nrow(prop_cell_types))) {
      proportion_grid_metrics <- calculate_cell_proportion_grid_metrics3D(df, 
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
      metric_df_list[["PBSAC"]][index, c("slice", "reference", "target")] <- c(i, prop_cell_types$ref[j], prop_cell_types$tar[j])
      metric_df_list[["PBSAC"]][index, "PBSAC"] <- PBSAC
      
      metric_df_list[["PBP"]][index, c("slice", "reference", "target")] <- c(i, prop_cell_types$ref[j], prop_cell_types$tar[j])
      metric_df_list[["PBP"]][index, thresholds_colnames] <- PBP_df$prevalence
    }
    
    # Get entropy grid metrics
    for (j in seq_len(nrow(entropy_cell_types))) {
      entropy_grid_metrics <- calculate_entropy_grid_metrics3D(df, 
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
      metric_df_list[["EBSAC"]][index, c("slice", "cell_types")] <- c(i, entropy_cell_types$cell_types[j])
      metric_df_list[["EBSAC"]][index, "EBSAC"] <- EBSAC
      
      metric_df_list[["EBP"]][index, c("slice", "cell_types")] <- c(i, entropy_cell_types$cell_types[j])
      metric_df_list[["EBP"]][index, thresholds_colnames] <- EBP_df$prevalence
    }  
  }
  
  
  ### 2D analysis -----------------------------
  else {
    
    minimum_distance_data <- calculate_minimum_distances_between_cell_types2D(df,
                                                                              cell_types,
                                                                              show_summary = F,
                                                                              plot_image = F)
    
    minimum_distance_data_summary <- summarise_distances_between_cell_types2D(minimum_distance_data)
    ## Fill in 4 rows at a time for AMD df (as we have A/A, A/B, B/A, B/B)
    index <- 4 * (i - 1) + 1 # index is 1, 5, 9, 13
    
    metric_df_list[["AMD"]][index:(index + 3), "slice"] <- i
    metric_df_list[["AMD"]][index:(index + 3), "reference"] <- minimum_distance_data_summary$reference
    metric_df_list[["AMD"]][index:(index + 3), "target"] <- minimum_distance_data_summary$target
    metric_df_list[["AMD"]][index:(index + 3), "AMD"] <- minimum_distance_data_summary$mean
    
    
    
    
    
    index1 <- 2 * (i - 1) + 1 # index1 is 1, 3, 5, ...
    index2 <- 4 * (i - 1) + 1 # index2 is 1, 5, 9, 13...
    for (reference_cell_type in cell_types) {
      gradient_data <- calculate_all_gradient_cc_metrics2D(df,
                                                           reference_cell_type,
                                                           cell_types,
                                                           radii,
                                                           plot_image = F)
      
      target_cell_type <- setdiff(cell_types, reference_cell_type)
      
      metric_df_list[["MS"]][index1, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
      metric_df_list[["MS"]][index1, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$mixing_score
      
      metric_df_list[["NMS"]][index1, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
      metric_df_list[["NMS"]][index1, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$normalised_mixing_score
      
      metric_df_list[["ACINP"]][index1, c("slice", "reference", "target")] <- c(i, reference_cell_type, "B")
      metric_df_list[["ACINP"]][index1, radii_colnames] <- gradient_data[["cells_in_neighbourhood_proportion"]][["B"]]
      
      metric_df_list[["AE"]][index1, c("slice", "reference", "target")] <- c(i, reference_cell_type, "A,B")
      metric_df_list[["AE"]][index1, radii_colnames] <- gradient_data[["entropy"]]$entropy
      
      index1 <- index1 + 1
      
      for (target_cell_type in cell_types) {
        ## Calculate ACIN, CKR, CKL, COO as target cell type can also be the reference cell type
        
        # ACIN
        metric_df_list[["ACIN"]][index2, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["ACIN"]][index2, radii_colnames] <- gradient_data[["cells_in_neighbourhood"]][[target_cell_type]]
        
        # CKR
        metric_df_list[["CKR"]][index2, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["CKR"]][index2, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]] / gradient_data[["cross_K"]][["expected"]]
        
        # CLR
        metric_df_list[["CLR"]][index2, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["CLR"]][index2, radii_colnames] <- gradient_data[["cross_L"]][[target_cell_type]] / gradient_data[["cross_L"]][["expected"]]
        
        # COO
        metric_df_list[["COO"]][index2, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["COO"]][index2, radii_colnames] <- gradient_data[["co_occurrence"]][[target_cell_type]]
        
        # CGR
        metric_df_list[["CGR"]][index2, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
        metric_df_list[["CGR"]][index2, radii_colnames] <- gradient_data[["cross_G"]][[target_cell_type]][["observed_cross_G"]] / gradient_data[["cross_G"]][[target_cell_type]][["expected_cross_G"]]
        
        index2 <- index2 + 1
      }
    }
    
    
    # Get proportion grid metrics
    for (j in seq_len(nrow(prop_cell_types))) {
      proportion_grid_metrics <- calculate_cell_proportion_grid_metrics2D(df, 
                                                                          n_splits,
                                                                          strsplit(prop_cell_types$ref[j], ",")[[1]], 
                                                                          strsplit(prop_cell_types$tar[j], ",")[[1]],
                                                                          plot_image = F)
      
      PBSAC <- calculate_spatial_autocorrelation2D(proportion_grid_metrics, 
                                                   "proportion",
                                                   weight_method = 0.1)
      
      PBP_df <- calculate_prevalence_gradient2D(proportion_grid_metrics,
                                                "proportion",
                                                show_AUC = F,
                                                plot_image = F)
      
      index <- nrow(prop_cell_types) * (i - 1) + j
      metric_df_list[["PBSAC"]][index, c("slice", "reference", "target")] <- c(i, prop_cell_types$ref[j], prop_cell_types$tar[j])
      metric_df_list[["PBSAC"]][index, "PBSAC"] <- PBSAC
      
      metric_df_list[["PBP"]][index, c("slice", "reference", "target")] <- c(i, prop_cell_types$ref[j], prop_cell_types$tar[j])
      metric_df_list[["PBP"]][index, thresholds_colnames] <- PBP_df$prevalence
    }
    
    # Get entropy grid metrics
    for (j in seq_len(nrow(entropy_cell_types))) {
      entropy_grid_metrics <- calculate_entropy_grid_metrics2D(df, 
                                                               n_splits,
                                                               strsplit(entropy_cell_types$cell_types[j], ",")[[1]], 
                                                               plot_image = F)
      
      EBSAC <- calculate_spatial_autocorrelation2D(entropy_grid_metrics, 
                                                   "entropy",
                                                   weight_method = 0.1)
      
      EBP_df <- calculate_prevalence_gradient2D(entropy_grid_metrics,
                                                "entropy",
                                                show_AUC = F,
                                                plot_image = F)
      
      index <- nrow(entropy_cell_types) * (i - 1) + j
      metric_df_list[["EBSAC"]][index, c("slice", "cell_types")] <- c(i, entropy_cell_types$cell_types[j])
      metric_df_list[["EBSAC"]][index, "EBSAC"] <- EBSAC
      
      metric_df_list[["EBP"]][index, c("slice", "cell_types")] <- c(i, entropy_cell_types$cell_types[j])
      metric_df_list[["EBP"]][index, thresholds_colnames] <- EBP_df$prevalence
    }  
  }
}
