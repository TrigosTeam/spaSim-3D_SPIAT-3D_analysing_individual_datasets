### Read data -----
setwd("~/R/data3D/openST_human_metastatic_lymph_node")
data3D <- read.csv("human_metastatic_lymph_node_df.csv")

### Set up data frames to contain results -----
n_slices <- length(unique(data3D$Cell.Z.Position))
cell_types <- unique(data3D$Cell.Type)
cell_types <- cell_types[cell_types != "unknown"]
n_cell_type_combinations <- length(cell_types)^2

# Define AMD data frames as well as constants

AMD_df_colnames <- c("slice", "reference", "target", "AMD")
AMD_df <- data.frame(matrix(nrow = (n_slices + 1) * n_cell_type_combinations, ncol = length(AMD_df_colnames)))
colnames(AMD_df) <- AMD_df_colnames


# Define MS, NMS, ACIN, ACINP, CKR, CLR, CGR, COO, AE data frames as well as constants
radii <- seq(20, 100, 10)
radii_colnames <- paste("r", radii, sep = "")

MS_df_colnames <- c("slice", "reference", "target", radii_colnames)
MS_df <- data.frame(matrix(nrow = (n_slices + 1) * n_cell_type_combinations, ncol = length(MS_df_colnames)))
colnames(MS_df) <- MS_df_colnames

# NMS has same data frame output as MS
NMS_df <- MS_df

# Only choose prop(A) as prop(A) = 1 - prop(B) always
ACINP_df_colnames <- c("slice", "reference", "target", radii_colnames)
ACINP_df <- data.frame(matrix(nrow = (n_slices + 1) * n_cell_type_combinations, ncol = length(ACINP_df_colnames)))
colnames(ACINP_df) <- ACINP_df_colnames

# AE has same data frame output as ACINP
AE_df <- ACINP_df

## ACIN and CKR are twice as large
# (ref A and tar A or B) OR (ref B and tar B or A)
ACIN_df_colnames <- c("slice", "reference", "target", radii_colnames)
ACIN_df <- data.frame(matrix(nrow = (n_slices + 1) * n_cell_type_combinations, ncol = length(ACIN_df_colnames)))
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


### Analysis -----
n_slices <- length(unique(data3D$Cell.Z.Position))
slice_z_coords <- unique(data3D$Cell.Z.Position)
n_cell_type_combinations <- length(cell_types)^2

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
        
        if (is.null(gradient_data)) {
          metric_df_list[["ACIN"]][gradient_index, radii_colnames] <- NA
          metric_df_list[["ACINP"]][gradient_index, radii_colnames] <- NA
          metric_df_list[["CKR"]][gradient_index, radii_colnames] <- NA
          metric_df_list[["CLR"]][gradient_index, radii_colnames] <- NA
          metric_df_list[["COO"]][gradient_index, radii_colnames] <- NA
          metric_df_list[["CGR"]][gradient_index, radii_colnames] <- NA
          metric_df_list[["MS"]][gradient_index, radii_colnames] <- NA
          metric_df_list[["NMS"]][gradient_index, radii_colnames] <- NA
          metric_df_list[["AE"]][gradient_index, radii_colnames] <- NA
        }
        else {
          metric_df_list[["ACIN"]][gradient_index, radii_colnames] <- gradient_data[["cells_in_neighbourhood"]][[target_cell_type]]
          metric_df_list[["ACINP"]][gradient_index, radii_colnames] <- gradient_data[["cells_in_neighbourhood_proportion"]][[target_cell_type]]
          metric_df_list[["CKR"]][gradient_index, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]] / gradient_data[["cross_K"]][["expected"]]
          metric_df_list[["CLR"]][gradient_index, radii_colnames] <- gradient_data[["cross_L"]][[target_cell_type]] / gradient_data[["cross_L"]][["expected"]]
          metric_df_list[["COO"]][gradient_index, radii_colnames] <- gradient_data[["co_occurrence"]][[target_cell_type]]
          metric_df_list[["CGR"]][gradient_index, radii_colnames] <- gradient_data[["cross_G"]][[target_cell_type]][["observed_cross_G"]] / gradient_data[["cross_G"]][[target_cell_type]][["expected_cross_G"]]
          
          if (reference_cell_type != target_cell_type) {
            metric_df_list[["MS"]][gradient_index, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$mixing_score
            metric_df_list[["NMS"]][gradient_index, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$normalised_mixing_score
            metric_df_list[["AE"]][gradient_index, radii_colnames] <- gradient_data[["entropy"]]$entropy
          }
          else {
            metric_df_list[["MS"]][gradient_index, radii_colnames] <- Inf
            metric_df_list[["NMS"]][gradient_index, radii_colnames] <- Inf
            metric_df_list[["AE"]][gradient_index, radii_colnames] <- Inf
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
        
        if (is.null(gradient_data)) {
          metric_df_list[["ACIN"]][gradient_index, radii_colnames] <- NA
          metric_df_list[["ACINP"]][gradient_index, radii_colnames] <- NA
          metric_df_list[["CKR"]][gradient_index, radii_colnames] <- NA
          metric_df_list[["CLR"]][gradient_index, radii_colnames] <- NA
          metric_df_list[["COO"]][gradient_index, radii_colnames] <- NA
          metric_df_list[["CGR"]][gradient_index, radii_colnames] <- NA
          metric_df_list[["MS"]][gradient_index, radii_colnames] <- NA
          metric_df_list[["NMS"]][gradient_index, radii_colnames] <- NA
          metric_df_list[["AE"]][gradient_index, radii_colnames] <- NA
        }
        else {
          metric_df_list[["ACIN"]][gradient_index, radii_colnames] <- gradient_data[["cells_in_neighbourhood"]][[target_cell_type]]
          metric_df_list[["ACINP"]][gradient_index, radii_colnames] <- gradient_data[["cells_in_neighbourhood_proportion"]][[target_cell_type]]
          metric_df_list[["CKR"]][gradient_index, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]] / gradient_data[["cross_K"]][["expected"]]
          metric_df_list[["CLR"]][gradient_index, radii_colnames] <- gradient_data[["cross_L"]][[target_cell_type]] / gradient_data[["cross_L"]][["expected"]]
          metric_df_list[["COO"]][gradient_index, radii_colnames] <- gradient_data[["co_occurrence"]][[target_cell_type]]
          metric_df_list[["CGR"]][gradient_index, radii_colnames] <- gradient_data[["cross_G"]][[target_cell_type]][["observed_cross_G"]] / gradient_data[["cross_G"]][[target_cell_type]][["expected_cross_G"]]
          
          if (reference_cell_type != target_cell_type) {
            metric_df_list[["MS"]][gradient_index, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$mixing_score
            metric_df_list[["NMS"]][gradient_index, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$normalised_mixing_score
            metric_df_list[["AE"]][gradient_index, radii_colnames] <- gradient_data[["entropy"]]$entropy
          }
          else {
            metric_df_list[["MS"]][gradient_index, radii_colnames] <- Inf
            metric_df_list[["NMS"]][gradient_index, radii_colnames] <- Inf
            metric_df_list[["AE"]][gradient_index, radii_colnames] <- Inf
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

setwd("~/R/SPIAT-3D_benchmarking/public_3D_data_analysis")
saveRDS(metric_df_list, "openST_human_metastatic_lymph_node_metric_df_list.RDS")

### Plot analysis of 2D and 3D data -----
setwd("~/R/SPIAT-3D_benchmarking/public_3D_data_analysis")
metric_df_list <- readRDS("openST_human_metastatic_lymph_node_metric_df_list.RDS")

get_gradient <- function(metric) {
  if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "CLR", "COO", "CGR")) {
    return("radius")
  }
  else if (metric %in% c("PBP", "EBP")) {
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

gradient_radii_metrics <- c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "CLR", "COO", "CGR")


for (metric in gradient_radii_metrics) {
  metric_AUC_name <- paste(metric, "AUC", sep = "_")
  
  if (metric %in% c("MS", "NMS", "ACIN", "CKR", "CLR", "COO", "CGR")) {
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

# PBP_AUC 3D
PBP_df <- metric_df_list[["PBP"]]
PBP_df$PBP_AUC <- apply(PBP_df[ , threshold_colnames], 1, sum) * 0.01
PBP_AUC_df <- PBP_df[ , c("slice", "reference", "target", "PBP_AUC")]
metric_df_list[["PBP_AUC"]] <- PBP_AUC_df

# EBP_AUC 3D
EBP_df <- metric_df_list[["EBP"]]
EBP_df$EBP_AUC <- apply(EBP_df[ , threshold_colnames], 1, sum) * 0.01
EBP_AUC_df <- EBP_df[ , c("slice", "cell_types", "EBP_AUC")]
metric_df_list[["EBP_AUC"]] <- EBP_AUC_df


## Functions to plot
plot_3D_vs_2D <- function(metric_df_list,
                          metric) {
  
  # Get metric_df for current metric
  metric_df <- metric_df_list[[metric]]
  
  metric_df <- metric_df[metric_df$reference %in% c("CAF", "CAM", "Tumor"), ]
  metric_df <- metric_df[metric_df$target %in% c("CAF", "CAM", "Tumor"), ]
  
  # Change and further subset columns of metric_df_subset
  colnames(metric_df)[colnames(metric_df) == metric] <- "value"

  # Add dummy column
  metric_df$dummy <- paste("dummy")
  
  fig <- ggplot(metric_df, aes(x = dummy, y = metric)) +
    geom_boxplot(data = metric_df[metric_df$slice != max(metric_df$slice), ],
                 outlier.shape = NA, fill = "lightgray") +
    geom_jitter(data = metric_df[metric_df$slice != max(metric_df$slice), ],
                width = 0.2, alpha = 0.5, color = "#0062c5") +
    geom_point(data = metric_df[metric_df$slice == max(metric_df$slice), ],
               shape = 8, color = "#bb0036", size = 3) +  # Red stars for 3D value
    labs(title = "", x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.line.x = element_blank()
    ) +
    facet_grid(reference ~ target, scales = "free_y")
    
  
  return(fig)
}

plot_3D_vs_error_box_plot <- function(metric_df_list,
                                      metrics) {
  
  # Put all data into this data frame
  combined_plot_df <- data.frame()
  
  
  for (metric in metrics) {
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
    
    # Calculate error for each slice, and remove 3D row
    value_3D <- metric_df_subset[["value"]][metric_df_subset[["slice"]] == 0]
    metric_df_subset[["error"]] <- ((metric_df_subset[["value"]] - value_3D) / value_3D) * 100
    metric_df_subset["value"] <- NULL
    metric_df_subset <- metric_df_subset[metric_df_subset[["slice"]] != 0, ]
    
    combined_plot_df <- rbind(combined_plot_df, metric_df_subset)
  }
  
  fig <- ggplot(combined_plot_df, aes(x = metric, y = error)) +
    geom_boxplot(outlier.shape = NA, fill = "lightgray") +  # Hide default outliers to avoid duplication
    geom_jitter(width = 0.2, alpha = 0.5, color = "#0062c5") +  # Add dots with slight horizontal jitter
    labs(title = "Error Distribution by Metric",
         x = "Metric",
         y = "Error (%)") +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
  
  return(fig)
}

## Get plot
metrics <- c("AMD", "ACIN_AUC", "ACINP_AUC", "AE_AUC", "MS_AUC", "NMS_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")


fig_3D_vs_2D <- plot_3D_vs_2D(metric_df_list,
                              metric)

fig_3D_vs_error_box_plot <- plot_3D_vs_error_box_plot(metric_df_list,
                                                      metrics)
methods::show(fig_3D_vs_2D)
methods::show(fig_3D_vs_error_box_plot)


setwd("~/R/plots/public_data")
pdf("openST_human_metastatic_lymph_node.pdf", width = 10, height = 8)

print(fig_3D_vs_2D)
print(fig_3D_vs_error_box_plot)

dev.off()
