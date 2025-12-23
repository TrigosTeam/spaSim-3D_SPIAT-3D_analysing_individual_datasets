### Read data and set file name (THE ONLY PART YOU NEED TO CHANGE) ------
setwd("***directory to your analysis***") # e.g. "~/R/SPIAT-3D_benchmarking/public_3D_data_analysis/metric_df_lists"
metric_df_list <- readRDS("***your metric_df_list.RDS***") # e.g. "CyCIF_colorectal_cancer_metric_df_list.RDS"

file_name <- "*** file_save_name.pdf ***" # e.g. openST_human_metastatic_lymph_node.pdf

### Libraries -----
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

### Format data -----

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

gradient_radii_metrics <- c("MS", "NMS", "ACIN", "ANE", "ANC", "CKR", "CLR", "COO", "CGR")


for (metric in gradient_radii_metrics) {
  metric_AUC_name <- paste(metric, "AUC", sep = "_")
  
  if (metric %in% c("MS", "NMS", "ANC", "CKR", "CLR", "COO", "CGR")) {
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


### Plot analysis -----

## Functions to plot
# Utility function to get metric cell types
get_metric_cell_types <- function(metric) {
  # Get metric_cell_types
  if (metric %in% c("AMD", "ANC", "ACIN", "CKR", "ANC_AUC", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC")) {
    metric_cell_types <- data.frame(ref = c("Tumour"), tar = c("Immune"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("MS", "NMS", "MS_AUC", "NMS_AUC")) {
    metric_cell_types <- data.frame(ref = c("Tumour"), tar = c("Immune"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("ANE", "ANE_AUC")) {
    metric_cell_types <- data.frame(ref = c("Tumour"), tar = c("Tumour,Immune"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("PBSAC", "PBP", "PBP_AUC")) {
    metric_cell_types <- data.frame(ref = c("Tumour"), tar = c("Immune"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("EBSAC", "EBP", "EBP_AUC")) {
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
                             index) {
  
  metric_cell_types <- get_metric_cell_types(metric)
  
  if (metric %in% c("AMD", "ANC", "ACIN", "CKR", "CLR", "COO", "CGR", "MS", "NMS", "ANC_AUC", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "MS_AUC", "NMS_AUC", "PBSAC", "PBP", "PBP_AUC")) {
    metric_df_subset <- metric_df[metric_df$reference == metric_cell_types[index, "ref"] & metric_df$target == metric_cell_types[index, "tar"], ] 
  }
  else if (metric %in% c("ANE", "ANE_AUC")) {
    metric_df_subset <- metric_df[metric_df$reference == metric_cell_types[index, "ref"] & metric_df$target == metric_cell_types[index, "tar"], ] 
  }
  else if (metric %in% c("EBSAC", "EBP", "EBP_AUC")) {
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
  
  # Change and further subset columns of metric_df_subset
  colnames(metric_df)[colnames(metric_df) == metric] <- "value"
  
  # Add dummy column
  metric_df$dummy <- paste("dummy")
  
  fig <- ggplot(metric_df, aes(x = dummy, y = value)) +
    geom_boxplot(data = metric_df[metric_df$slice != as.character(max(as.numeric(metric_df$slice))), ],
                 outlier.shape = NA, fill = "lightgray") +
    geom_jitter(data = metric_df[metric_df$slice != as.character(max(as.numeric(metric_df$slice))), ],
                width = 0.2, alpha = 0.5, color = "#0062c5") +
    geom_point(data = metric_df[metric_df$slice == as.character(max(as.numeric(metric_df$slice))), ],
               shape = 8, color = "#bb0036", size = 3) +  # Red stars for 3D value
    labs(title = metric, x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.line.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.line.y = element_blank()
    ) +
    facet_grid(reference ~ target, scales = "free_y")
  
  
  return(fig)
}

plot_3D_vs_error_by_pair_box_plot <- function(metric_df_list,
                                              metric) {
  
  # Get metric_df for current metric
  metric_df <- metric_df_list[[metric]]
  
  # Change and further subset columns of metric_df
  colnames(metric_df)[colnames(metric_df) == metric] <- "value"
  
  # Obtain 3D value
  value_3D <- metric_df[["value"]][metric_df[["slice"]] == as.character(max(as.numeric(metric_df$slice)))]
  
  # Calculate error
  metric_df[["error"]] <- ((metric_df[["value"]] - value_3D) / value_3D) * 100
  
  # Remove value column (only using error now)
  metric_df["value"] <- NULL
  
  # Remove 3D data (integrated into error)
  metric_df <- metric_df[metric_df[["slice"]] != max(as.integer(metric_df$slice)), ]
  
  # Add pair column
  metric_df$pair <- paste(metric_df$reference, metric_df$target, sep = "/")
  
  fig <- ggplot(metric_df, aes(x = pair, y = error)) +
    geom_boxplot(outlier.shape = NA, fill = "lightgray") +  # Hide default outliers to avoid duplication
    geom_jitter(width = 0.2, alpha = 0.5, color = "#0062c5") +  # Add dots with slight horizontal jitter
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    labs(title = "Error Distribution by Pair",
         x = "Pair",
         y = "Error (%)") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)  # Rotate x-axis labels vertically
    )
  
  return(fig)
}

plot_3D_vs_error_by_slice_box_plot <- function(metric_df_list,
                                               metric) {
  
  # Get metric_df for current metric
  metric_df <- metric_df_list[[metric]]
  
  # Change and further subset columns of metric_df
  colnames(metric_df)[colnames(metric_df) == metric] <- "value"
  
  # Obtain 3D value
  value_3D <- metric_df[["value"]][metric_df[["slice"]] == as.character(max(as.numeric(metric_df$slice)))]
  
  # Calculate error
  metric_df[["error"]] <- ((metric_df[["value"]] - value_3D) / value_3D) * 100
  
  # Remove value column (only using error now)
  metric_df["value"] <- NULL
  
  # Remove 3D data (integrated into error)
  metric_df <- metric_df[metric_df[["slice"]] != max(as.integer(metric_df$slice)), ]
  
  # Factor slice column so it is from 1, 2, 3, ...
  metric_df$slice <- factor(metric_df$slice, as.character(sort(as.integer(unique(metric_df$slice)))))
  
  fig <- ggplot(metric_df, aes(x = slice, y = error)) +
    geom_boxplot(outlier.shape = NA, fill = "lightgray") +  # Hide default outliers to avoid duplication
    geom_jitter(width = 0.2, alpha = 0.5, color = "#0062c5") +  # Add dots with slight horizontal jitter
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    labs(title = "Error Distribution by Slice",
         x = "Slice Index",
         y = "Error (%)") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)  # Rotate x-axis labels vertically
    )
  
  return(fig)
}


plot_3D_vs_error_all_metrics_by_pair_box_plot <- function(metric_df_list,
                                                          metrics) {
  
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change and further subset columns of metric_df
    colnames(metric_df)[colnames(metric_df) == metric] <- "value"
    
    # Obtain 3D value
    value_3D <- metric_df[["value"]][metric_df[["slice"]] == as.character(max(as.numeric(metric_df$slice)))]
    
    # Calculate error
    metric_df[["error"]] <- ((metric_df[["value"]] - value_3D) / value_3D) * 100
    
    # Remove value column (only using error now)
    metric_df["value"] <- NULL
    
    # Remove 3D data (integrated into error)
    metric_df <- metric_df[metric_df[["slice"]] != max(as.integer(metric_df$slice)), ]
    
    if (!(metric %in% c("EBSAC", "EBP_AUC"))) {
      # Add pair column
      metric_df$pair <- paste(metric_df$reference, metric_df$target, sep = "/")
      
      # Extract median values for error for each pair
      median_df <- metric_df %>%
        group_by(pair) %>%
        dplyr::summarize(median_error = median(error, na.rm = TRUE), .groups = "drop")
      
      median_df$pair <- NULL
    }
    else {
      # Extract median values for error for each pair
      median_df <- metric_df %>%
        group_by(cell_types) %>%
        dplyr::summarize(median_error = median(error, na.rm = TRUE), .groups = "drop")
      
      median_df$cell_types <- NULL
    }
    
    median_df$metric <- metric  
    
    # Add median_df to plot_df
    plot_df <- rbind(plot_df, median_df)
  }
  
  fig <- ggplot(plot_df, aes(x = metric, y = median_error)) +
    geom_boxplot(outlier.shape = NA, fill = "lightgray") +  # Hide default outliers to avoid duplication
    geom_jitter(width = 0.2, alpha = 0.5, color = "#0062c5") +  # Add dots with slight horizontal jitter
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    labs(title = "Error Distribution by Metric, showing Median Error for each Pair",
         x = "Metric",
         y = "Error (%)") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  return(fig)
}

plot_3D_vs_error_all_metrics_by_all_pairs_box_plot <- function(metric_df_list,
                                                               metrics) {
  
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change and further subset columns of metric_df
    colnames(metric_df)[colnames(metric_df) == metric] <- "value"
    
    # Obtain 3D value
    value_3D <- metric_df[["value"]][metric_df[["slice"]] == as.character(max(as.numeric(metric_df$slice)))]
    
    # Calculate error
    metric_df[["error"]] <- ((metric_df[["value"]] - value_3D) / value_3D) * 100
    
    # Remove value column (only using error now)
    metric_df["value"] <- NULL
    
    # Remove 3D data (integrated into error)
    metric_df <- metric_df[metric_df[["slice"]] != max(as.integer(metric_df$slice)), ]
    
    if (metric %in% c("EBSAC", "EBP_AUC")) {
      # For EBSAC and EBP_AUC, assume pair is the same as cell_types for consistency
      metric_df$pair <- gsub(',', '/', metric_df$cell_types)
    }
    else if (metric %in% c("AE_AUC")) {
      # For AE_AUC, assume pair is the same as target_cell_type for consistency (as target is of form A,B already)
      metric_df$pair <- gsub(',', '/', metric_df$target)
    }
    else {
      # Add reference-target column
      metric_df$pair <- paste(metric_df$reference, metric_df$target, sep = "/")
    }
    
    metric_df$metric <- metric  
    
    # Add metric_df to plot_df
    plot_df <- rbind(plot_df, metric_df[ , c("error", "slice", "pair", "metric")])
  }
  
  fig <- ggplot(plot_df, aes(x = metric, y = error, color = pair)) +
    geom_boxplot(outlier.shape = NA, fill = "lightgray") +  # Hide default outliers to avoid duplication
    geom_jitter(width = 0.2, alpha = 0.5, aes(color = pair)) +  # Add dots with slight horizontal jitter
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    labs(title = "Error Distribution by Metric, showing Error for each Pair",
         x = "Metric",
         y = "Error (%)") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    scale_color_manual(values = c(
      "Tumour/Tumour" = "#007128",
      "Tumour/Immune" = "#0062c5",
      "Immune/Tumour" = "#f77e3b",
      "Immune/Immune" = "#bb0036"
    )) +
    scale_fill_manual(values = c(
      "Tumour/Tumour" = "#007128",
      "Tumour/Immune" = "#0062c5",
      "Immune/Tumour" = "#f77e3b",
      "Immune/Immune" = "#bb0036"
    ))
  
  return(fig)
}

plot_3D_vs_error_all_metrics_by_slice_box_plot <- function(metric_df_list,
                                                           metrics) {
  
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change and further subset columns of metric_df
    colnames(metric_df)[colnames(metric_df) == metric] <- "value"
    
    # Obtain 3D value
    value_3D <- metric_df[["value"]][metric_df[["slice"]] == as.character(max(as.numeric(metric_df$slice)))]
    
    # Calculate error
    metric_df[["error"]] <- ((metric_df[["value"]] - value_3D) / value_3D) * 100
    
    # Remove value column (only using error now)
    metric_df["value"] <- NULL
    
    # Remove 3D data (integrated into error)
    metric_df <- metric_df[metric_df[["slice"]] != max(as.integer(metric_df$slice)), ]
    
    # Factor slice column so it is from 1, 2, 3, ...
    metric_df$slice <- factor(metric_df$slice, as.character(sort(as.integer(unique(metric_df$slice)))))
    
    
    median_df <- metric_df %>%
      group_by(slice) %>%
      dplyr::summarize(median_error = median(error, na.rm = TRUE), .groups = "drop")
    
    median_df$metric <- metric  
    
    # Add median_df to plot_df
    plot_df <- rbind(plot_df, median_df)
  }
  
  fig <- ggplot(plot_df, aes(x = metric, y = median_error)) +
    geom_boxplot(outlier.shape = NA, fill = "lightgray") +  # Hide default outliers to avoid duplication
    geom_jitter(width = 0.2, alpha = 0.5, color = "#0062c5") +  # Add dots with slight horizontal jitter
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    labs(title = "Error Distribution by Metric, showing Median Error for each Slice",
         x = "Metric",
         y = "Error (%)") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  return(fig)
}


plot_3D_vs_error_all_metrics_by_pair_and_slice_box_plot <- function(metric_df_list,
                                                                    metrics) {
  
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change and further subset columns of metric_df
    colnames(metric_df)[colnames(metric_df) == metric] <- "value"
    
    # Obtain 3D value
    value_3D <- metric_df[["value"]][metric_df[["slice"]] == as.character(max(as.numeric(metric_df$slice)))]
    
    # Calculate error
    metric_df[["error"]] <- ((metric_df[["value"]] - value_3D) / value_3D) * 100
    
    # Remove value column (only using error now)
    metric_df["value"] <- NULL
    
    # Remove 3D data (integrated into error)
    metric_df <- metric_df[metric_df[["slice"]] != max(as.integer(metric_df$slice)), ]
    
    # Add metric column
    metric_df$metric <- metric
    
    # Keep error and metric column
    metric_df <- metric_df[ , c("error", "metric")]
    
    # Add median_df to plot_df
    plot_df <- rbind(plot_df, metric_df)
  }
  
  fig <- ggplot(plot_df, aes(x = metric, y = error)) +
    geom_boxplot(outlier.shape = NA, fill = "lightgray") +  # Hide default outliers to avoid duplication
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    labs(title = "Error Distribution by Metric, showing all Error Values",
         x = "Metric",
         y = "Error (%)") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  return(fig)
}

plot_3D_vs_2D_all_metrics_for_one_pair_and_by_slice_box_plot <- function(metric_df_list,
                                                                         metrics) {
  
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    metric_df <- subset_metric_df(metric,
                                  metric_df,
                                  index = 1)
    
    # Change and further subset columns of metric_df
    colnames(metric_df)[colnames(metric_df) == metric] <- "value"
    
    # Add metric column
    metric_df$metric <- metric
    
    # Keep value, metric and slice column
    metric_df <- metric_df[ , c("value", "metric", "slice")]
    
    # Add median_df to plot_df
    plot_df <- rbind(plot_df, metric_df)
  }
  
  # Add dummy column
  plot_df$dummy <- paste("dummy")
  
  fig <- ggplot(plot_df, aes(x = dummy, y = value)) +
    geom_boxplot(data = plot_df[plot_df$slice != as.character(max(as.numeric(plot_df$slice))), ],
                 outlier.shape = NA, fill = "lightgray") +
    geom_jitter(data = plot_df[plot_df$slice != as.character(max(as.numeric(plot_df$slice))), ],
                width = 0.2, alpha = 0.5, color = "#0062c5") +
    geom_point(data = plot_df[plot_df$slice == as.character(max(as.numeric(plot_df$slice))), ],
               shape = 8, color = "#bb0036", size = 3) +
    facet_wrap(~ metric, scales = "free_y") +  # Facet by metric with independent y-axes
    labs(title = "Value of Metric Distribution by Metric, for one cell pair and for each slice",
         x = "",
         y = "") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.line.x = element_blank(),
    )
  
  
  
  return(fig)
}

plot_3D_vs_error_all_metrics_for_one_pair_and_by_slice_box_plot <- function(metric_df_list,
                                                                            metrics) {
  
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    metric_df <- subset_metric_df(metric,
                                  metric_df,
                                  index = 1)
    
    # Change and further subset columns of metric_df
    colnames(metric_df)[colnames(metric_df) == metric] <- "value"
    
    # Obtain 3D value
    value_3D <- metric_df[["value"]][metric_df[["slice"]] == as.character(max(as.numeric(metric_df$slice)))]
    
    # Calculate error
    metric_df[["error"]] <- ((metric_df[["value"]] - value_3D) / value_3D) * 100
    
    # Remove value column (only using error now)
    metric_df["value"] <- NULL
    
    # Remove 3D data (integrated into error)
    metric_df <- metric_df[metric_df[["slice"]] != max(as.integer(metric_df$slice)), ]
    
    # Add metric column
    metric_df$metric <- metric
    
    # Keep error and metric column
    metric_df <- metric_df[ , c("error", "metric")]
    
    # Add median_df to plot_df
    plot_df <- rbind(plot_df, metric_df)
  }
  
  fig <- ggplot(plot_df, aes(x = metric, y = error)) +
    geom_boxplot(outlier.shape = NA, fill = "lightgray") +  # Hide default outliers to avoid duplication
    geom_jitter(width = 0.2, alpha = 0.5, color = "#0062c5") +  # Add dots with slight horizontal jitter
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    labs(title = "Error Distribution by Metric, for one cell pair and for each slice",
         x = "Metric",
         y = "Error (%)") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  return(fig)
}


## Get plot
metric <- "AMD"
metrics <- c("AMD", "ANC_AUC", "ACIN_AUC", "ANE_AUC", "MS_AUC", "NMS_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")


# This is for a SINGLE metric
fig_3D_vs_2D <- plot_3D_vs_2D(metric_df_list,
                              metric)

# This is for a SINGLE metric
fig_3D_vs_error_by_pair_box_plot <- plot_3D_vs_error_by_pair_box_plot(metric_df_list,
                                                                      metric)

# This is for a SINGLE metric
fig_3D_vs_error_by_slice_box_plot <- plot_3D_vs_error_by_slice_box_plot(metric_df_list,
                                                                        metric)


fig_3D_vs_error_all_metrics_by_pair_box_plot <- 
  plot_3D_vs_error_all_metrics_by_pair_box_plot(metric_df_list,
                                                metrics)

fig_3D_vs_error_all_metrics_by_all_pairs_box_plot <-
  plot_3D_vs_error_all_metrics_by_all_pairs_box_plot(metric_df_list,
                                                     metrics)

fig_3D_vs_error_all_metrics_by_slice_box_plot <- 
  plot_3D_vs_error_all_metrics_by_slice_box_plot(metric_df_list,
                                                 metrics)

fig_3D_vs_error_all_metrics_by_pair_and_slice_box_plot <- 
  plot_3D_vs_error_all_metrics_by_pair_and_slice_box_plot(metric_df_list,
                                                          metrics)

# fig_3D_vs_2D_all_metrics_for_one_pair_and_by_slice_box_plot <-
#   plot_3D_vs_2D_all_metrics_for_one_pair_and_by_slice_box_plot(metric_df_list,
#                                                                            metrics)
# 
# fig_3D_vs_error_all_metrics_for_one_pair_and_by_slice_box_plot <-
#   plot_3D_vs_error_all_metrics_for_one_pair_and_by_slice_box_plot(metric_df_list,
#                                                                               metrics)


### Plotting and upload ------
setwd("~/R/plots/public_data")
pdf(file_name, width = 12, height = 8)

print(fig_3D_vs_2D)
print(fig_3D_vs_error_by_pair_box_plot)
print(fig_3D_vs_error_by_slice_box_plot)
print(fig_3D_vs_error_all_metrics_by_pair_box_plot)
# print(fig_3D_vs_error_all_metrics_by_all_pairs_box_plot)
print(fig_3D_vs_error_all_metrics_by_slice_box_plot)
print(fig_3D_vs_error_all_metrics_by_pair_and_slice_box_plot)
# print(fig_3D_vs_2D_all_metrics_for_one_pair_and_by_slice_box_plot)
# print(fig_3D_vs_error_all_metrics_for_one_pair_and_by_slice_box_plot)

dev.off()
