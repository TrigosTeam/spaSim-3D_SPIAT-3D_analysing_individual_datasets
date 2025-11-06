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

### Read data -----
setwd("~/R/SPIAT-3D_benchmarking/public_3D_data_analysis")
metric_df_list <- readRDS("merfish_mouse_hypothalamus_metric_df_list.RDS")

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


### Plot analysis -----

## Functions to plot
plot_3D_vs_2D <- function(metric_df_list,
                          metric) {
  
  # Get metric_df for current metric
  metric_df <- metric_df_list[[metric]]
  
  # Change and further subset columns of metric_df_subset
  colnames(metric_df)[colnames(metric_df) == metric] <- "value"
  
  # Add dummy column
  metric_df$dummy <- paste("dummy")
  
  fig <- ggplot(metric_df, aes(x = dummy, y = value)) +
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
      axis.line.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.line.y = element_blank()
    ) +
    facet_grid(reference ~ target, scales = "free_y")
  
  
  return(fig)
}

plot_3D_vs_error_by_cell_combination_box_plot <- function(metric_df_list,
                                                          metric) {
  
  # Get metric_df for current metric
  metric_df <- metric_df_list[[metric]]
  
  # Change and further subset columns of metric_df
  colnames(metric_df)[colnames(metric_df) == metric] <- "value"
  
  # Obtain 3D value
  value_3D <- metric_df[["value"]][metric_df[["slice"]] == max(metric_df$slice)]
  
  # Calculate error
  metric_df[["error"]] <- ((metric_df[["value"]] - value_3D) / value_3D) * 100
  
  # Remove value column (only using error now)
  metric_df["value"] <- NULL
  
  # Remove 3D data (integrated into error)
  metric_df <- metric_df[metric_df[["slice"]] != max(as.integer(metric_df$slice)), ]
  
  # Add reference-target column
  metric_df$reference_target <- paste(metric_df$reference, metric_df$target, sep = "/")
  
  fig <- ggplot(metric_df, aes(x = reference_target, y = error)) +
    geom_boxplot(outlier.shape = NA, fill = "lightgray") +  # Hide default outliers to avoid duplication
    # geom_jitter(width = 0.2, alpha = 0.5, color = "#0062c5") +  # Add dots with slight horizontal jitter
    labs(title = "Error Distribution by Reference/Target combination",
         x = "Reference/Target combination",
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
  value_3D <- metric_df[["value"]][metric_df[["slice"]] == max(metric_df$slice)]
  
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
    # geom_jitter(width = 0.2, alpha = 0.5, color = "#0062c5") +  # Add dots with slight horizontal jitter
    labs(title = "Error Distribution by Slice",
         x = "Slice index",
         y = "Error (%)") +
    lims(y = c(-200, 200)) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)  # Rotate x-axis labels vertically
    )
  
  return(fig)
}


plot_3D_vs_error_all_metrics_by_cell_combination_box_plot <- function(metric_df_list,
                                                                      metrics) {
  
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change and further subset columns of metric_df
    colnames(metric_df)[colnames(metric_df) == metric] <- "value"
    
    # Obtain 3D value
    value_3D <- metric_df[["value"]][metric_df[["slice"]] == max(metric_df$slice)]
    
    # Calculate error
    metric_df[["error"]] <- ((metric_df[["value"]] - value_3D) / value_3D) * 100
    
    # Remove value column (only using error now)
    metric_df["value"] <- NULL
    
    # Remove 3D data (integrated into error)
    metric_df <- metric_df[metric_df[["slice"]] != max(as.integer(metric_df$slice)), ]
    
    if (!(metric %in% c("EBSAC", "EBP_AUC"))) {
      # Add reference-target column
      metric_df$reference_target <- paste(metric_df$reference, metric_df$target, sep = "/")
      
      # Extract median values for error for each reference-target
      median_df <- metric_df %>%
        group_by(reference_target) %>%
        dplyr::summarize(median_error = median(error, na.rm = TRUE), .groups = "drop")
      
      median_df$reference_target <- NULL
    }
    else {
      # Extract median values for error for each reference-target
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
    labs(title = "Error Distribution by metric, showing median error for each Reference/Target combination",
         x = "Metric",
         y = "Error (%)") +
    lims(y = c(-200, 200)) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
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
    value_3D <- metric_df[["value"]][metric_df[["slice"]] == max(metric_df$slice)]
    
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
    labs(title = "Error Distribution by Metric, showing median error for each Slice",
         x = "Metric",
         y = "Error (%)") +
    lims(y = c(-200, 200)) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  return(fig)
}


plot_3D_vs_error_all_metrics_by_combination_and_slice_box_plot <- function(metric_df_list,
                                                                           metrics) {
  
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change and further subset columns of metric_df
    colnames(metric_df)[colnames(metric_df) == metric] <- "value"
    
    # Obtain 3D value
    value_3D <- metric_df[["value"]][metric_df[["slice"]] == max(metric_df$slice)]
    
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
    labs(title = "Error Distribution by Metric, showing all error values",
         x = "Metric",
         y = "Error (%)") +
    lims(y = c(-200, 200)) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  return(fig)
}

## Get plot
metrics <- c("AMD", "ACIN_AUC", "ACINP_AUC", "AE_AUC", "MS_AUC", "NMS_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")


# This is for a SINGLE metric
fig_3D_vs_2D <- plot_3D_vs_2D(metric_df_list,
                              metric)

# This is for a SINGLE metric
fig_3D_vs_error_by_combination_box_plot <- plot_3D_vs_error_by_cell_combination_box_plot(metric_df_list,
                                                                                         metric)

# This is for a SINGLE metric
fig_3D_vs_error_by_slice_box_plot <- plot_3D_vs_error_by_slice_box_plot(metric_df_list,
                                                                        metric)


fig_3D_vs_error_all_metrics_by_cell_combination_box_plot <- 
  plot_3D_vs_error_all_metrics_by_cell_combination_box_plot(metric_df_list,
                                                            metrics)

fig_3D_vs_error_all_metrics_by_combination_and_slice_box_plot <- 
  plot_3D_vs_error_all_metrics_by_combination_and_slice_box_plot(metric_df_list,
                                                                 metrics)


fig_3D_vs_error_all_metrics_by_slice_box_plot <- 
  plot_3D_vs_error_all_metrics_by_slice_box_plot(metric_df_list,
                                                 metrics)


methods::show(fig_3D_vs_2D)
methods::show(fig_3D_vs_error_by_combination_box_plot)
methods::show(fig_3D_vs_error_by_slice_box_plot)
methods::show(fig_3D_vs_error_all_metrics_by_cell_combination_box_plot)
methods::show(fig_3D_vs_error_all_metrics_by_slice_box_plot)
methods::show(fig_3D_vs_error_all_metrics_by_combination_and_slice_box_plot)



setwd("~/R/plots/public_data")
pdf("merfish_mouse_hypothalamus.pdf", width = 12, height = 8)

print(fig_3D_vs_2D)
print(fig_3D_vs_error_by_combination_box_plot)
print(fig_3D_vs_error_by_slice_box_plot)
print(fig_3D_vs_error_all_metrics_by_cell_combination_box_plot)
print(fig_3D_vs_error_all_metrics_by_slice_box_plot)
print(fig_3D_vs_error_all_metrics_by_combination_and_slice_box_plot)

dev.off()
