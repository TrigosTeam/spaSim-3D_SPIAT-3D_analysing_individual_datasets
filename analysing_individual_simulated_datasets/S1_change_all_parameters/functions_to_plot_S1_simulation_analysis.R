# Libraries ----
library(cowplot)
library(ggplot2)

# Read data and set up ----
setwd("~/R/spaSim-3D_SPIAT-3D_analysing_individual_datasets/analysing_individual_simulated_datasets/S1_change_all_parameters")
metric_df_list <- readRDS("S1_metric_df_list_5001_to_10000.RDS")

file_name <- "S1_plots.pdf"

# Functions ----

# Plot 2D vs 3D for each metric and pair
plot_2D_vs_3D_by_metric_and_pair_scatter_plot <- function(metric_df_list,
                                                          metrics) {
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Modify the metric_df so that there is a column of 3D values, and a column of 2D values
    # 2D values are in the second half of the metric_df
    metric_values2D <- metric_df[[metric]][((nrow(metric_df) / 2 ) + 1):nrow(metric_df)]
    metric_df <- metric_df[1:(nrow(metric_df) / 2), ]
    metric_df[["value2D"]] <- metric_values2D
    colnames(metric_df)[colnames(metric_df) == metric] <- "value3D"
    metric_df[["dimension"]] <- NULL
    
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
    plot_df <- rbind(plot_df, metric_df[ , c("value3D", "value2D", "pair", "metric")])
    
  }
  
  pairs <- unique(plot_df$pair)
  
  fig <- ggplot(plot_df, aes(x = value3D, y = value2D)) +
    geom_point(alpha = 0.5, color = "#0062c5") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "#bb0036", linewidth = 1) + # Dotted line of the equation y = x
    labs(title = "Scatterplots showing 2D vs 3D, for each Metric and Pair",
         x = "3D value",
         y = "2D value") +
    facet_wrap(~ interaction(metric, pair), scales = "free", nrow = length(pairs), ncol = length(metrics)) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_blank(),            # remove x-axis numbers
      axis.text.y = element_blank(),            # remove y-axis numbers
      axis.ticks.x = element_blank(),           # remove x-axis ticks
      axis.ticks.y = element_blank()
    )
  
  return(fig)
}

plot_error_vs_3D_by_metric_and_pair_scatter_plot <- function(metric_df_list,
                                                             metrics) {
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Modify the metric_df so that there is a column of 3D values, and a column of 2D values
    # 2D values are in the second half of the metric_df
    metric_values2D <- metric_df[[metric]][((nrow(metric_df) / 2 ) + 1):nrow(metric_df)]
    metric_df <- metric_df[1:(nrow(metric_df) / 2), ]
    metric_df[["value2D"]] <- metric_values2D
    colnames(metric_df)[colnames(metric_df) == metric] <- "value3D"
    metric_df[["dimension"]] <- NULL
    
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
    plot_df <- rbind(plot_df, metric_df[ , c("value3D", "value2D", "pair", "metric")])
    
  }
  
  # Get error column
  plot_df$error <- (plot_df$value2D - plot_df$value3D) / plot_df$value3D * 100
  
  pairs <- unique(plot_df$pair)
  
  fig <- ggplot(plot_df, aes(x = value3D, y = error)) +
    geom_point(alpha = 0.5, color = "#0062c5") +
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    labs(title = "Scatterplots showing Error vs 3D, for each Metric and Pair",
         x = "3D value",
         y = "Error (%)") +
    facet_grid(pair ~ metric, scales = "free") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  return(fig)
}

plot_error_vs_3D_by_metric_and_pair_box_plot <- function(metric_df_list,
                                                         metrics) {
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Modify the metric_df so that there is a column of 3D values, and a column of 2D values
    # 2D values are in the second half of the metric_df
    metric_values2D <- metric_df[[metric]][((nrow(metric_df) / 2 ) + 1):nrow(metric_df)]
    metric_df <- metric_df[1:(nrow(metric_df) / 2), ]
    metric_df[["value2D"]] <- metric_values2D
    colnames(metric_df)[colnames(metric_df) == metric] <- "value3D"
    metric_df[["dimension"]] <- NULL
    
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
    plot_df <- rbind(plot_df, metric_df[ , c("value3D", "value2D", "pair", "metric")])
    
  }
  
  # Get error column
  plot_df$error <- (plot_df$value2D - plot_df$value3D) / plot_df$value3D * 100
  
  pairs <- unique(plot_df$pair)
  
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
      "A/A" = "#bb0036",
      "A/B" = "#0062c5",
      "B/A" = "#007128",
      "B/B" = "#f77e3b"
    )) +
    scale_fill_manual(values = c(
      "A/A" = "#bb0036",
      "A/B" = "#0062c5",
      "B/A" = "#007128",
      "B/B" = "#f77e3b"
    ))
  
  return(fig)
}


# Running the functions ----

metrics <- c("AMD", "ACIN_AUC", "ACINP_AUC", "AE_AUC", "MS_AUC", "NMS_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")


fig_2D_vs_3D_by_metric_and_pair_scatter_plot <- plot_2D_vs_3D_by_metric_and_pair_scatter_plot(metric_df_list,
                                                                                          metrics)

fig_error_vs_3D_by_metric_and_pair_scatter_plot <- plot_error_vs_3D_by_metric_and_pair_scatter_plot(metric_df_list,
                                                                                               metrics)

fig_error_vs_3D_by_metric_and_pair_box_plot <- plot_error_vs_3D_by_metric_and_pair_box_plot(metric_df_list,
                                                                                            metrics)

setwd("~/R/plots/S1")
pdf(file_name, width = 18, height = 10)

print(fig_2D_vs_3D_by_metric_and_pair_scatter_plot)
print(fig_error_vs_3D_by_metric_and_pair_scatter_plot)
print(fig_error_vs_3D_by_metric_and_pair_box_plot)

dev.off()
