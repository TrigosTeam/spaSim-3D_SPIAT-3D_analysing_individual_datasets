# Libraries ----
library(cowplot)
library(ggplot2)

# Read data and set up ----
setwd("~/R/spaSim-3D_SPIAT-3D_analysing_individual_datasets/analysing_individual_simulated_datasets/S2_change_each_parameter")
metric_df_lists <- list(
  a = readRDS("S2_metric_df_list_1_to_6000.RDS"),
  b = readRDS("S2_metric_df_list_6001_to_12000.RDS"),
  c = readRDS("S2_metric_df_list_12001_to_18000.RDS"),
  d = readRDS("S2_metric_df_list_18001_to_24000.RDS")
)


metrics <- c("AMD", "ACIN_AUC", "ACINP_AUC", "AE_AUC", "MS_AUC", "NMS_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")

metric_df_list <- list()

# Combine datasets together
for (metric in metrics) {
  
  metric_df <- data.frame()
  i <- 0
  
  for (temp_metric_df_list in metric_df_lists) {
    
    temp_metric_df <- temp_metric_df_list[[metric]]
  
    temp_metric_df[["simulation"]] <- as.numeric(temp_metric_df[["simulation"]]) + i * 6000
      
    metric_df <- rbind(metric_df, temp_metric_df)
  
    i <- i + 1
  }
  metric_df_list[[metric]] <- metric_df
}

# Functions -----
plot_3D_vs_parameters_for_non_gradient_metrics_scatter_plot <- function(metric_df_list, metric) {
  
}

plot_3D_vs_parameters_for_gradient_metrics_line_graph <- function(metric_df_list, metric) {
  
}

plot_3D_and_2D_vs_parameters_for_non_gradient_metrics_scatter_plot <- function(metric_df_list, metric) {
  
}

plot_error_vs_parameters_for_non_gradient_metrics_scatter_plot <- function(metric_df_list, metric) {
  
}

plot_2D_vs_slice_for_non_gradient_metrics_violin_plot <- function(metric_df_list, metric) {
  
}

plot_3D_and_2D_vs_slice_for_non_gradient_metrics_violin_plot <- function(metric_df_list, metric) {
  
}



# Running the functions ----

metrics <- c("AMD", "ACIN_AUC", "ACINP_AUC", "AE_AUC", "MS_AUC", "NMS_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")