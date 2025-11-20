# Libraries ----
library(cowplot)
library(ggplot2)
library(S4Vectors)
library(stringr)
library(dplyr)
library(scales)

# Read data and set up ----
setwd("~/R/spaSim-3D_SPIAT-3D_analysing_individual_datasets/analysing_individual_simulated_datasets/S2_change_each_parameter")
metric_df_list <- readRDS("S2_metric_df_list.RDS")
parameters_df <- readRDS("S2_parameters_df.RDS")


# Functions -----
# Subset metric_df_list and parameters_df
add_pair_to_metric_df <- function(metric_df, metric) {
  
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
  return(metric_df)
}

get_parameters <- function(arrangement, shape) {
  parameters <- c()
  
  # Add other parameters specific to arrangement and shape
  if (arrangement == "mixed") {
    parameters <- c(parameters, "cluster_prop_A")
  }
  else if (arrangement == "ringed") {
    parameters <- c(parameters, "ring_width_factor")        
  }
  else if (arrangement == "separated") {
    parameters <- c(parameters, "distance")
  }
  
  if (shape == "ellipsoid") {
    parameters <- c(parameters, "E_volume")
  }
  else if (shape == "network") {
    parameters <- c(parameters, "N_width")
  }
  
  # All plots include bg_prop_A and bg_prop_B as parameters
  parameters <- c(parameters, "bg_prop_A", "bg_prop_B")
  
  return(parameters)
}

plot_3D_vs_parameters_for_non_gradient_metrics_scatter_plot <- function(metric_df_list, parameters_df, metric) {
  
  fig_list <- list()
  
  # Define parameters
  arrangements <- c("mixed", "ringed", "separated")
  shapes <- c("ellipsoid", "network")
  
  # Subset the metric_df_list
  metric_df <- metric_df_list[[metric]]
  metric_df <- metric_df[metric_df$slice == 0, ] # slice == 0 refers to 3D values
  
  # Add 'pair' column to metric_df
  metric_df <- add_pair_to_metric_df(metric_df, metric)
  
  pairs <- unique(metric_df$pair)
  
  # Add to parameters_df
  parameters_df$distance <- 450 - parameters_df$cluster1_x_coord # 450 is the x-coordinate of cluster2
  parameters_df$E_volume <- (4 / 3) * pi * parameters_df$E_radius_x * parameters_df$E_radius_y * parameters_df$E_radius_z
  
  # Update parameters_df
  parameters_df$variable_parameter[parameters_df$variable_parameter == "E_radius_x"] <- "E_volume"
  parameters_df$variable_parameter[parameters_df$variable_parameter == "cluster1_x_coord"] <- "distance"
  
  # Combine metric_df with parameters_df
  metric_df <- cbind(metric_df, parameters_df[rep(1:nrow(parameters_df), each = length(pairs)), ])
  
  for (arrangement in arrangements) {
    for (shape in shapes) {
      
      # Get a merged arrangement shape variable
      arrangement_shape <- paste(arrangement, shape, sep = "_")
      fig_list[[arrangement_shape]] <- list()
      
      # Subset metric_df for arrangement and shape
      metric_arrangement_shape_df <- metric_df[metric_df$arrangement == arrangement & metric_df$shape == shape, ]
      
      # Get parameters
      parameters <- get_parameters(arrangement, shape)
      
      for (pair in pairs) {
        # Subset for pair
        metric_arrangement_shape_pair_df <- metric_arrangement_shape_df[metric_arrangement_shape_df$pair == pair, ]
        fig_list[[arrangement_shape]][[pair]] <- list()
        
        for (parameter in parameters) {
          # Subset for parameter
          plot_df <- metric_arrangement_shape_pair_df[metric_arrangement_shape_pair_df$variable_parameter == parameter, ]
          
          # Get correlation and p-value
          correlation_test <- cor.test(plot_df[[parameter]], plot_df[[metric]], method = "spearman")
          correlation <- correlation_test$estimate
          p_value <- correlation_test$p.value
          
          # Format correlation and p-value
          formatCustomSci <- function(x) {
            x_sci <- str_split_fixed(formatC(x, format = "e"), "e", 2)
            alpha <- round(as.numeric(x_sci[ , 1]), 1)
            power <- as.integer(x_sci[ , 2])
            paste(alpha, power, sep = "e")
          }
          if (p_value == 0) p_value <- 2.2e-308
          if (0 < p_value && p_value < 1e-3)  {
            p_value <- formatCustomSci(p_value)
          }
          else {
            p_value  <- round(p_value, 3)
          }
          
          fig <- ggplot(plot_df, aes_string(parameter, metric)) +
            geom_point(size = 1) +
            ggtitle(paste("r:", round(correlation, 3), ", p:", p_value)) +
            theme_minimal() +
            theme(
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              plot.title = element_text(size = 10),
              axis.text.x = element_text(size = 8),  # make x-axis text smaller
              axis.text.y = element_text(size = 8)   # make y-axis text smaller
            ) +
            scale_x_continuous(breaks = pretty_breaks(n = 3)) + 
            scale_y_continuous(breaks = pretty_breaks(n = 3)) 
          
          fig_list[[arrangement_shape]][[pair]][[parameter]] <- fig
        }
      }
    }
  }
  
  arrangement_shape_figs <- list()
  
  for (shape in shapes) {  
    for (arrangement in arrangements) {

      # Get a merged arrangement shape variable
      arrangement_shape <- paste(arrangement, shape, sep = "_")
      
      pair_figs <- list()
      
      for (pair in pairs) {

        parameters <- get_parameters(arrangement, shape)
        
        pair_fig <- plot_grid(plotlist = fig_list[[arrangement_shape]][[pair]],
                              ncol = length(fig_list[[arrangement_shape]][[pair]]))
        title <- ggdraw() + draw_label(paste("pair:", pair))
        pair_fig <- plot_grid(title, pair_fig, ncol = 1, rel_heights = c(0.1, 1))
        
        pair_figs[[pair]] <- pair_fig
      }
      arrangement_shape_fig <- plot_grid(plotlist = pair_figs, ncol = 1)
      
      title <- ggdraw() + draw_label(paste("arrangement-shape: ", arrangement, "-", shape, sep = ""), fontface = "bold")
      arrangement_shape_fig <- plot_grid(title, arrangement_shape_fig, ncol = 1, rel_heights = c(0.02, 1))  + 
        theme(plot.margin = margin(10, 10, 10, 10),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1))  
      
      arrangement_shape_figs[[arrangement_shape]] <- arrangement_shape_fig
    }
  }
  
  fig <- plot_grid(plotlist = arrangement_shape_figs,
                   nrow = length(shapes),
                   ncol = length(arrangements))
  
  return(fig)
}

plot_3D_vs_parameters_for_gradient_metrics_line_graph <- function(metric_df_list, parameters_df, metric) {
  
}

plot_3D_and_2D_vs_parameters_for_non_gradient_metrics_scatter_plot <- function(metric_df_list, parameters_df, metric) {
  fig_list <- list()
  
  # Define parameters
  arrangements <- c("mixed", "ringed", "separated")
  shapes <- c("ellipsoid", "network")
  
  # Subset the metric_df_list
  metric_df <- metric_df_list[[metric]]
  
  # Add 'pair' column to metric_df
  metric_df <- add_pair_to_metric_df(metric_df, metric)
  
  pairs <- unique(metric_df$pair)
  
  # Make 'slice' column categorical
  metric_df$slice <- as.character(metric_df$slice)
  
  # Add to parameters_df
  parameters_df$distance <- 450 - parameters_df$cluster1_x_coord # 450 is the x-coordinate of cluster2
  parameters_df$E_volume <- (4 / 3) * pi * parameters_df$E_radius_x * parameters_df$E_radius_y * parameters_df$E_radius_z
  
  # Update parameters_df
  parameters_df$variable_parameter[parameters_df$variable_parameter == "E_radius_x"] <- "E_volume"
  parameters_df$variable_parameter[parameters_df$variable_parameter == "cluster1_x_coord"] <- "distance"
  
  # Combine metric_df with parameters_df, use '4' as 3 slices + 3D value
  metric_df <- cbind(metric_df, do.call(rbind, replicate(4, parameters_df[rep(1:nrow(parameters_df), each = length(pairs)), ], simplify = FALSE)))

  for (arrangement in arrangements) {
    for (shape in shapes) {
      
      # Get a merged arrangement shape variable
      arrangement_shape <- paste(arrangement, shape, sep = "_")
      fig_list[[arrangement_shape]] <- list()
      
      # Subset metric_df for arrangement and shape
      metric_arrangement_shape_df <- metric_df[metric_df$arrangement == arrangement & metric_df$shape == shape, ]
      
      # Get parameters
      parameters <- get_parameters(arrangement, shape)
      
      for (pair in pairs) {
        # Subset for pair
        metric_arrangement_shape_pair_df <- metric_arrangement_shape_df[metric_arrangement_shape_df$pair == pair, ]
        fig_list[[arrangement_shape]][[pair]] <- list()
        
        for (parameter in parameters) {
          # Subset for parameter
          plot_df <- metric_arrangement_shape_pair_df[metric_arrangement_shape_pair_df$variable_parameter == parameter, ]
          
          fig <- ggplot(plot_df, aes_string(parameter, metric, color = "slice")) +
            geom_point(size = 1) +
            theme_minimal() +
            theme(
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              plot.title = element_text(size = 10),
              axis.text.x = element_text(size = 8),  # make x-axis text smaller
              axis.text.y = element_text(size = 8),   # make y-axis text smaller
              legend.position = "none"
            ) +
            scale_x_continuous(breaks = pretty_breaks(n = 3)) + 
            scale_y_continuous(breaks = pretty_breaks(n = 3)) 
          
          fig_list[[arrangement_shape]][[pair]][[parameter]] <- fig
        }
      }
    }
  }
  
  arrangement_shape_figs <- list()
  
  for (shape in shapes) {  
    for (arrangement in arrangements) {
      
      # Get a merged arrangement shape variable
      arrangement_shape <- paste(arrangement, shape, sep = "_")
      
      pair_figs <- list()
      
      for (pair in pairs) {
        
        parameters <- get_parameters(arrangement, shape)
        
        pair_fig <- plot_grid(plotlist = fig_list[[arrangement_shape]][[pair]],
                              ncol = length(fig_list[[arrangement_shape]][[pair]]))
        title <- ggdraw() + draw_label(paste("pair:", pair))
        pair_fig <- plot_grid(title, pair_fig, ncol = 1, rel_heights = c(0.1, 1))
        
        pair_figs[[pair]] <- pair_fig
      }
      arrangement_shape_fig <- plot_grid(plotlist = pair_figs, ncol = 1)
      
      title <- ggdraw() + draw_label(paste("arrangement-shape: ", arrangement, "-", shape, sep = ""), fontface = "bold")
      arrangement_shape_fig <- plot_grid(title, arrangement_shape_fig, ncol = 1, rel_heights = c(0.02, 1))  + 
        theme(plot.margin = margin(10, 10, 10, 10),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1))  
      
      arrangement_shape_figs[[arrangement_shape]] <- arrangement_shape_fig
    }
  }
  
  fig <- plot_grid(plotlist = arrangement_shape_figs,
                   nrow = length(shapes),
                   ncol = length(arrangements))
  
  return(fig) 
}

plot_error_vs_parameters_for_non_gradient_metrics_scatter_plot <- function(metric_df_list, parameters_df, metric) {
  
}

plot_2D_vs_slice_for_non_gradient_metrics_violin_plot <- function(metric_df_list, parameters_df, metric) {
  
}

plot_3D_and_2D_vs_slice_for_non_gradient_metrics_violin_plot <- function(metric_df_list, parameters_df, metric) {
  
}



# Running the functions ----

metrics <- c("AMD", "ACIN_AUC", "ACINP_AUC", "AE_AUC", "MS_AUC", "NMS_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")

fig_3D_vs_parameters_for_non_gradient_metrics_scatter_plot <- plot_3D_vs_parameters_for_non_gradient_metrics_scatter_plot(metric_df_list,
                                                                                                                          parameters_df,
                                                                                                                          "AMD")


setwd("~/R/plots/S2")
pdf("fig_3D_vs_parameters_for_non_gradient_metrics_scatter_plot.pdf", width = 25, height = 25)

print(fig_3D_vs_parameters_for_non_gradient_metrics_scatter_plot)

dev.off()



fig_3D_and_2D_vs_parameters_for_non_gradient_metrics_scatter_plot <- plot_3D_and_2D_vs_parameters_for_non_gradient_metrics_scatter_plot(metric_df_list,
                                                                                                                                        parameters_df,
                                                                                                                                        "AMD")

setwd("~/R/plots/S2")
pdf("fig_3D_and_2D_vs_parameters_for_non_gradient_metrics_scatter_plot.pdf", width = 25, height = 25)

print(fig_3D_and_2D_vs_parameters_for_non_gradient_metrics_scatter_plot)

dev.off()

