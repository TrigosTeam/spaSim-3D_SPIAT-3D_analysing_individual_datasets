# Libraries ----
library(cowplot)
library(ggplot2)
library(S4Vectors)
library(stringr)
library(dplyr)

# Read data and set up ----
setwd("~/R/S1_data")

metric_df_list <- readRDS("S1_metric_df_list.RDS")

file_name <- "S1_plots.pdf"

# Functions ----

# Plot 2D vs 3D for each metric and pair
plot_2D_vs_3D_by_metric_and_pair_scatter_plot <- function(metric_df_list,
                                                          metrics) {
  plot_df <- data.frame()
  pval_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Modify the metric_df so that there is a column of 3D values, and a column of 2D values
    metric_values2D <- metric_df[[metric]][metric_df[["dimension"]] == "2D"]
    metric_df <- metric_df[metric_df[["dimension"]] == "3D", ]
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
    
    # Add to pval_df
    p_values <- c()
    for (pair in unique(metric_df$pair)) {
      
      # Ignore when pair is invalid
      if (sum(is.finite(metric_df[["value3D"]][metric_df$pair == pair])) == 0) {
        next
      } 
      
      wilcox_test  <- wilcox.test(metric_df[["value3D"]][metric_df$pair == pair], 
                                  metric_df[["value2D"]][metric_df$pair == pair], 
                                  paired = TRUE)
      p_value <- wilcox_test$p.value
      
      # Format p-value
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
      
      p_values <- c(p_values, p_value)
    }
    pval_df <- rbind(pval_df, 
                     data.frame(metric = metric,
                                pair = unique(metric_df$pair),
                                p_value = p_values))
  }
  
  pairs <- unique(plot_df$pair)

  fig <- ggplot(plot_df, aes(x = value3D, y = value2D)) +
    geom_point(alpha = 0.25, color = "#0062c5", size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "#bb0036", linewidth = 1) + # Dotted line of the equation y = x
    labs(title = "Scatterplots showing 2D vs 3D, for each Metric and Pair",
         x = "3D value",
         y = "2D value") +
    facet_wrap(~ interaction(metric, pair), scales = "free", nrow = length(pairs), ncol = length(metrics)) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(size = 4),  # make x-axis text smaller
      axis.text.y = element_text(size = 4)   # make y-axis text smaller
    ) +
    # Add p-value text
    geom_text(
      data = pval_df,
      aes(
        x = Inf, y = -Inf,   # bottom-right corner
        label = paste0("p: ", p_value)
      ),
      inherit.aes = FALSE,
      hjust = 1.1, vjust = -0.5,
      size = 2,
      color = "black"
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
    metric_values2D <- metric_df[[metric]][metric_df[["dimension"]] == "2D"]
    metric_df <- metric_df[metric_df[["dimension"]] == "3D", ]
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
    geom_point(alpha = 0.25, color = "#0062c5", size = 1) +
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    labs(title = "Scatterplots showing Error vs 3D, for each Metric and Pair",
         x = "3D value",
         y = "Error (%)") +
    facet_wrap(~ interaction(metric, pair), scales = "free", nrow = length(pairs), ncol = length(metrics)) +  
    scale_y_continuous(limits = c(-100, 500)) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(size = 6),  # make x-axis text smaller
      axis.text.y = element_text(size = 6)   # make y-axis text smaller
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
    metric_values2D <- metric_df[[metric]][metric_df[["dimension"]] == "2D"]
    metric_df <- metric_df[metric_df[["dimension"]] == "3D", ]
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
    geom_boxplot(fill = "lightgray") +
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    labs(title = "Error Distribution by Metric, showing Error for each Pair",
         x = "Metric",
         y = "Error (%)") +
    scale_y_continuous(limits = c(-100, 400)) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    scale_color_manual(values = c(
      "A/A" = "#bb0036",
      "A/B" = "#0062c5",
      "B/A" = "#007128",
      "B/B" = "#f77e3b"
    ))
  
  return(fig)
}




# Running the functions ----

# # Subset...
# random_numbers <- sample(1:40000, 1000, replace = FALSE)
# metric_df_list_subset <- list()
# for (metric in names(metric_df_list)) {
#   metric_df <- metric_df_list[[metric]]
#   metric_df_subset <- metric_df[c(random_numbers, random_numbers + 40000), ]
#   metric_df_list_subset[[metric]] <- metric_df_subset
# }

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
