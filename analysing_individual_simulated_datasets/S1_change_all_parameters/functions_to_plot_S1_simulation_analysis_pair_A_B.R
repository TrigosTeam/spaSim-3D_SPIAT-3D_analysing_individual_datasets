# Libraries ----
library(cowplot)
library(ggplot2)
library(S4Vectors)
library(stringr)
library(dplyr)

# Read data and set up ----
setwd("~/R/S1_data")

metric_df_list <- readRDS("S1_metric_df_list.RDS")

parameters_df <- readRDS("parameters_df.RDS")

# Functions ----

# Plot 2D vs 3D for each metric and pair for a random slice
plot_2D_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot <- function(metric_df_list,
                                                                               metrics) {
  plot_df <- data.frame()
  pval_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change slice column to numeric for safety
    metric_df$slice <- as.numeric(metric_df$slice)
    
    # Modify the metric_df so that there is a column of 3D values, and a column of 2D values
    # Choose a random slice for 2D values
    metric_values2D <- metric_df[[metric]][metric_df[["slice"]] == sample(unique(metric_df[["slice"]]), 1)]
    # Slice == 0 is the 3D value
    metric_df <- metric_df[metric_df[["slice"]] == 0, ]
    metric_df[["value2D"]] <- metric_values2D
    colnames(metric_df)[colnames(metric_df) == metric] <- "value3D"
    
    if (metric %in% c("EBSAC", "EBP_AUC")) {
      # For EBSAC and EBP_AUC, assume pair is the same as cell_types for consistency
      metric_df$pair <- gsub(',', '/', metric_df$cell_types)
    }
    else if (metric %in% c("ANE_AUC")) {
      # For ANE_AUC, assume pair is the same as target_cell_type for consistency (as target is of form A,B already)
      metric_df$pair <- gsub(',', '/', metric_df$target)
    }
    else {
      # Add reference-target column
      metric_df$pair <- paste(metric_df$reference, metric_df$target, sep = "/")
    }
    
    # Subset for A/B pair
    metric_df <- metric_df[metric_df$pair == "A/B", ]
    
    metric_df$metric <- metric  
    
    # Add metric_df to plot_df
    plot_df <- rbind(plot_df, metric_df[ , c("value3D", "value2D", "metric")])
    
    # Add to pval_df
    p_values <- c()
    
    # Ignore when pair is invalid
    if (sum(is.finite(metric_df[["value3D"]])) == 0) {
      next
    } 
    
    wilcox_test  <- wilcox.test(metric_df[["value3D"]], 
                                metric_df[["value2D"]], 
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
    
    pval_df <- rbind(pval_df, 
                     data.frame(metric = metric,
                                p_value = p_values))
  }
  
  sci_if_big <- function(x) { 
    out <- ifelse( 
      abs(x) >= 5e3, 
      scales::scientific_format()(x), 
      scales::number_format(big.mark = "", decimal.mark = ".")(x) ) # Remove leading zeros in exponent (e.g., 3e03 → 3e3) 
    
    out <- gsub("e([+-])0+", "e\\1", out) 
    
    out
  } 
  
  fig <- ggplot(plot_df, aes(x = value3D, y = value2D)) +
    geom_point(alpha = 0.25, color = "#0062c5", size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                color = "#bb0036", linewidth = 1) +
    labs(
      title = "Scatterplots showing 2D vs 3D, for each metric, a random slice and pair A/B",
      x = "3D value",
      y = "2D value"
    ) +
    facet_wrap(~ interaction(metric), scales = "free", ncol = length(metrics)) +
    
    scale_x_continuous(n.breaks = 3, labels = sci_if_big) + 
    scale_y_continuous(n.breaks = 3, labels = sci_if_big) +
    
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      
      # Axis tick labels
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      
      # Axis titles
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      
      # Plot title
      plot.title = element_text(size = 15),
      
      # Facet strip titles
      strip.text = element_text(size = 12)
    ) +
    geom_text(
      data = pval_df,
      aes(x = Inf, y = -Inf, label = paste0("p: ", p_value)),
      inherit.aes = FALSE,
      hjust = 1.1, vjust = -0.5,
      size = 4,   # p-value font size
      color = "black"
    )
  
  
  return(fig)
}

# Plot error vs 3D for each metric and pair for a random slice
plot_error_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot <- function(metric_df_list,
                                                                                  metrics) {
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Choose a random slice for 2D values
    metric_values2D <- metric_df[[metric]][metric_df[["slice"]] == sample(unique(metric_df[["slice"]]), 1)]
    # Slice == 0 is the 3D value
    metric_df <- metric_df[metric_df[["slice"]] == 0, ]
    metric_df[["value2D"]] <- metric_values2D
    colnames(metric_df)[colnames(metric_df) == metric] <- "value3D"
    
    if (metric %in% c("EBSAC", "EBP_AUC")) {
      # For EBSAC and EBP_AUC, assume pair is the same as cell_types for consistency
      metric_df$pair <- gsub(',', '/', metric_df$cell_types)
    }
    else if (metric %in% c("ANE_AUC")) {
      # For ANE_AUC, assume pair is the same as target_cell_type for consistency (as target is of form A,B already)
      metric_df$pair <- gsub(',', '/', metric_df$target)
    }
    else {
      # Add reference-target column
      metric_df$pair <- paste(metric_df$reference, metric_df$target, sep = "/")
    }
    
    # Subset for A/B pair
    metric_df <- metric_df[metric_df$pair == "A/B", ]
    
    metric_df$metric <- metric  
    
    # Add metric_df to plot_df
    plot_df <- rbind(plot_df, metric_df[ , c("value3D", "value2D", "pair", "metric")])
    
  }
  
  # Get error column
  plot_df$error <- (plot_df$value2D - plot_df$value3D) / plot_df$value3D * 100
  
  fig <- ggplot(plot_df, aes(x = value3D, y = error)) +
    geom_point(alpha = 0.25, color = "#0062c5", size = 1) +
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    labs(title = "Scatterplots showing Error vs 3D, for each Metric, for a random slice and pair A/B",
         x = "3D value",
         y = "Error (%)") +
    facet_wrap(~ interaction(metric, pair), scales = "free", ncol = length(metrics)) +  
    scale_y_continuous(limits = c(-100, 500)) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(size = 6),  # make x-axis text smaller
      axis.text.y = element_text(size = 6)   # make y-axis text smaller
    )
  
  return(fig)
}

# Plot error vs 3D for each metric and pair for a random slice
plot_error_vs_3D_by_metric_and_pair_A_B_for_random_slice_box_plot <- function(metric_df_list,
                                                                              metrics) {
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Choose a random slice for 2D values
    metric_values2D <- metric_df[[metric]][metric_df[["slice"]] == sample(unique(metric_df[["slice"]]), 1)]
    # Slice == 0 is the 3D value
    metric_df <- metric_df[metric_df[["slice"]] == 0, ]
    metric_df[["value2D"]] <- metric_values2D
    colnames(metric_df)[colnames(metric_df) == metric] <- "value3D"
    
    if (metric %in% c("EBSAC", "EBP_AUC")) {
      # For EBSAC and EBP_AUC, assume pair is the same as cell_types for consistency
      metric_df$pair <- gsub(',', '/', metric_df$cell_types)
    }
    else if (metric %in% c("ANE_AUC")) {
      # For ANE_AUC, assume pair is the same as target_cell_type for consistency (as target is of form A,B already)
      metric_df$pair <- gsub(',', '/', metric_df$target)
    }
    else {
      # Add reference-target column
      metric_df$pair <- paste(metric_df$reference, metric_df$target, sep = "/")
    }
    
    # Subset for A/B pair
    metric_df <- metric_df[metric_df$pair == "A/B", ]
    
    metric_df$metric <- metric  
    
    # Add metric_df to plot_df
    plot_df <- rbind(plot_df, metric_df[ , c("value3D", "value2D", "pair", "metric")])
    
  }
  
  # Get error column
  plot_df$error <- (plot_df$value2D - plot_df$value3D) / plot_df$value3D * 100
  
  # Get correlation
  corr_df <- plot_df %>% 
    group_by(metric) %>% 
    summarise(
      corr = cor(value3D, value2D, method = "spearman", use = "complete.obs"), 
      .groups = "drop")
  
  # map corr ∈ [-1,1] → error ∈ [-100,400]
  corr_to_error <- function(c) scales::rescale(c, to = c(-100, 400), from = c(0, 1))
  
  # inverse transform for the secondary axis
  error_to_corr <- function(e) scales::rescale(e, to = c(0, 1), from = c(-100, 400))
  
  corr_df$y_trans <- corr_to_error(corr_df$corr)
  
  
  fig <- ggplot(plot_df, aes(x = metric, y = error)) +
    geom_boxplot(fill = "lightgray") +
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) +
    
    # median labels
    stat_summary(
      fun = function(z) median(z, na.rm = TRUE),
      geom = "text",
      aes(label = after_stat(sprintf("%.1f", y))),
      vjust = -0.5,
      size = 3.5,
      color = "black"
    ) +
    
    # ⭐ correlation stars
    geom_point(data = corr_df, 
               aes(x = metric, y = y_trans), 
              shape = 8, # star 
              size = 5, 
              color = "#0062c5") +
    
    labs(
      title = "Error Distribution by Metric, for a random slice and pair A/B, with Spearman Correlation",
      x = "Metric",
      y = "Error (%)"
    ) +

    scale_y_continuous(limits = c(-100, 400), 
                      sec.axis = sec_axis(~ error_to_corr(.), name = "Spearman Correlation")) +
    
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  
  return(fig)
}



# Plot 2D vs 3D for each metric and pair for averaged slice
plot_2D_vs_3D_by_metric_and_pair_A_B_for_averaged_slice_scatter_plot <- function(metric_df_list,
                                                                                 metrics) {
  plot_df <- data.frame()
  pval_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change slice column to numeric for safety
    metric_df$slice <- as.numeric(metric_df$slice)
    
    # Add pair column
    if (metric %in% c("EBSAC", "EBP_AUC")) {
      # For EBSAC and EBP_AUC, assume pair is the same as cell_types for consistency
      metric_df$pair <- gsub(',', '/', metric_df$cell_types)
    }
    else if (metric %in% c("ANE_AUC")) {
      # For ANE_AUC, assume pair is the same as target_cell_type for consistency (as target is of form A,B already)
      metric_df$pair <- gsub(',', '/', metric_df$target)
    }
    else {
      # Add reference-target column
      metric_df$pair <- paste(metric_df$reference, metric_df$target, sep = "/")
    }
    
    # Subset for A/B pair
    metric_df <- metric_df[metric_df$pair == "A/B", ]
    
    # Modify the metric_df so that there is a column of 3D values, and a column of 2D values
    colnames(metric_df)[colnames(metric_df) == metric] <- "metric_value"
    metric_df$sim_pair <- paste(metric_df$simulation, metric_df$pair, sep = "_")
    
    metric_df3D <- metric_df[metric_df[["slice"]] == 0, ]
    colnames(metric_df3D)[colnames(metric_df3D) == "metric_value"] <- "value3D"
    
    # Take the average of all slices for each simulation and pair to get 2D values
    metric_df2D <- metric_df[metric_df[["slice"]] != 0, ] %>%
      group_by(sim_pair) %>%
      summarise(value2D = mean(metric_value, na.rm = T))
    
    # Merge metric_df3D and metric_df2D
    metric_df <- merge(metric_df3D, metric_df2D, by = "sim_pair")
    
    metric_df$metric <- metric  
    
    # Add metric_df to plot_df
    plot_df <- rbind(plot_df, metric_df[ , c("value3D", "value2D", "pair", "metric")])
    
    # Add to pval_df
    p_values <- c()
    
    # Ignore when pair is invalid
    if (sum(is.finite(metric_df[["value3D"]])) == 0) {
      next
    } 
    
    wilcox_test  <- wilcox.test(metric_df[["value3D"]], 
                                metric_df[["value2D"]], 
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
    
    pval_df <- rbind(pval_df, 
                     data.frame(metric = metric,
                                p_value = p_values))
  }
  
  
  
  fig <- ggplot(plot_df, aes(x = value3D, y = value2D)) +
    geom_point(alpha = 0.25, color = "#0062c5", size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "#bb0036", linewidth = 1) + # Dotted line of the equation y = x
    labs(title = "Scatterplots showing 2D vs 3D, for each Metric, for averaged slices and pair A/B",
         x = "3D value",
         y = "2D value") +
    facet_wrap(~ interaction(metric), scales = "free", ncol = length(metrics)) +
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

# Plot error vs 3D for each metric and pair for averaged slice
plot_error_vs_3D_by_metric_and_pair_A_B_for_averaged_slice_scatter_plot <- function(metric_df_list,
                                                                                    metrics) {
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change slice column to numeric for safety
    metric_df$slice <- as.numeric(metric_df$slice)
    
    # Add pair column
    if (metric %in% c("EBSAC", "EBP_AUC")) {
      # For EBSAC and EBP_AUC, assume pair is the same as cell_types for consistency
      metric_df$pair <- gsub(',', '/', metric_df$cell_types)
    }
    else if (metric %in% c("ANE_AUC")) {
      # For ANE_AUC, assume pair is the same as target_cell_type for consistency (as target is of form A,B already)
      metric_df$pair <- gsub(',', '/', metric_df$target)
    }
    else {
      # Add reference-target column
      metric_df$pair <- paste(metric_df$reference, metric_df$target, sep = "/")
    }
    
    # Subset for A/B pair
    metric_df <- metric_df[metric_df$pair == "A/B", ]
    
    # Modify the metric_df so that there is a column of 3D values, and a column of 2D values
    colnames(metric_df)[colnames(metric_df) == metric] <- "metric_value"
    metric_df$sim_pair <- paste(metric_df$simulation, metric_df$pair, sep = "_")
    
    metric_df3D <- metric_df[metric_df[["slice"]] == 0, ]
    colnames(metric_df3D)[colnames(metric_df3D) == "metric_value"] <- "value3D"
    
    # Take the average of all slices for each simulation and pair to get 2D values
    metric_df2D <- metric_df[metric_df[["slice"]] != 0, ] %>%
      group_by(sim_pair) %>%
      summarise(value2D = mean(metric_value, na.rm = T))
    
    # Merge metric_df3D and metric_df2D
    metric_df <- merge(metric_df3D, metric_df2D, by = "sim_pair")
    
    metric_df$metric <- metric  
    
    # Add metric_df to plot_df
    plot_df <- rbind(plot_df, metric_df[ , c("value3D", "value2D", "pair", "metric")])
    
  }
  
  # Get error column
  plot_df$error <- (plot_df$value2D - plot_df$value3D) / plot_df$value3D * 100
  
  fig <- ggplot(plot_df, aes(x = value3D, y = error)) +
    geom_point(alpha = 0.25, color = "#0062c5", size = 1) +
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    labs(title = "Scatterplots showing Error vs 3D, for each Metric, for averaged slices and pair A/B",
         x = "3D value",
         y = "Error (%)") +
    facet_wrap(~ interaction(metric), scales = "free", ncol = length(metrics)) +  
    scale_y_continuous(limits = c(-100, 500)) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(size = 6),  # make x-axis text smaller
      axis.text.y = element_text(size = 6)   # make y-axis text smaller
    )
  
  return(fig)
}

# Plot error vs 3D for each metric and pair for averaged slice
plot_error_vs_3D_by_metric_and_pair_A_B_for_averaged_slice_box_plot <- function(metric_df_list,
                                                                                metrics) {
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change slice column to numeric for safety
    metric_df$slice <- as.numeric(metric_df$slice)
    
    # Add pair column
    if (metric %in% c("EBSAC", "EBP_AUC")) {
      # For EBSAC and EBP_AUC, assume pair is the same as cell_types for consistency
      metric_df$pair <- gsub(',', '/', metric_df$cell_types)
    }
    else if (metric %in% c("ANE_AUC")) {
      # For ANE_AUC, assume pair is the same as target_cell_type for consistency (as target is of form A,B already)
      metric_df$pair <- gsub(',', '/', metric_df$target)
    }
    else {
      # Add reference-target column
      metric_df$pair <- paste(metric_df$reference, metric_df$target, sep = "/")
    }
    
    # Subset for A/B pair
    metric_df <- metric_df[metric_df$pair == "A/B", ]
    
    # Modify the metric_df so that there is a column of 3D values, and a column of 2D values
    colnames(metric_df)[colnames(metric_df) == metric] <- "metric_value"
    metric_df$sim_pair <- paste(metric_df$simulation, metric_df$pair, sep = "_")
    
    metric_df3D <- metric_df[metric_df[["slice"]] == 0, ]
    colnames(metric_df3D)[colnames(metric_df3D) == "metric_value"] <- "value3D"
    
    # Take the average of all slices for each simulation and pair to get 2D values
    metric_df2D <- metric_df[metric_df[["slice"]] != 0, ] %>%
      group_by(sim_pair) %>%
      summarise(value2D = mean(metric_value, na.rm = T))
    
    # Merge metric_df3D and metric_df2D
    metric_df <- merge(metric_df3D, metric_df2D, by = "sim_pair")
    
    metric_df$metric <- metric  
    
    # Add metric_df to plot_df
    plot_df <- rbind(plot_df, metric_df[ , c("value3D", "value2D", "pair", "metric")])
    
  }
  
  # Get error column
  plot_df$error <- (plot_df$value2D - plot_df$value3D) / plot_df$value3D * 100
  
  # Get correlation
  corr_df <- plot_df %>% 
    group_by(metric) %>% 
    summarise(
      corr = cor(value3D, value2D, method = "spearman", use = "complete.obs"), 
      .groups = "drop")
  
  # map corr ∈ [-1,1] → error ∈ [-100,400]
  corr_to_error <- function(c) scales::rescale(c, to = c(-100, 400), from = c(0, 1))
  
  # inverse transform for the secondary axis
  error_to_corr <- function(e) scales::rescale(e, to = c(0, 1), from = c(-100, 400))
  
  corr_df$y_trans <- corr_to_error(corr_df$corr)
  
  fig <- ggplot(plot_df, aes(x = metric, y = error)) +
    geom_boxplot(fill = "lightgray") +
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    
    # Median labels
    stat_summary(
      fun = function(z) median(z, na.rm = TRUE),
      geom = "text",
      aes(label = after_stat(sprintf("%.1f", y))),
      vjust = -0.5,
      size = 3.5,
      color = "black"
    ) +
    
    # ⭐ correlation stars
    geom_point(data = corr_df, 
               aes(x = metric, y = y_trans), 
               shape = 8, # star 
               size = 5, 
               color = "#0062c5") +
    
    labs(title = "Error Distribution by Metric, for averaged slices and pair A/B",
         x = "Metric",
         y = "Error (%)") +
    
    scale_y_continuous(limits = c(-100, 400), 
                       sec.axis = sec_axis(~ error_to_corr(.), name = "Spearman Correlation")) +
    
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  return(fig)
}



# Plot 2D vs 3D for each metric and pair for 3 slices
plot_2D_vs_3D_by_metric_and_pair_A_B_for_three_slices_scatter_plot <- function(metric_df_list,
                                                                               metrics) {
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change slice column to numeric for safety
    metric_df$slice <- as.numeric(metric_df$slice)
    
    # Add pair column
    if (metric %in% c("EBSAC", "EBP_AUC")) {
      # For EBSAC and EBP_AUC, assume pair is the same as cell_types for consistency
      metric_df$pair <- gsub(',', '/', metric_df$cell_types)
    }
    else if (metric %in% c("ANE_AUC")) {
      # For ANE_AUC, assume pair is the same as target_cell_type for consistency (as target is of form A,B already)
      metric_df$pair <- gsub(',', '/', metric_df$target)
    }
    else {
      # Add reference-target column
      metric_df$pair <- paste(metric_df$reference, metric_df$target, sep = "/")
    }
    
    # Subset for A/B pair
    metric_df <- metric_df[metric_df$pair == "A/B", ]
    
    # Modify the metric_df so that there is a column of 3D values, and a column of 2D values
    colnames(metric_df)[colnames(metric_df) == metric] <- "metric_value"
    metric_df$sim_pair <- paste(metric_df$simulation, metric_df$pair, sep = "_")
    
    metric_df3D <- metric_df[metric_df[["slice"]] == 0, ]
    colnames(metric_df3D)[colnames(metric_df3D) == "metric_value"] <- "value3D"
    
    # Take 3 slices for each simulation and pair to get 2D values
    # Slice indices 7, 10, 13 correspond with 145, 175, 205 z coord
    metric_df2D <- metric_df[metric_df[["slice"]] %in% c(7, 10, 13), ] 
    
    # Merge metric_df3D and metric_df2D
    metric_df <- metric_df2D %>% left_join(metric_df3D[, c("sim_pair", "value3D")], by = "sim_pair")
    colnames(metric_df)[colnames(metric_df) == "metric_value"] <- "value2D"
    
    metric_df$metric <- metric  
    
    # Add metric_df to plot_df
    plot_df <- rbind(plot_df, metric_df[ , c("value3D", "value2D", "pair", "metric", "slice")])
    
  }
  
  # Change slice column back to character for plotting
  plot_df$slice <- as.character(metric_df$slice)
  
  plot_df$slice <- factor(plot_df$slice, levels = c('7', '10', '13'))
  
  fig <- ggplot(plot_df, aes(x = value3D, y = value2D, color = slice)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "#bb0036", linewidth = 1) + # Dotted line of the equation y = x
    labs(title = "Scatterplots showing 2D vs 3D, for each Metric, 3 slices and pair A/B",
         x = "3D value",
         y = "2D value") +
    facet_wrap(~ interaction(metric), scales = "free", ncol = length(metrics)) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(size = 4),  # make x-axis text smaller
      axis.text.y = element_text(size = 4)   # make y-axis text smaller
    ) +
    scale_color_manual(
      values = c(
        "7" = "#9437a8",
        "10" = "#007128",
        "13" = "#b8db50"
        
      )
    )
  
  return(fig)
}

# Plot error vs 3D for each metric and pair for 3 slices
plot_error_vs_3D_by_metric_and_pair_A_B_for_three_slices_scatter_plot <- function(metric_df_list,
                                                                                  metrics) {
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change slice column to numeric for safety
    metric_df$slice <- as.numeric(metric_df$slice)
    
    # Add pair column
    if (metric %in% c("EBSAC", "EBP_AUC")) {
      # For EBSAC and EBP_AUC, assume pair is the same as cell_types for consistency
      metric_df$pair <- gsub(',', '/', metric_df$cell_types)
    }
    else if (metric %in% c("ANE_AUC")) {
      # For ANE_AUC, assume pair is the same as target_cell_type for consistency (as target is of form A,B already)
      metric_df$pair <- gsub(',', '/', metric_df$target)
    }
    else {
      # Add reference-target column
      metric_df$pair <- paste(metric_df$reference, metric_df$target, sep = "/")
    }
    
    # Subset for A/B pair
    metric_df <- metric_df[metric_df$pair == "A/B", ]
    
    # Modify the metric_df so that there is a column of 3D values, and a column of 2D values
    colnames(metric_df)[colnames(metric_df) == metric] <- "metric_value"
    metric_df$sim_pair <- paste(metric_df$simulation, metric_df$pair, sep = "_")
    
    metric_df3D <- metric_df[metric_df[["slice"]] == 0, ]
    colnames(metric_df3D)[colnames(metric_df3D) == "metric_value"] <- "value3D"
    
    # Take 3 slices for each simulation and pair to get 2D values
    # Slice indices 7, 10, 13 correspond with 145, 175, 205 z coord
    metric_df2D <- metric_df[metric_df[["slice"]] %in% c(7, 10, 13), ] 
    
    # Merge metric_df3D and metric_df2D
    metric_df <- metric_df2D %>% left_join(metric_df3D[, c("sim_pair", "value3D")], by = "sim_pair")
    colnames(metric_df)[colnames(metric_df) == "metric_value"] <- "value2D"
    
    metric_df$metric <- metric  
    
    # Add metric_df to plot_df
    plot_df <- rbind(plot_df, metric_df[ , c("value3D", "value2D", "pair", "metric", "slice")])
    
  }
  
  # Change slice column back to character for plotting
  plot_df$slice <- as.character(metric_df$slice)
  
  plot_df$slice <- factor(plot_df$slice, levels = c('7', '10', '13'))
  
  # Get error column
  plot_df$error <- (plot_df$value2D - plot_df$value3D) / plot_df$value3D * 100
  
  fig <- ggplot(plot_df, aes(x = value3D, y = error, color = slice)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    labs(title = "Scatterplots showing Error vs 3D, for each Metric, 3 slices and pair A/B",
         x = "3D value",
         y = "Error (%)") +
    facet_wrap(~ interaction(metric), scales = "free", ncol = length(metrics)) +  
    scale_y_continuous(limits = c(-100, 500)) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(size = 6),  # make x-axis text smaller
      axis.text.y = element_text(size = 6)   # make y-axis text smaller
    ) +
    scale_color_manual(
      values = c(
        "7" = "#9437a8",
        "10" = "#007128",
        "13" = "#b8db50"
        
      )
    )
  
  return(fig)
}



# Plot 2D vs 3D for each metric and pair for a random slice, annotating for tissue structure
plot_2D_vs_3D_by_metric_and_pair_A_B_for_random_slice_showing_structure_scatter_plot <- function(metric_df_list,
                                                                                                 metrics,
                                                                                                 parameters_df) {
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change slice column to numeric for safety
    metric_df$slice <- as.numeric(metric_df$slice)
    
    # Add pair column
    if (metric %in% c("EBSAC", "EBP_AUC")) {
      # For EBSAC and EBP_AUC, assume pair is the same as cell_types for consistency
      metric_df$pair <- gsub(',', '/', metric_df$cell_types)
    }
    else if (metric %in% c("ANE_AUC")) {
      # For ANE_AUC, assume pair is the same as target_cell_type for consistency (as target is of form A,B already)
      metric_df$pair <- gsub(',', '/', metric_df$target)
    }
    else {
      # Add reference-target column
      metric_df$pair <- paste(metric_df$reference, metric_df$target, sep = "/")
    }
    
    # Subset for A/B pair
    metric_df <- metric_df[metric_df$pair == "A/B", ]
    
    # Modify the metric_df so that there is a column of 3D values, and a column of 2D values
    # Choose a random slice for 2D values
    metric_values2D <- metric_df[[metric]][metric_df[["slice"]] == sample(unique(metric_df[["slice"]]), 1)]
    # Slice == 0 is the 3D value
    metric_df <- metric_df[metric_df[["slice"]] == 0, ]
    metric_df[["value2D"]] <- metric_values2D
    colnames(metric_df)[colnames(metric_df) == metric] <- "value3D"
    
    # Merge metric_df and parameters_df
    parameters_df$simulation <- seq(nrow(parameters_df))
    metric_df <- metric_df %>% left_join(parameters_df[, c("simulation", "arrangement", "shape")], by = "simulation")
    
    # Get structure column
    metric_df$structure <- paste(metric_df$arrangement, metric_df$shape, sep = "_")
    
    metric_df$metric <- metric  
    
    # Add metric_df to plot_df
    plot_df <- rbind(plot_df, metric_df[ , c("value3D", "value2D", "pair", "metric", "structure")])
    
  }
  
  
  
  fig <- ggplot(plot_df, aes(x = value3D, y = value2D, color = structure)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "#bb0036", linewidth = 1) + # Dotted line of the equation y = x
    labs(title = "Scatterplots showing 2D vs 3D, for each Metric and pair A/B, showing structure",
         x = "3D value",
         y = "2D value") +
    facet_wrap(~ interaction(metric), scales = "free", ncol = length(metrics)) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(size = 4),  # make x-axis text smaller
      axis.text.y = element_text(size = 4)   # make y-axis text smaller
    ) +
    scale_color_manual(
      values = c(
        "mixed_ellipsoid" = "#007128",
        "mixed_network" = "#b8db50",
        "ringed_ellipsoid" = "#9437a8",
        "ringed_network" = "#d99dff",
        "separated_ellipsoid" = "#770026",
        "separated_network" = "#bb0036"
      )
    )
  
  return(fig)
}

# Plot error vs 3D for each metric and pair for a random slice, annotating for tissue structure
plot_error_vs_3D_by_metric_and_pair_A_B_for_random_slice_showing_structure_scatter_plot <- function(metric_df_list,
                                                                                                    metrics,
                                                                                                    parameters_df) {
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change slice column to numeric for safety
    metric_df$slice <- as.numeric(metric_df$slice)
    
    # Add pair column
    if (metric %in% c("EBSAC", "EBP_AUC")) {
      # For EBSAC and EBP_AUC, assume pair is the same as cell_types for consistency
      metric_df$pair <- gsub(',', '/', metric_df$cell_types)
    }
    else if (metric %in% c("ANE_AUC")) {
      # For ANE_AUC, assume pair is the same as target_cell_type for consistency (as target is of form A,B already)
      metric_df$pair <- gsub(',', '/', metric_df$target)
    }
    else {
      # Add reference-target column
      metric_df$pair <- paste(metric_df$reference, metric_df$target, sep = "/")
    }
    
    # Subset for A/B pair
    metric_df <- metric_df[metric_df$pair == "A/B", ]
    
    # Modify the metric_df so that there is a column of 3D values, and a column of 2D values
    # Choose a random slice for 2D values
    metric_values2D <- metric_df[[metric]][metric_df[["slice"]] == sample(unique(metric_df[["slice"]]), 1)]
    # Slice == 0 is the 3D value
    metric_df <- metric_df[metric_df[["slice"]] == 0, ]
    metric_df[["value2D"]] <- metric_values2D
    colnames(metric_df)[colnames(metric_df) == metric] <- "value3D"
    
    # Merge metric_df and parameters_df
    parameters_df$simulation <- seq(nrow(parameters_df))
    metric_df <- metric_df %>% left_join(parameters_df[, c("simulation", "arrangement", "shape")], by = "simulation")
    
    # Get structure column
    metric_df$structure <- paste(metric_df$arrangement, metric_df$shape, sep = "_")
    
    metric_df$metric <- metric  
    
    # Add metric_df to plot_df
    plot_df <- rbind(plot_df, metric_df[ , c("value3D", "value2D", "pair", "metric", "structure")])
    
  }
  
  # Get error column
  plot_df$error <- (plot_df$value2D - plot_df$value3D) / plot_df$value3D * 100
  
  fig <- ggplot(plot_df, aes(x = value3D, y = value2D, color = structure)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    labs(title = "Scatterplots showing Error vs 3D, for each Metric and pair A/B, showing structure",
         x = "3D value",
         y = "Error (%)") +
    facet_wrap(~ interaction(metric), scales = "free", ncol = length(metrics)) +  
    scale_y_continuous(limits = c(-100, 500)) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(size = 6),  # make x-axis text smaller
      axis.text.y = element_text(size = 6)   # make y-axis text smaller
    ) +
    scale_color_manual(
      values = c(
        "mixed_ellipsoid" = "#007128",
        "mixed_network" = "#b8db50",
        "ringed_ellipsoid" = "#9437a8",
        "ringed_network" = "#d99dff",
        "separated_ellipsoid" = "#770026",
        "separated_network" = "#bb0036"
      )
    )
  
  return(fig)
}

# Plot 2D vs 3D correlation vs tissue structure for each metric and pair for a random slice
plot_2D_vs_3D_correlation_vs_tissue_structure_by_metric_and_pair_A_B_for_random_slice_bar_plot <- function(metric_df_list,
                                                                                                           metrics,
                                                                                                           parameters_df) {
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change slice column to numeric for safety
    metric_df$slice <- as.numeric(metric_df$slice)
    
    # Add pair column
    if (metric %in% c("EBSAC", "EBP_AUC")) {
      # For EBSAC and EBP_AUC, assume pair is the same as cell_types for consistency
      metric_df$pair <- gsub(',', '/', metric_df$cell_types)
    }
    else if (metric %in% c("ANE_AUC")) {
      # For ANE_AUC, assume pair is the same as target_cell_type for consistency (as target is of form A,B already)
      metric_df$pair <- gsub(',', '/', metric_df$target)
    }
    else {
      # Add reference-target column
      metric_df$pair <- paste(metric_df$reference, metric_df$target, sep = "/")
    }
    
    # Subset for A/B pair
    metric_df <- metric_df[metric_df$pair == "A/B", ]
    
    # Modify the metric_df so that there is a column of 3D values, and a column of 2D values
    # Choose a random slice for 2D values
    metric_values2D <- metric_df[[metric]][metric_df[["slice"]] == sample(unique(metric_df[["slice"]]), 1)]
    # Slice == 0 is the 3D value
    metric_df <- metric_df[metric_df[["slice"]] == 0, ]
    metric_df[["value2D"]] <- metric_values2D
    colnames(metric_df)[colnames(metric_df) == metric] <- "value3D"
    
    # Merge metric_df and parameters_df
    parameters_df$simulation <- seq(nrow(parameters_df))
    metric_df <- metric_df %>% left_join(parameters_df[, c("simulation", "arrangement", "shape")], by = "simulation")
    
    # Get structure column
    metric_df$structure <- paste(metric_df$arrangement, metric_df$shape, sep = "_")
    
    metric_df$metric <- metric  
    
    # Add metric_df to plot_df
    plot_df <- rbind(plot_df, metric_df[ , c("value3D", "value2D", "pair", "metric", "structure")])
    
  }
  
  # Compute spearman correlation
  corr_df <- plot_df %>%
    group_by(pair, metric, structure) %>%
    summarise(
      spearman_corr = cor(value3D, value2D, method = "spearman", use = "complete.obs"),
      .groups = "drop"
    )
  
  
  
  fig <- ggplot(corr_df, aes(x = structure, y = spearman_corr, fill = structure)) + 
    geom_col(alpha = 0.8) + 
    labs(title = "Bar plots showing Spearman Correlation vs Structure, for each Metric, a random slice and pair A/B", x = "Structure", y = "Spearman Correlation") + 
    facet_wrap(~ interaction(metric), scales = "free", ncol = length(metrics)) +
    theme_minimal() + 
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
          axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_manual(
      values = c(
        "mixed_ellipsoid" = "#007128",
        "mixed_network" = "#b8db50",
        "ringed_ellipsoid" = "#9437a8",
        "ringed_network" = "#d99dff",
        "separated_ellipsoid" = "#770026",
        "separated_network" = "#bb0036"
      )
    )
  
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

metrics <- c("AMD", "ANC_AUC", "ACIN_AUC", "ANE_AUC", "MS_AUC", "NMS_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "CK_AUC", "CL_AUC", "CG_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")


fig_2D_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot <- plot_2D_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot(metric_df_list,
                                                                                                                                        metrics)

fig_error_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot <- plot_error_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot(metric_df_list,
                                                                                                                                              metrics)

fig_error_vs_3D_by_metric_and_pair_A_B_for_random_slice_box_plot <- plot_error_vs_3D_by_metric_and_pair_A_B_for_random_slice_box_plot(metric_df_list,
                                                                                                                                      metrics)



fig_2D_vs_3D_by_metric_and_pair_A_B_for_averaged_slice_scatter_plot <- plot_2D_vs_3D_by_metric_and_pair_A_B_for_averaged_slice_scatter_plot(metric_df_list,
                                                                                                                                            metrics)

fig_error_vs_3D_by_metric_and_pair_A_B_for_averaged_slice_scatter_plot <- plot_error_vs_3D_by_metric_and_pair_A_B_for_averaged_slice_scatter_plot(metric_df_list,
                                                                                                                                                  metrics)

fig_error_vs_3D_by_metric_and_pair_A_B_for_averaged_slice_box_plot <- plot_error_vs_3D_by_metric_and_pair_A_B_for_averaged_slice_box_plot(metric_df_list,
                                                                                                                                          metrics)



fig_2D_vs_3D_by_metric_and_pair_A_B_for_three_slices_scatter_plot <- plot_2D_vs_3D_by_metric_and_pair_A_B_for_three_slices_scatter_plot(metric_df_list,
                                                                                                                                        metrics)

fig_error_vs_3D_by_metric_and_pair_A_B_for_three_slices_scatter_plot <- plot_error_vs_3D_by_metric_and_pair_A_B_for_three_slices_scatter_plot(metric_df_list,
                                                                                                                                              metrics)



fig_2D_vs_3D_by_metric_and_pair_A_B_for_random_slice_showing_structure_scatter_plot <- plot_2D_vs_3D_by_metric_and_pair_A_B_for_random_slice_showing_structure_scatter_plot(metric_df_list,
                                                                                                                                                                            metrics,
                                                                                                                                                                            parameters_df)

fig_error_vs_3D_by_metric_and_pair_A_B_for_random_slice_showing_structure_scatter_plot <- plot_error_vs_3D_by_metric_and_pair_A_B_for_random_slice_showing_structure_scatter_plot(metric_df_list,
                                                                                                                                                                                  metrics,
                                                                                                                                                                                  parameters_df)

fig_2D_vs_3D_correlation_vs_tissue_structure_by_metric_and_pair_A_B_for_random_slice_bar_plot <- plot_2D_vs_3D_correlation_vs_tissue_structure_by_metric_and_pair_A_B_for_random_slice_bar_plot(metric_df_list,
                                                                                                                                                                                                metrics,
                                                                                                                                                                                                parameters_df)

setwd("~/R/plots/S1")
pdf("random_slice_pair_A_B.pdf", width = 28, height = 2)

print(fig_2D_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot)
print(fig_error_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot)
print(fig_error_vs_3D_by_metric_and_pair_A_B_for_random_slice_box_plot)

dev.off()


setwd("~/R/plots/S1")
pdf("averaged_slices_pair_A_B.pdf", width = 24, height = 3)

print(fig_2D_vs_3D_by_metric_and_pair_A_B_for_averaged_slice_scatter_plot)
print(fig_error_vs_3D_by_metric_and_pair_A_B_for_averaged_slice_scatter_plot)
print(fig_error_vs_3D_by_metric_and_pair_A_B_for_averaged_slice_box_plot)

dev.off()


setwd("~/R/plots/S1")
pdf("three_slices_pair_A_B.pdf", width = 24, height = 3)

print(fig_2D_vs_3D_by_metric_and_pair_A_B_for_three_slices_scatter_plot)
print(fig_error_vs_3D_by_metric_and_pair_A_B_for_three_slices_scatter_plot)

dev.off()


setwd("~/R/plots/S1")
pdf("random_slice_showing_structure_pair_A_B.pdf", width = 24, height = 3)

print(fig_2D_vs_3D_by_metric_and_pair_A_B_for_random_slice_showing_structure_scatter_plot)
print(fig_error_vs_3D_by_metric_and_pair_A_B_for_random_slice_showing_structure_scatter_plot)
print(fig_2D_vs_3D_correlation_vs_tissue_structure_by_metric_and_pair_A_B_for_random_slice_bar_plot)

dev.off()


