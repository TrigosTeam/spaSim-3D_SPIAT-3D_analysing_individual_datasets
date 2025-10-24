library(cowplot)
library(ggplot2)
library(S4Vectors)
library(stringr)
library(dplyr)
library(DescTools)

### Utility function to get metric cell types -----
get_metric_cell_types <- function(metric) {
  # Get metric_cell_types
  if (metric %in% c("AMD", "ACIN", "CKR", "CLR", "COO", "CGR", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC")) {
    metric_cell_types <- data.frame(ref = c("A", "A", "B", "B"), tar = c("A", "B", "A", "B"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("MS", "NMS", "MS_AUC", "NMS_AUC")) {
    metric_cell_types <- data.frame(ref = c("A", "B"), tar = c("B", "A"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("ACINP", "ACINP_AUC")) {
    metric_cell_types <- data.frame(ref = c("A", "B"), tar = c("A", "A"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("AE", "AE_AUC")) {
    metric_cell_types <- data.frame(ref = c("A", "B"), tar = c("A,B", "A,B"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("PBSAC", "PBP", "PBP_AUC")) {
    metric_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("EBSAC", "EBP", "EBP_AUC")) {
    metric_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))
  }
  else {
    stop("metric not found")
  }
  return(metric_cell_types)
}


### Utility function to subset metric_df -----
subset_metric_df <- function(metric,
                             metric_df,
                             metric_cell_types,
                             index) {
  if (metric %in% c("AMD", "ACIN", "CKR", "CLR", "COO", "CGR", "MS", "NMS", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "MS_AUC", "NMS_AUC", "PBSAC", "PBP", "PBP_AUC")) {
    metric_df <- metric_df[metric_df$reference == metric_cell_types[index, "ref"] & metric_df$target == metric_cell_types[index, "tar"], ] 
  }
  else if (metric %in% c("ACINP", "AE", "ACINP_AUC", "AE_AUC")) {
    metric_df <- metric_df[metric_df$reference == metric_cell_types[index, "ref"], ] 
  }
  else if (metric %in% c("EBSAC", "EBP", "EBP_AUC")) {
    metric_df <- metric_df[metric_df$cell_types == metric_cell_types[index, "cell_types"], ]
  }
  else {
    stop("metric not found")
  }
  return(metric_df)
}

### Utility function to duplicate df by *rows* -----
duplicate_df <- function(df, n_times) {
  df <- df %>%
    mutate(row_num = row_number())
  df <- do.call(bind_rows, replicate(n_times, df, simplify = FALSE)) %>%
    arrange(row_num)
  df$row_num <- NULL
  
  return(df)
}
### Utility function to get title --------
get_metric_cell_types_title <- function(metric, metric_cell_types, index) {
  
  if (metric %in% c("AMD", "ACIN", "ACINP", "AE", "CKR", "CLR", "COO", "CGR", "MS", "NMS", "ACIN_AUC", "ACINP_AUC", "AE_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "MS_AUC", "NMS_AUC", "PBSAC", "PBP", "PBP_AUC")) {
    title <- ggdraw() +
      draw_label(paste("Reference:", metric_cell_types[index, "ref"], "Target:", metric_cell_types[index, "tar"]),
                 fontface = 'bold')
  }
  else if (metric %in% c("EBSAC", "EBP", "EBP_AUC")) {
    title <- ggdraw() + 
      draw_label(paste("Cell types of interest:", metric_cell_types[index, "cell_types"]), 
                 fontface='bold')
  }
  else {
    stop("metric not found")
  }
  return(title)
}
### Function for non-gradient output ------------------------------------------
plot_non_gradient_metric <- function(spes_table, 
                                     metric, 
                                     metric_df, 
                                     arrangement, 
                                     plots_metadata) {
  
  ### Modify plots_metadata
  # Change plots_metadata arrangement to inputted arrangement
  plots_metadata$arrangement$x_aes <- arrangement
  
  # Change plots_metadata y_aes to inputted metric
  for (i in seq_along(plots_metadata)) {
    # Modify the y_aes element
    plots_metadata[[i]]$y_aes <- metric
  }
  
  # Get metric_cell_types
  metric_cell_types <- get_metric_cell_types(metric)
  
  # Define plotting function
  formatCustomSci <- function(x) {
    x_sci <- str_split_fixed(formatC(x, format = "e"), "e", 2)
    alpha <- round(as.numeric(x_sci[ , 1]), 1)
    power <- as.integer(x_sci[ , 2])
    paste(alpha, power, sep = "e")
  }
  
  create_plot <- function(data, x_aes, y_aes, title = "") {
    
    size <- 0.5
    
    data <- data[data$variable_parameter == x_aes, ]
    
    if (sum(!is.na(data[[y_aes]])) >= 2 && any(data[[y_aes]] != 0)) {
      correlation_test <- cor.test(data[[x_aes]], data[[y_aes]], method = "spearman")
      correlation <- correlation_test$estimate
      p_value <- correlation_test$p.value
      if (p_value == 0) p_value <- 2.2e-308
      if (0 < p_value && p_value < 1e-3)  {
        p_value <- formatCustomSci(p_value)
      }
      else {
        p_value  <- round(p_value, 3)
      }
      title <- paste("r:", round(correlation, 3), ", p:", p_value)
    }
    
    plot <- ggplot(data, aes_string(x = x_aes, y = y_aes)) +
      labs(x = x_aes, y = y_aes) +
      ggtitle(title) +
      theme_bw() +
      theme(plot.title = element_text(size = 9))
    
    # Use scientific notation for ellipsoid volume
    if (x_aes == "E_volume") {
      plot <- plot + 
        geom_point(size = size) +
        scale_x_continuous(labels = formatCustomSci)
    }
    else if (typeof(data[[x_aes]]) == "double") {
      breaks <- pretty(c(min(data[[x_aes]]), max(data[[x_aes]])), n = 2)
      plot <- plot + 
        geom_point(size = size) + 
        scale_x_continuous(breaks = breaks)
    }
    # Factored character is an integer
    else if (typeof(data[[x_aes]]) %in% c("integer", "character")) {
      plot <- plot + 
        geom_violin()
    }
    return(plot)
  }
  
  # Put plots into an organised list
  plots_list <- list()
  
  for (i in seq(nrow(metric_cell_types))) {
    
    # Subset metric_df for chosen pair/cell types
    metric_df_subset <- subset_metric_df(metric, metric_df, metric_cell_types, i)
    
    # Combine spes_table and metric_df
    plot_df <- cbind(spes_table, metric_df_subset)
    
    # Factor
    if (!is.null(plot_df$shape)) plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    if (!is.null(plot_df$slice)) plot_df$slice <- as.character(plot_df$slice)
    
    # Generate plots based on plots_metadata, use final column of metric_cell_types
    plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]] <- lapply(plots_metadata, function(plot_def) {
      x_aes <- plot_def$x_aes
      y_aes <- plot_def$y_aes
      title <- " "
      plot <- create_plot(data = plot_df, x_aes = x_aes, y_aes = y_aes, title = title)
      return(plot)
    })
  }
  
  # Combine the plots together using metric_cell_types
  combined_plots_list <- list()
  for (i in seq(nrow(metric_cell_types))) {
    
    # Get final column
    cells <- metric_cell_types[i, ncol(metric_cell_types)]
    
    plots <- plot_grid(plotlist = plots_list[[cells]], nrow = 1, ncol = length(plots_list[[cells]]))
    
    title <- get_metric_cell_types_title(metric, metric_cell_types, i)
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    combined_plots_list[[cells]] <- fig
  }
  
  # Combine the combined plots into one big plot
  non_gradient_metric_plot <- plot_grid(plotlist = combined_plots_list,
                                        nrow = length(combined_plots_list), 
                                        ncol = 1)
  
  return(non_gradient_metric_plot)
}

### Function for gradient output ----------------------------------------------
plot_gradient_metric <- function(spes_table, 
                                 metric, 
                                 metric_df, 
                                 arrangement, 
                                 gradient_type,
                                 plots_metadata) {
  ### Modify plots_metadata
  # Change plots_metadata arrangement to inputted arrangement
  plots_metadata$arrangement$color_aes <- arrangement
  
  # Change plots_metadata x_aes to gradient_type and y_aes to inputted metric
  for (i in seq_along(plots_metadata)) {
    plots_metadata[[i]]$x_aes <- gradient_type
    plots_metadata[[i]]$y_aes <- metric
  }
  
  # Get gradient/threshold values
  if (gradient_type == "radius") {
    gradient <- seq(20, 100, 10)
    gradient_colnames <- paste("r", gradient, sep = "")    
  }
  else if (gradient_type == "threshold") {
    gradient <- seq(0.01, 1, 0.01)
    gradient_colnames <- paste("t", gradient, sep = "")
  }
  
  # Get metric_cell_types
  metric_cell_types <- get_metric_cell_types(metric)
  
  # Define plotting function
  formatCustomSci <- function(x) {
    x_sci <- str_split_fixed(formatC(x, format = "e"), "e", 2)
    alpha <- as.numeric(x_sci[ , 1])
    power <- as.integer(x_sci[ , 2])
    paste(alpha, power, sep = "e")
  }
  create_plot <- function(data, x_aes, y_aes, color_aes, group_aes, title = "") {
    
    data <- data[data$variable_parameter == color_aes, ]
    
    breaks <- pretty(c(min(data[[color_aes]]), max(data[[color_aes]])), n = 3)
    
    plot <- ggplot(data, aes_string(x = x_aes, y = y_aes, group = group_aes, color = color_aes)) +
      labs(title = title, x = x_aes, y = y_aes) +
      theme_bw() +
      geom_line()
    
    # Use scientific notation for ellipsoid volume
    if (color_aes == "E_volume") {
      plot <- plot + scale_color_continuous(breaks = breaks, labels = formatCustomSci(breaks))
    }
    else {
      plot <- plot + scale_color_continuous(breaks = breaks)
    }
    
    return(plot)
  }
  
  # Put plots into an organised list
  plots_list <- list()
  
  for (i in seq(nrow(metric_cell_types))) {
    
    # Subset metric_df for chosen pair/cell types
    metric_df_subset <- subset_metric_df(metric, metric_df, metric_cell_types, i)
    
    # Combine spes_table and metric_df
    plot_df <- cbind(spes_table, metric_df_subset)
    
    # Melt
    plot_df <- reshape2::melt(plot_df, , gradient_colnames)
    
    # Change last 2 column names
    colnames(plot_df)[c(ncol(plot_df) - 1, ncol(plot_df))] <- c(gradient_type, metric)
    
    # Extract radius value from radius strings (r1 -> 1, r2 -> 2...)
    plot_df[[gradient_type]] <- unfactor(plot_df[[gradient_type]])
    plot_df[[gradient_type]] <- as.numeric(substr(plot_df[[gradient_type]], 2, nchar(plot_df[[gradient_type]])))
    
    # Factor
    if (!is.null(plot_df$shape)) plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    
    if (!is.null(plot_df$slice)) {
      plot_df$slice <- as.character(plot_df$slice)
      plot_df$key <- paste(plot_df$spe, plot_df$slice, sep = "_")
      group_aes = "key"
    }
    else {
      group_aes = "spe"
    }
    
    # Generate plots based on plots_metadata, use final column of metric_cell_types
    plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]] <- lapply(plots_metadata, function(plot_def) {
      x_aes <- plot_def$x_aes
      y_aes <- plot_def$y_aes
      color_aes <- plot_def$color_aes
      title <- plot_def$title
      plot <- create_plot(data = plot_df, x_aes = x_aes, y_aes = y_aes, group_aes = group_aes, color_aes = color_aes, title = title)
      return(plot)
    })
  }
  
  # Extract legends from first set of plots
  legends_list <- lapply(plots_list[[1]], function(plot) {
    plot_legend <- get_legend(plot + theme(legend.direction = "horizontal"))
    return(plot_legend)
  })
  legends <- plot_grid(plotlist = legends_list, nrow = 1)
  
  # Combine the plots together using metric_cell_types
  combined_plots_list <- list()
  for (i in seq(nrow(metric_cell_types))) {
    
    # Remove legend from base plots
    for (j in seq(length(plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]]))) {
      plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]][[j]] <- 
        plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]][[j]] + theme(legend.position = "none")
    }
    
    # Getting current set of cell types from metric_cell_types
    cells <- metric_cell_types[i, ncol(metric_cell_types)]
    
    plots <- plot_grid(plotlist = plots_list[[cells]], nrow = 1, ncol = length(plots_list[[cells]]))
    
    title <- get_metric_cell_types_title(metric, metric_cell_types, i)
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    combined_plots_list[[cells]] <- fig
  }
  
  # Combine the combined plots into one big plot
  gradient_metric_plot <- plot_grid(plotlist = combined_plots_list,
                                    nrow = length(combined_plots_list), ncol = 1)
  
  # Add legends
  gradient_metric_with_legends_plot <- plot_grid(gradient_metric_plot, legends,
                                                 nrow = 2, ncol = 1,
                                                 rel_heights = c(1, 0.1))
  
  return(gradient_metric_with_legends_plot)
}


### Function to compare 3D vs 2D one slice ------------------------------------------
plot_3D_vs_2D_metric_one_slice <- function(spes_table, 
                                           metric, 
                                           metric_df3D,
                                           metric_df2D,
                                           arrangement, 
                                           plots_metadata) {
  
  ### Modify plots_metadata
  # Change plots_metadata arrangement to inputted arrangement
  plots_metadata$arrangement$color_aes <- arrangement
  
  # Get metric_cell_types
  metric_cell_types <- get_metric_cell_types(metric)
  
  # Define plotting function
  formatCustomSci <- function(x) {
    x_sci <- str_split_fixed(formatC(x, format = "e"), "e", 2)
    alpha <- as.numeric(x_sci[ , 1])
    power <- as.integer(x_sci[ , 2])
    paste(alpha, power, sep = "e")
  }
  
  create_plot <- function(data, x_aes, y_aes, color_aes, title = "") {
    
    size <- 0.5
    
    data <- data[data$variable_parameter == color_aes, ]
    
    plot <- ggplot(data, aes_string(x = x_aes, y = y_aes, color = color_aes)) +
      labs(title = title, x = x_aes, y = y_aes) +
      theme_bw() +
      xlim(min(c(data[[x_aes]], data[[y_aes]])), max(c(data[[x_aes]], data[[y_aes]]))) +
      ylim(min(c(data[[x_aes]], data[[y_aes]])), max(c(data[[x_aes]], data[[y_aes]]))) +
      geom_abline(intercept = 0, slope = 1, color = "black", linetype = "longdash")
    
    
    # Use scientific notation for ellipsoid volume
    if (color_aes == "E_volume") {
      plot <- plot + 
        geom_point(size = size) +
        scale_color_continuous(labels = formatCustomSci)
    }
    else if (typeof(data[[color_aes]]) == "double") {
      breaks <- pretty(c(min(data[[color_aes]]), max(data[[color_aes]])), n = 3)
      plot <- plot + 
        geom_point(size = size) + 
        scale_color_continuous(breaks = breaks)
    }
    return(plot)
  }
  
  # Put plots into an organised list
  plots_list <- list()
  
  for (i in seq(nrow(metric_cell_types))) {
    
    # Subset metric_df for chosen pair/cell types
    metric_df3D_subset <- subset_metric_df(metric, metric_df3D, metric_cell_types, i)
    metric_df2D_subset <- subset_metric_df(metric, metric_df2D, metric_cell_types, i)
    
    # Combine with spes table
    plot_df <- spes_table
    plot_df[[paste(metric, "3D", sep = "_")]] <- metric_df3D_subset[[metric]]
    plot_df[[paste(metric, "2D", sep = "_")]] <- metric_df2D_subset[[metric]]
    
    # Factor
    if (!is.null(plot_df$shape)) plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    if (!is.null(plot_df$slice)) plot_df$slice <- as.character(plot_df$slice)
    
    # Generate plots based on plots_metadata, use final column of metric_cell_types
    plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]] <- lapply(plots_metadata, function(plot_def) {
      x_aes <- paste(metric, "3D", sep = "_")
      y_aes <- paste(metric, "2D", sep = "_")
      color_aes <- plot_def$color_aes
      title <- plot_def$title
      plot <- create_plot(data = plot_df, x_aes = x_aes, y_aes = y_aes, color_aes = color_aes, title = title)
      return(plot)
    })
  }
  
  # Extract legends from first set of plots
  legends_list <- lapply(plots_list[[1]], function(plot) {
    plot_legend <- get_legend(plot + theme(legend.direction = "horizontal"))
    return(plot_legend)
  })
  legends <- plot_grid(plotlist = legends_list, nrow = 1)
  
  
  # Combine the plots together using metric_cell_types
  combined_plots_list <- list()
  for (i in seq(nrow(metric_cell_types))) {
    
    # Remove legend from base plots
    for (j in seq(length(plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]]))) {
      plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]][[j]] <- 
        plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]][[j]] + theme(legend.position = "none")
    }
    
    # Getting current set of cell types from metric_cell_types
    cells <- metric_cell_types[i, ncol(metric_cell_types)]
    
    plots <- plot_grid(plotlist = plots_list[[cells]], nrow = 1, ncol = length(plots_list[[cells]]))
    
    title <- get_metric_cell_types_title(metric, metric_cell_types, i)
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    combined_plots_list[[cells]] <- fig
  }
  
  # Combine the combined plots into one big plot
  combined_plots <- plot_grid(plotlist = combined_plots_list,
                              nrow = length(combined_plots_list), 
                              ncol = 1)
  
  combined_plots_with_legends <- plot_grid(combined_plots, 
                                           legends,
                                           nrow = 2, ncol = 1,
                                           rel_heights = c(1, 0.1))
  return(combined_plots_with_legends)
}


### Function to compare 3D vs 2D all slices ------------------------------------------
plot_3D_vs_2D_metric_all_slices <- function(spes_table, 
                                            metric, 
                                            metric_df3D,
                                            metric_df2D,
                                            arrangement, 
                                            plots_metadata) {
  
  ### Modify plots_metadata
  # Change plots_metadata arrangement to inputted arrangement
  plots_metadata$arrangement$label <- arrangement  
  
  # Get metric_cell_types
  metric_cell_types <- get_metric_cell_types(metric)
  
  # Get number of slices
  n_slices <- length(unique(metric_df2D[["slice"]]))
  
  # Define plotting function
  formatCustomSci <- function(x) {
    x_sci <- str_split_fixed(formatC(x, format = "e"), "e", 2)
    alpha <- as.numeric(x_sci[ , 1])
    power <- as.integer(x_sci[ , 2])
    paste(alpha, power, sep = "e")
  }
  
  create_plot <- function(data, x_aes, y_aes, color_aes, label, title = "") {
    
    size <- 0.5
    
    data <- data[data$variable_parameter == label, ]
    
    plot <- ggplot(data, aes_string(x = x_aes, y = y_aes, color = color_aes)) +
      labs(title = title, x = x_aes, y = y_aes) +
      theme_bw() +
      xlim(min(c(data[[x_aes]], data[[y_aes]]), na.rm = T), max(c(data[[x_aes]], data[[y_aes]]), na.rm = T)) +
      ylim(min(c(data[[x_aes]], data[[y_aes]]), na.rm = T), max(c(data[[x_aes]], data[[y_aes]]), na.rm = T)) +
      geom_point(size = 0.5) +
      geom_abline(intercept = 0, slope = 1, color = "black", linetype = "longdash")
    
    return(plot)
  }
  
  # Put plots into an organised list
  plots_list <- list()
  
  for (i in seq(nrow(metric_cell_types))) {
    
    # Subset metric_df for chosen pair/cell types
    metric_df3D_subset <- subset_metric_df(metric, metric_df3D, metric_cell_types, i)
    metric_df2D_subset <- subset_metric_df(metric, metric_df2D, metric_cell_types, i)
    
    # Combine with spes table
    plot_df <- spes_table
    plot_df[[paste(metric, "3D", sep = "_")]] <- metric_df3D_subset[[metric]]
    
    plot_df <- duplicate_df(plot_df, n_slices)
    
    plot_df[[paste(metric, "2D", sep = "_")]] <- metric_df2D_subset[[metric]]
    plot_df[["slice"]] <- metric_df2D_subset[["slice"]]
    
    # Factor
    if (!is.null(plot_df$shape)) plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    if (!is.null(plot_df$slice)) plot_df$slice <- as.character(plot_df$slice)
    
    # Generate plots based on plots_metadata, use final column of metric_cell_types
    plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]] <- lapply(plots_metadata, function(plot_def) {
      x_aes <- paste(metric, "3D", sep = "_")
      y_aes <- paste(metric, "2D", sep = "_")
      color_aes <- plot_def$color_aes
      label <- plot_def$label
      title <- plot_def$title
      plot <- create_plot(data = plot_df, x_aes = x_aes, y_aes = y_aes, color_aes = color_aes, label = label, title = title)
      return(plot)
    })
  }
  
  # Extract legends from first set of plots
  legends_list <- lapply(plots_list[[1]], function(plot) {
    plot_legend <- get_legend(plot + theme(legend.direction = "horizontal"))
    return(plot_legend)
  })
  legends <- plot_grid(plotlist = legends_list, nrow = 1)
  
  
  # Combine the plots together using metric_cell_types
  combined_plots_list <- list()
  for (i in seq(nrow(metric_cell_types))) {
    
    # Remove legend from base plots
    for (j in seq(length(plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]]))) {
      plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]][[j]] <- 
        plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]][[j]] + theme(legend.position = "none")
    }
    
    # Getting current set of cell types from metric_cell_types
    cells <- metric_cell_types[i, ncol(metric_cell_types)]
    
    plots <- plot_grid(plotlist = plots_list[[cells]], nrow = 1, ncol = length(plots_list[[cells]]))
    
    title <- get_metric_cell_types_title(metric, metric_cell_types, i)
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    combined_plots_list[[cells]] <- fig
  }
  
  # Combine the combined plots into one big plot
  combined_plots <- plot_grid(plotlist = combined_plots_list,
                              nrow = length(combined_plots_list), 
                              ncol = 1)
  
  # Get labels, if specified
  labels_vec <- unlist(lapply(plots_metadata, function(x) {
    return(x$label)
  }))
  labels <- list()
  for (label in labels_vec) {
    label <- ggdraw() +
      draw_label(label)
    labels[[length(labels) + 1]] <- label
  }
  labels <- plot_grid(plotlist = labels, nrow = 1)
  combined_plots_with_legends_and_labels <- plot_grid(combined_plots, 
                                                      legends,
                                                      labels,
                                                      nrow = 3, ncol = 1,
                                                      rel_heights = c(1, 0.1, 0.1))
  return(combined_plots_with_legends_and_labels)
}

### Function for non_gradient 3D vs 2D ERROR output ----
plot_error_non_gradient_metric <- function(spes_table, 
                                           metric, 
                                           metric_df3D,
                                           metric_df2D,
                                           arrangement, 
                                           plots_metadata) {
  
  ### Modify plots_metadata
  # Change plots_metadata arrangement to inputted arrangement
  plots_metadata$arrangement$x_aes <- arrangement
  
  # Get metric_cell_types
  metric_cell_types <- get_metric_cell_types(metric)
  
  # Get number of slices
  n_slices <- length(unique(metric_df2D[["slice"]]))
  
  # Define plotting function
  formatCustomSci <- function(x) {
    x_sci <- str_split_fixed(formatC(x, format = "e"), "e", 2)
    alpha <- as.numeric(x_sci[ , 1])
    power <- as.integer(x_sci[ , 2])
    paste(alpha, power, sep = "e")
  }
  
  create_plot <- function(data, x_aes, y_aes, color_aes, title = "") {
    
    size <- 0.5
    
    data <- data[data$variable_parameter == x_aes, ]
    
    plot <- ggplot(data, aes_string(x = x_aes, y = y_aes, color = color_aes)) +
      labs(title = title, x = x_aes, y = y_aes) +
      theme_bw() +
      ylim(min(c(0, data[[y_aes]]), na.rm = T), max(c(0, data[[y_aes]]), na.rm = T))
    
    # Use scientific notation for ellipsoid volume
    if (x_aes == "E_volume") {
      plot <- plot + 
        geom_point(size = size) +
        geom_abline(intercept = 0, slope = 0, color = "black", linetype = "longdash") +
        scale_x_continuous(labels = formatCustomSci)
    }
    else if (typeof(data[[x_aes]]) == "double") {
      breaks <- pretty(c(min(data[[x_aes]]), max(data[[x_aes]])), n = 3)
      plot <- plot + 
        geom_point(size = size) + 
        geom_abline(intercept = 0, slope = 0, color = "black", linetype = "longdash") +
        scale_x_continuous(breaks = breaks)
    }
    return(plot)
  }
  
  # Put plots into an organised list
  plots_list <- list()
  
  for (i in seq(nrow(metric_cell_types))) {
    
    # Subset metric_df for chosen pair/cell types
    metric_df3D_subset <- subset_metric_df(metric, metric_df3D, metric_cell_types, i)
    metric_df2D_subset <- subset_metric_df(metric, metric_df2D, metric_cell_types, i)
    
    # Combine with spes table
    plot_df <- spes_table
    plot_df[[paste(metric, "3D", sep = "_")]] <- metric_df3D_subset[[metric]]
    
    plot_df <- duplicate_df(plot_df, n_slices)
    
    plot_df[[paste(metric, "2D", sep = "_")]] <- metric_df2D_subset[[metric]]
    plot_df[["slice"]] <- metric_df2D_subset[["slice"]]
    
    # Error: ((2D - 3D) / 3D) * 100%
    plot_df[[paste(metric, "error", sep = "_")]] <- 100 * (plot_df[[paste(metric, "2D", sep = "_")]] - plot_df[[paste(metric, "3D", sep = "_")]]) / (plot_df[[paste(metric, "3D", sep = "_")]])
    
    # Factor
    if (!is.null(plot_df$shape)) plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    if (!is.null(plot_df$slice)) plot_df$slice <- as.character(plot_df$slice)
    
    # Generate plots based on plots_metadata, use final column of metric_cell_types
    plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]] <- lapply(plots_metadata, function(plot_def) {
      x_aes <- plot_def$x_aes
      y_aes <- paste(metric, "error", sep = "_")
      color_aes <- plot_def$color_aes
      title <- plot_def$title
      plot <- create_plot(data = plot_df, x_aes = x_aes, y_aes = y_aes, color_aes = color_aes, title = title)
      return(plot)
    })
  }
  
  # Extract legends from first set of plots
  legends_list <- lapply(plots_list[[1]], function(plot) {
    plot_legend <- get_legend(plot + theme(legend.direction = "horizontal"))
    return(plot_legend)
  })
  legends <- plot_grid(plotlist = legends_list, nrow = 1)
  
  
  # Combine the plots together using metric_cell_types
  combined_plots_list <- list()
  for (i in seq(nrow(metric_cell_types))) {
    
    # Remove legend from base plots
    for (j in seq(length(plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]]))) {
      plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]][[j]] <- 
        plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]][[j]] + theme(legend.position = "none")
    }
    
    # Getting current set of cell types from metric_cell_types
    cells <- metric_cell_types[i, ncol(metric_cell_types)]
    
    plots <- plot_grid(plotlist = plots_list[[cells]], nrow = 1, ncol = length(plots_list[[cells]]))
    
    title <- get_metric_cell_types_title(metric, metric_cell_types, i)
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    combined_plots_list[[cells]] <- fig
  }
  
  # Combine the combined plots into one big plot
  combined_plots <- plot_grid(plotlist = combined_plots_list,
                              nrow = length(combined_plots_list), 
                              ncol = 1)
  
  combined_plots_with_legends <- plot_grid(combined_plots, 
                                           legends,
                                           nrow = 2, ncol = 1,
                                           rel_heights = c(1, 0.1))
  return(combined_plots_with_legends)
}

### Function for gradient 3D vs 2D ERROR output one slice ----
plot_error_gradient_metric_one_slice <- function(spes_table, 
                                                 metric, 
                                                 metric_df3D,
                                                 metric_df2D, 
                                                 arrangement, 
                                                 gradient_type,
                                                 plots_metadata) {
  ### Modify plots_metadata
  # Change plots_metadata arrangement to inputted arrangement
  plots_metadata$arrangement$color_aes <- arrangement
  
  # Change plots_metadata x_aes to gradient_type and y_aes to inputted metric
  for (i in seq_along(plots_metadata)) {
    plots_metadata[[i]]$x_aes <- gradient_type
  }
  
  # Get gradient/threshold values
  if (gradient_type == "radius") {
    gradient <- seq(20, 100, 10)
    gradient_colnames <- paste("r", gradient, sep = "")    
  }
  else if (gradient_type == "threshold") {
    gradient <- seq(0.01, 1, 0.01)
    gradient_colnames <- paste("t", gradient, sep = "")
  }
  
  # Get metric_cell_types
  metric_cell_types <- get_metric_cell_types(metric)
  
  
  # Define plotting function
  formatCustomSci <- function(x) {
    x_sci <- str_split_fixed(formatC(x, format = "e"), "e", 2)
    alpha <- as.numeric(x_sci[ , 1])
    power <- as.integer(x_sci[ , 2])
    paste(alpha, power, sep = "e")
  }
  create_plot <- function(data, x_aes, y_aes, color_aes, group_aes, title = "") {
    
    data <- data[data$variable_parameter == color_aes, ]
    
    breaks <- pretty(c(min(data[[color_aes]]), max(data[[color_aes]])), n = 3)
    
    plot <- ggplot(data, aes_string(x = x_aes, y = y_aes, group = group_aes, color = color_aes)) +
      labs(title = title, x = x_aes, y = y_aes) +
      theme_bw() +
      geom_line() +
      geom_abline(intercept = 0, slope = 0, color = "black", linetype = "longdash")
    
    # Use scientific notation for ellipsoid volume
    if (color_aes == "E_volume") {
      plot <- plot + scale_color_continuous(labels = formatCustomSci)
    }
    else {
      plot <- plot + scale_color_continuous(breaks = breaks)
    }
    
    return(plot)
  }
  
  # Put plots into an organised list
  plots_list <- list()
  
  for (i in seq(nrow(metric_cell_types))) {
    
    # Subset metric_df for chosen pair/cell types
    metric_df3D_subset <- subset_metric_df(metric, metric_df3D, metric_cell_types, i)
    metric_df2D_subset <- subset_metric_df(metric, metric_df2D, metric_cell_types, i)
    
    # Get error df
    error_df <- ((metric_df2D_subset[, gradient_colnames] - metric_df3D_subset[, gradient_colnames]) / metric_df3D_subset[, gradient_colnames]) * 100
    error_df[["spe"]] <- metric_df3D_subset[["spe"]]
    
    # Combine spes_table and error_df
    plot_df <- cbind(spes_table, error_df)
    
    # Melt
    plot_df <- reshape2::melt(plot_df, , gradient_colnames)
    
    # Change last 2 column names
    colnames(plot_df)[c(ncol(plot_df) - 1, ncol(plot_df))] <- c(gradient_type, paste(metric, "error", sep = "_"))
    
    # Extract radius value from radius strings (r1 -> 1, r2 -> 2...)
    plot_df[[gradient_type]] <- unfactor(plot_df[[gradient_type]])
    plot_df[[gradient_type]] <- as.numeric(substr(plot_df[[gradient_type]], 2, nchar(plot_df[[gradient_type]])))
    
    # Factor
    if (!is.null(plot_df$shape)) plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    
    if (!is.null(plot_df$slice)) {
      plot_df$slice <- as.character(plot_df$slice)
      plot_df$key <- paste(plot_df$spe, plot_df$slice, sep = "_")
      group_aes = "key"
    }
    else {
      group_aes = "spe"
    }
    
    # Generate plots based on plots_metadata, use final column of metric_cell_types
    plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]] <- lapply(plots_metadata, function(plot_def) {
      x_aes <- plot_def$x_aes
      y_aes <- paste(metric, "error", sep = "_")
      color_aes <- plot_def$color_aes
      title <- plot_def$title
      plot <- create_plot(data = plot_df, x_aes = x_aes, y_aes = y_aes, group_aes = group_aes, color_aes = color_aes, title = title)
      return(plot)
    })
  }
  
  # Extract legends from first set of plots
  legends_list <- lapply(plots_list[[1]], function(plot) {
    plot_legend <- get_legend(plot + theme(legend.direction = "horizontal"))
    return(plot_legend)
  })
  legends <- plot_grid(plotlist = legends_list, nrow = 1)
  
  # Combine the plots together using metric_cell_types
  combined_plots_list <- list()
  for (i in seq(nrow(metric_cell_types))) {
    
    # Remove legend from base plots
    for (j in seq(length(plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]]))) {
      plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]][[j]] <- 
        plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]][[j]] + theme(legend.position = "none")
    }
    
    # Getting current set of cell types from metric_cell_types
    cells <- metric_cell_types[i, ncol(metric_cell_types)]
    
    plots <- plot_grid(plotlist = plots_list[[cells]], nrow = 1, ncol = length(plots_list[[cells]]))
    
    title <- get_metric_cell_types_title(metric, metric_cell_types, i)
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    combined_plots_list[[cells]] <- fig
  }
  
  # Combine the combined plots into one big plot
  gradient_metric_plot <- plot_grid(plotlist = combined_plots_list,
                                    nrow = length(combined_plots_list), ncol = 1)
  
  # Add legends
  gradient_metric_with_legends_plot <- plot_grid(gradient_metric_plot, legends,
                                                 nrow = 2, ncol = 1,
                                                 rel_heights = c(1, 0.1))
  
  return(gradient_metric_with_legends_plot)
}

### Function for non-gradient output all slices with ground truth ----
plot_non_gradient_metric_all_slices_ground_truth <- function(spes_table, 
                                                             metric, 
                                                             metric_df3D, 
                                                             metric_df2D, 
                                                             arrangement, 
                                                             plots_metadata) {
  
  ### Modify plots_metadata
  # Change plots_metadata arrangement to inputted arrangement
  plots_metadata$arrangement$x_aes <- arrangement
  
  # Change plots_metadata y_aes to inputted metric
  for (i in seq_along(plots_metadata)) {
    # Modify the y_aes element
    plots_metadata[[i]]$y_aes <- metric
  }
  
  # Get metric_cell_types
  metric_cell_types <- get_metric_cell_types(metric)
  
  # Define plotting function
  formatCustomSci <- function(x) {
    x_sci <- str_split_fixed(formatC(x, format = "e"), "e", 2)
    alpha <- as.numeric(x_sci[ , 1])
    power <- as.integer(x_sci[ , 2])
    paste(alpha, power, sep = "e")
  }
  
  create_plot <- function(data2D, data3D, x_aes, y_aes, color_aes, title = "") {
    
    size <- 0.5
    
    data3D <- data3D[data3D$variable_parameter == x_aes, ]
    data2D <- data2D[data2D$variable_parameter == x_aes, ]
    
    plot <- ggplot() +
      labs(title = title, x = x_aes, y = y_aes) +
      geom_point(data = data3D, aes_string(x = x_aes, y = y_aes), size = size) +
      theme_bw()
    
    # Use scientific notation for ellipsoid volume
    if (x_aes == "E_volume") {
      plot <- plot + 
        geom_point(data = data2D, aes_string(x = x_aes, y = y_aes, color = color_aes), size = size) +
        scale_x_continuous(labels = formatCustomSci)
    }
    else if (typeof(data2D[[x_aes]]) == "double") {
      breaks <- pretty(c(min(data2D[[x_aes]]), max(data2D[[x_aes]])), n = 2)
      plot <- plot + 
        geom_point(data = data2D, aes_string(x = x_aes, y = y_aes, color = color_aes), size = size) +
        scale_x_continuous(breaks = breaks)
    }
    
    return(plot)
  }
  
  # Put plots into an organised list
  plots_list <- list()
  
  # Get number of slices
  n_slices <- length(unique(metric_df2D[["slice"]]))
  
  for (i in seq(nrow(metric_cell_types))) {
    
    # Subset metric_df for chosen pair/cell types
    metric_df3D_subset <- subset_metric_df(metric, metric_df3D, metric_cell_types, i)
    metric_df2D_subset <- subset_metric_df(metric, metric_df2D, metric_cell_types, i)
    
    # Combine with spes table
    plot_df3D <- spes_table
    plot_df3D[[metric]] <- metric_df3D_subset[[metric]]
    
    plot_df2D <- duplicate_df(spes_table, n_slices)
    
    plot_df2D[[metric]] <- metric_df2D_subset[[metric]]
    plot_df2D[["slice"]] <- metric_df2D_subset[["slice"]]
    
    # Factor
    if (!is.null(plot_df2D$shape)) plot_df2D$shape <- factor(plot_df2D$shape, c("Ellipsoid", "Network"))
    if (!is.null(plot_df2D$slice)) plot_df2D$slice <- as.character(plot_df2D$slice)
    
    if (!is.null(plot_df3D$shape)) plot_df3D$shape <- factor(plot_df3D$shape, c("Ellipsoid", "Network"))
    if (!is.null(plot_df3D$slice)) plot_df3D$slice <- as.character(plot_df3D$slice)
    
    # Generate plots based on plots_metadata, use final column of metric_cell_types
    plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]] <- lapply(plots_metadata, function(plot_def) {
      x_aes <- plot_def$x_aes
      y_aes <- plot_def$y_aes
      color_aes <- plot_def$color_aes
      title <- plot_def$title
      plot <- create_plot(data2D = plot_df2D, data3D = plot_df3D, x_aes = x_aes, y_aes = y_aes, color_aes = color_aes, title = title)
      return(plot)
    })
  }
  
  # Extract legends from first set of plots
  legends_list <- lapply(plots_list[[1]], function(plot) {
    plot_legend <- get_legend(plot + theme(legend.direction = "horizontal"))
    return(plot_legend)
  })
  legends <- plot_grid(plotlist = legends_list, nrow = 1)
  
  # Combine the plots together using metric_cell_types
  combined_plots_list <- list()
  for (i in seq(nrow(metric_cell_types))) {
    
    # Remove legend from base plots
    for (j in seq(length(plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]]))) {
      
      plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]][[j]] <- 
        plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]][[j]] + theme(legend.position = "none")
    }
    
    # Get final column
    cells <- metric_cell_types[i, ncol(metric_cell_types)]
    
    plots <- plot_grid(plotlist = plots_list[[cells]], nrow = 1, ncol = length(plots_list[[cells]]))
    
    title <- get_metric_cell_types_title(metric, metric_cell_types, i)
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    combined_plots_list[[cells]] <- fig
  }
  
  # Combine the combined plots into one big plot
  non_gradient_metric_plot <- plot_grid(plotlist = combined_plots_list,
                                        nrow = length(combined_plots_list), 
                                        ncol = 1)
  
  non_gradient_metric_with_legends_plot <- plot_grid(non_gradient_metric_plot,
                                                     legends,
                                                     nrow = 2,
                                                     ncol = 1,
                                                     rel_heights = c(1, 0.1))
  return(non_gradient_metric_with_legends_plot)
}

### Function for slice violin plots -------
plot_violin_all_slices <- function(spes_table,
                                   metric,
                                   metric_df2D,
                                   arrangement,
                                   plots_metadata) {
  
  ### Modify plots_metadata
  # Change plots_metadata arrangement to inputted arrangement
  plots_metadata$arrangement$label <- arrangement
  
  # Change plots_metadata y_aes to inputted metric
  for (i in seq_along(plots_metadata)) {
    # Modify the y_aes element
    plots_metadata[[i]]$y_aes <- metric
  }
  
  # Get metric_cell_types
  metric_cell_types <- get_metric_cell_types(metric)
  
  # Define plotting function
  formatCustomSci <- function(x) {
    x_sci <- str_split_fixed(formatC(x, format = "e"), "e", 2)
    alpha <- round(as.numeric(x_sci[ , 1]), 1)
    power <- as.integer(x_sci[ , 2])
    paste(alpha, power, sep = "e")
  }
  
  create_plot <- function(data, x_aes, y_aes, label, title = "") {
    
    data <- data[data$variable_parameter == label, ]
    
    do_JT_test <- FALSE
    
    count <- 0
    for (i in unique(data[[x_aes]])) {
      if (sum(!is.na(data[[y_aes]][data[[x_aes]] == i])) >= 2) count <- count + 1
    }
    if (count >= 2) do_JT_test <- TRUE
    if (all(data[[y_aes]] == 0, na.rm = T)) do_JT_test <- FALSE
    
    # JT test
    if (do_JT_test) {
      data_means <- aggregate(data[[y_aes]], by = list(data[[x_aes]]), FUN = mean, na.rm = T)
      
      colnames(data_means) <- c("slice", "mean_value")
      data_means <- data_means %>%
        mutate(diff = mean_value - lag(mean_value))
      trend_direction <- ifelse(mean(data_means$diff, na.rm = T) > 0, "increasing", ifelse(mean(data_means$diff, na.rm = T) < 0, "decreasing", "two.sided"))
      jt_results <- JonckheereTerpstraTest(data[[x_aes]], data[[y_aes]], alternative = trend_direction)
      p_value <- jt_results$p.value
      if (p_value == 0) p_value <- 2.2e-308
      if (0 < p_value && p_value < 1e-3)  {
        p_value <- formatCustomSci(p_value)
      }
      else {
        p_value  <- round(p_value, 3)
      }
      title <- paste("p:", p_value)
    }
    
    plot <- ggplot() +
      labs(x = x_aes, y = y_aes) +
      theme_bw() +
      ggtitle(title) +
      theme(plot.title = element_text(size = 9)) +
      geom_violin(data = data, aes_string(x = x_aes, y = y_aes))
  }
  
  # Put plots into an organised list
  plots_list <- list()
  
  # Get number of slices
  n_slices <- length(unique(metric_df2D[["slice"]]))
  
  for (i in seq(nrow(metric_cell_types))) {
    
    # Subset metric_df for chosen pair/cell types
    metric_df2D_subset <- subset_metric_df(metric, metric_df2D, metric_cell_types, i)
    
    # Combine with spes table
    plot_df <- duplicate_df(spes_table, n_slices)
    
    plot_df[[metric]] <- metric_df2D_subset[[metric]]
    plot_df[["slice"]] <- metric_df2D_subset[["slice"]]
    
    
    # Factor
    if (!is.null(plot_df$shape)) plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    if (!is.null(plot_df$slice)) plot_df$slice <- as.character(plot_df$slice)
    
    # Generate plots based on plots_metadata, use final column of metric_cell_types
    plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]] <- lapply(plots_metadata, function(plot_def) {
      x_aes <- plot_def$x_aes
      y_aes <- plot_def$y_aes
      label <- plot_def$label
      title <- " "
      plot <- create_plot(data = plot_df, x_aes = x_aes, y_aes = y_aes, label = label, title = title)
      return(plot)
    })
  }
  
  # Combine the plots together using metric_cell_types
  combined_plots_list <- list()
  for (i in seq(nrow(metric_cell_types))) {
    
    # Get final column
    cells <- metric_cell_types[i, ncol(metric_cell_types)]
    
    plots <- plot_grid(plotlist = plots_list[[cells]], nrow = 1, ncol = length(plots_list[[cells]]))
    
    title <- get_metric_cell_types_title(metric, metric_cell_types, i)
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    combined_plots_list[[cells]] <- fig
  }
  
  # Combine the combined plots into one big plot
  combined_plots <- plot_grid(plotlist = combined_plots_list,
                              nrow = length(combined_plots_list), 
                              ncol = 1)
  
  # Get labels, if specified
  labels_vec <- unlist(lapply(plots_metadata, function(x) {
    return(x$label)
  }))
  labels <- list()
  for (label in labels_vec) {
    label <- ggdraw() +
      draw_label(label)
    labels[[length(labels) + 1]] <- label
  }
  labels <- plot_grid(plotlist = labels, nrow = 1)
  combined_plots_with_labels <- plot_grid(combined_plots, 
                                          labels,
                                          nrow = 3, ncol = 1,
                                          rel_heights = c(1, 0.1, 0.1))
  return(combined_plots_with_labels)
}

### Function for slice violin plots with ground truth -------
plot_violin_all_slices_ground_truth <- function(spes_table,
                                                metric,
                                                metric_df3D,
                                                metric_df2D,
                                                arrangement,
                                                plots_metadata) {
  
  ### Modify plots_metadata
  # Change plots_metadata arrangement to inputted arrangement
  plots_metadata$arrangement$label <- arrangement
  
  # Change plots_metadata y_aes to inputted metric
  for (i in seq_along(plots_metadata)) {
    # Modify the y_aes element
    plots_metadata[[i]]$y_aes <- metric
  }
  
  # Get metric_cell_types
  metric_cell_types <- get_metric_cell_types(metric)
  
  # Define plotting function
  formatCustomSci <- function(x) {
    x_sci <- str_split_fixed(formatC(x, format = "e"), "e", 2)
    alpha <- as.numeric(x_sci[ , 1])
    power <- as.integer(x_sci[ , 2])
    paste(alpha, power, sep = "e")
  }
  
  create_plot <- function(data, x_aes, y_aes, label, title = "") {
    
    data <- data[data$variable_parameter == label, ]
    
    plot <- ggplot() +
      labs(title = title, x = x_aes, y = y_aes) +
      geom_violin(data = data, aes_string(x = x_aes, y = y_aes)) +
      theme_bw()
    
    return(plot)
  }
  
  # Put plots into an organised list
  plots_list <- list()
  
  # Get number of slices
  n_slices <- length(unique(metric_df2D[["slice"]]))
  
  for (i in seq(nrow(metric_cell_types))) {
    
    # Subset metric_df for chosen pair/cell types
    metric_df3D_subset <- subset_metric_df(metric, metric_df3D, metric_cell_types, i)
    metric_df2D_subset <- subset_metric_df(metric, metric_df2D, metric_cell_types, i)
    
    # Combine with spes table
    plot_df3D <- spes_table
    plot_df3D[[metric]] <- metric_df3D_subset[[metric]]
    plot_df3D[["slice"]] <- "GT"
    
    plot_df2D <- duplicate_df(spes_table, n_slices)
    
    plot_df2D[[metric]] <- metric_df2D_subset[[metric]]
    plot_df2D[["slice"]] <- metric_df2D_subset[["slice"]]
    
    plot_df <- rbind(plot_df3D, plot_df2D)
    
    # Factor
    if (!is.null(plot_df$shape)) plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    if (!is.null(plot_df$slice)) plot_df$slice <- as.character(plot_df$slice)
    
    # Generate plots based on plots_metadata, use final column of metric_cell_types
    plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]] <- lapply(plots_metadata, function(plot_def) {
      x_aes <- plot_def$x_aes
      y_aes <- plot_def$y_aes
      label <- plot_def$label
      title <- plot_def$title
      plot <- create_plot(data = plot_df, x_aes = x_aes, y_aes = y_aes, label = label, title = title)
      return(plot)
    })
  }
  
  # Combine the plots together using metric_cell_types
  combined_plots_list <- list()
  for (i in seq(nrow(metric_cell_types))) {
    
    # Get final column
    cells <- metric_cell_types[i, ncol(metric_cell_types)]
    
    plots <- plot_grid(plotlist = plots_list[[cells]], nrow = 1, ncol = length(plots_list[[cells]]))
    
    title <- get_metric_cell_types_title(metric, metric_cell_types, i)
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    combined_plots_list[[cells]] <- fig
  }
  
  # Combine the combined plots into one big plot
  combined_plots <- plot_grid(plotlist = combined_plots_list,
                              nrow = length(combined_plots_list), 
                              ncol = 1)
  
  # Get labels, if specified
  labels_vec <- unlist(lapply(plots_metadata, function(x) {
    return(x$label)
  }))
  labels <- list()
  for (label in labels_vec) {
    label <- ggdraw() +
      draw_label(label)
    labels[[length(labels) + 1]] <- label
  }
  labels <- plot_grid(plotlist = labels, nrow = 1)
  combined_plots_with_labels <- plot_grid(combined_plots, 
                                          labels,
                                          nrow = 3, ncol = 1,
                                          rel_heights = c(1, 0.1, 0.1))
  return(combined_plots_with_labels)
}

library(cowplot)
library(ggplot2)
library(S4Vectors)
library(stringr)
library(dplyr)

### Utility functions -------------------------
get_gradient <- function(metric) {
  if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "CLR", "COO", "CGR", "CGR")) return("radius")
  return("threshold")
}

### Read metric_df_lists --------------------------------------------------------------
setwd("~/R/SPIAT-3D_benchmarking/simulations_and_analysis_S1/S1_data")
metric_df_lists3D <- readRDS("metric_df_lists3D.RDS")
metric_df_lists2D <- readRDS("metric_df_lists2D.RDS")

### Turn gradient radii metrics into AUC and add to metric_df list ------------------
get_AUC_for_radii_gradient_metrics <- function(y) {
  x <- radii
  h <- diff(x)[1]
  n <- length(x)
  
  AUC <- (h / 2) * (y[1] + 2 * sum(y[2:(n - 1)]) + y[n])
  
  return(AUC)
}

arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")

radii <- seq(20, 100, 10)
radii_colnames <- paste("r", radii, sep = "")

gradient_radii_metrics <- c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "CLR", "COO", "CGR", "CGR")


for (arrangement in arrangements) {
  
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    for (metric in gradient_radii_metrics) {
      metric_AUC_name <- paste(metric, "AUC", sep = "_")
      
      if (metric %in% c("MS", "NMS", "ACIN", "CKR", "CLR", "COO", "CGR", "CGR")) {
        subset_colnames <- c("spe", "reference", "target", metric_AUC_name)
      }
      else {
        subset_colnames <- c("spe", "reference", metric_AUC_name)
      }
      
      # 3D
      df <- metric_df_lists3D[[spes_metadata_index]][[metric]]
      df[[metric_AUC_name]] <- apply(df[ , radii_colnames], 1, get_AUC_for_radii_gradient_metrics)
      
      df <- df[ , subset_colnames]
      metric_df_lists3D[[spes_metadata_index]][[metric_AUC_name]] <- df
      
      # 2D
      subset_colnames <- c(subset_colnames, "slice")
      df <- metric_df_lists2D[[spes_metadata_index]][[metric]]
      df[[metric_AUC_name]] <- apply(df[ , radii_colnames], 1, get_AUC_for_radii_gradient_metrics)
      
      df <- df[ , subset_colnames]
      metric_df_lists2D[[spes_metadata_index]][[metric_AUC_name]] <- df
      
    }
  }
}


### Turn threshold radii metrics into AUC and add to metric_df list--------------
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")

thresholds <- seq(0.01, 1, 0.01)
threshold_colnames <- paste("t", thresholds, sep = "")

for (arrangement in arrangements) {
  
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    # PBP_AUC 3D
    PBP_df <- metric_df_lists3D[[spes_metadata_index]][["PBP"]]
    PBP_df$PBP_AUC <- apply(PBP_df[ , threshold_colnames], 1, sum) * 0.01
    PBP_AUC_df <- PBP_df[ , c("spe", "reference", "target", "PBP_AUC")]
    metric_df_lists3D[[spes_metadata_index]][["PBP_AUC"]] <- PBP_AUC_df
    
    # EBP_AUC 3D
    EBP_df <- metric_df_lists3D[[spes_metadata_index]][["EBP"]]
    EBP_df$EBP_AUC <- apply(EBP_df[ , threshold_colnames], 1, sum) * 0.01
    EBP_AUC_df <- EBP_df[ , c("spe", "cell_types", "EBP_AUC")]
    metric_df_lists3D[[spes_metadata_index]][["EBP_AUC"]] <- EBP_AUC_df
    
    # PBP_AUC 2D
    PBP_df <- metric_df_lists2D[[spes_metadata_index]][["PBP"]]
    PBP_df$PBP_AUC <- apply(PBP_df[ , threshold_colnames], 1, sum) * 0.01
    PBP_AUC_df <- PBP_df[ , c("spe", "slice", "reference", "target", "PBP_AUC")]
    metric_df_lists2D[[spes_metadata_index]][["PBP_AUC"]] <- PBP_AUC_df
    
    # EBP_AUC 2D
    EBP_df <- metric_df_lists2D[[spes_metadata_index]][["EBP"]]
    EBP_df$EBP_AUC <- apply(EBP_df[ , threshold_colnames], 1, sum) * 0.01
    EBP_AUC_df <- EBP_df[ , c("spe", "slice", "cell_types", "EBP_AUC")]
    metric_df_lists2D[[spes_metadata_index]][["EBP_AUC"]] <- EBP_AUC_df
  }
}


### Get plots for 3D metric analysis (non-gradient) ------------------------------------------
# Read spes_table
setwd("~/R/SPIAT-3D_benchmarking/simulations_and_analysis_S1/S1_data")
spes_table <- read.table("spes_table.csv")
spes_table$cluster_prop_B <- 1 - spes_table$cluster_prop_A
spes_table[spes_table$variable_parameter == "cluster_prop_A", "variable_parameter"] <- "cluster_prop_B"
spes_table$distance <- 450 - spes_table$cluster_x_coord
spes_table[spes_table$variable_parameter == "cluster_x_coord", "variable_parameter"] <- "distance" 
spes_table$E_volume <- (4/3) * pi * spes_table$E_radius_x * spes_table$E_radius_y * spes_table$E_radius_z
spes_table[spes_table$variable_parameter == "E_radius_z", "variable_parameter"] <- "E_volume" 

# Set up plots metadata
non_gradient_plots_metadata <- list(
  ellipsoid = list(
    arrangement = list(x_aes = "temp_arrangement", y_aes = "metric"),
    bg_prop_A = list(x_aes = "bg_prop_A", y_aes = "metric"),
    bg_prop_B = list(x_aes = "bg_prop_B", y_aes = "metric"),
    E_volume = list(x_aes = "E_volume", y_aes = "metric")
  ),
  network = list(
    arrangement = list(x_aes = "temp_arrangement", y_aes = "metric"),
    bg_prop_A = list(x_aes = "bg_prop_A", y_aes = "metric"),
    bg_prop_B = list(x_aes = "bg_prop_B", y_aes = "metric"),
    N_width = list(x_aes = "N_width", y_aes = "metric")
  )
)

# Generate plots and plots into a list
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "CGR_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")

background_parameters <- c("bg_prop_A", "bg_prop_B")

arrangement_parameters <- list(mixed = "cluster_prop_B",
                               ringed = "ring_width_factor",
                               separated = "distance")

shape_parameters <- list(ellipsoid = c("E_volume"),
                         network = c("N_width"))


metric_plots3D <- list(mixed_ellipsoid = list(),
                       mixed_network = list(),
                       ringed_ellipsoid = list(),
                       ringed_network = list(),
                       separated_ellipsoid = list(),
                       separated_network = list())

for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    spes_table_subset <- spes_table[spes_table$variable_parameter %in% c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]]), 
                                    c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]], "variable_parameter")]
    
    for (metric in metrics) {
      metric_plots3D[[spes_metadata_index]][[metric]] <- plot_non_gradient_metric(spes_table_subset, 
                                                                                  metric, 
                                                                                  metric_df_lists3D[[spes_metadata_index]][[metric]], 
                                                                                  arrangement_parameters[[arrangement]], 
                                                                                  non_gradient_plots_metadata[[shape]])
      
    }
  }
}

# Put plots into a pdf
setwd("~/R/plots/S1")
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")

metrics_set1 <- c("AMD",  "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "CGR_AUC")
metrics_set2 <- c("MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")

pdf("plots3D_non_gradient.pdf", width = 25, height = 12)

for (metric in metrics_set1) {
  for (shape in shapes) {
    curr_metric_plots <- list()
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[arrangement]] <- metric_plots3D[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
    }
    plot <- plot_grid(plotlist = curr_metric_plots,
                      nrow = 1, 
                      ncol = length(arrangements))
    print(plot)
  }
}

for (metric in metrics_set2) {
  curr_metric_plots <- list()
  for (shape in shapes) {
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[spes_metadata_index]] <- metric_plots3D[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
    }
  }
  plot <- plot_grid(plotlist = curr_metric_plots,
                    nrow = length(shapes), 
                    ncol = length(arrangements))
  print(plot)
}


dev.off()

### Get plots for 3D metric analysis (gradient) ------------------------------------------
# Read spes_table
setwd("~/R/SPIAT-3D_benchmarking/simulations_and_analysis_S1/S1_data")
spes_table <- read.table("spes_table.csv")
spes_table$cluster_prop_B <- 1 - spes_table$cluster_prop_A
spes_table[spes_table$variable_parameter == "cluster_prop_A", "variable_parameter"] <- "cluster_prop_B"
spes_table$distance <- 450 - spes_table$cluster_x_coord
spes_table[spes_table$variable_parameter == "cluster_x_coord", "variable_parameter"] <- "distance" 
spes_table$E_volume <- (4/3) * pi * spes_table$E_radius_x * spes_table$E_radius_y * spes_table$E_radius_z
spes_table[spes_table$variable_parameter == "E_radius_z", "variable_parameter"] <- "E_volume" 

gradient_plots_metadata <- list(
  ellipsoid = list(
    arrangement = list(x_aes = "gradient", y_aes = "metric", color_aes = "temp_arrangement"),
    bg_prop_A = list(x_aes = "gradient", y_aes = "metric", color_aes = "bg_prop_A"),
    bg_prop_B = list(x_aes = "gradient", y_aes = "metric", color_aes = "bg_prop_B"),
    E_volume = list(x_aes = "gradient", y_aes = "metric", color_aes = "E_volume")
  ),
  network = list(
    arrangement = list(x_aes = "gradient", y_aes = "metric", color_aes = "temp_arrangement"),
    bg_prop_A = list(x_aes = "gradient", y_aes = "metric", color_aes = "bg_prop_A"),
    bg_prop_B = list(x_aes = "gradient", y_aes = "metric", color_aes = "bg_prop_B"),
    N_width = list(x_aes = "gradient", y_aes = "metric", color_aes = "N_width")
  )
)

# Generate plots and plots into a list
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics <- c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "CLR", "COO", "CGR", "PBP", "EBP")

background_parameters <- c("bg_prop_A", "bg_prop_B")

arrangement_parameters <- list(mixed = "cluster_prop_B",
                               ringed = "ring_width_factor",
                               separated = "distance")

shape_parameters <- list(ellipsoid = c("E_volume"),
                         network = c("N_width"))


metric_plots3D <- list(mixed_ellipsoid = list(),
                       mixed_network = list(),
                       ringed_ellipsoid = list(),
                       ringed_network = list(),
                       separated_ellipsoid = list(),
                       separated_network = list())

for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    spes_table_subset <- spes_table[spes_table$variable_parameter %in% c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]]), 
                                    c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]], "variable_parameter")]
    
    for (metric in metrics) {
      
      metric_plots3D[[spes_metadata_index]][[metric]] <- plot_gradient_metric(spes_table_subset, 
                                                                              metric,
                                                                              metric_df_lists3D[[spes_metadata_index]][[metric]], 
                                                                              arrangement_parameters[[arrangement]], 
                                                                              get_gradient(metric),
                                                                              gradient_plots_metadata[[shape]])
      
    }
  }
}

# Put plots into a pdf
setwd("~/R/plots/S1")
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")

metrics_set1 <- c("ACIN", "CKR", "CLR", "COO", "CGR")
metrics_set2 <- c("MS", "NMS", "ACINP", "AE", "PBP", "EBP")

pdf("plots3D_gradient.pdf", width = 25, height = 12)

for (metric in metrics_set1) {
  for (shape in shapes) {
    curr_metric_plots <- list()
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[arrangement]] <- metric_plots3D[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
    }
    plot <- plot_grid(plotlist = curr_metric_plots,
                      nrow = 1, 
                      ncol = length(arrangements))
    print(plot)
  }
}

for (metric in metrics_set2) {
  curr_metric_plots <- list()
  for (shape in shapes) {
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[spes_metadata_index]] <- metric_plots3D[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
    }
  }
  plot <- plot_grid(plotlist = curr_metric_plots,
                    nrow = length(shapes), 
                    ncol = length(arrangements))
  print(plot)
}


dev.off()

### Get plots for 2D metric analysis (middle slice only) ------------------------------------------
# Read spes_table
setwd("~/R/SPIAT-3D_benchmarking/simulations_and_analysis_S1/S1_data")
spes_table <- read.table("spes_table.csv")
spes_table$cluster_prop_B <- 1 - spes_table$cluster_prop_A
spes_table[spes_table$variable_parameter == "cluster_prop_A", "variable_parameter"] <- "cluster_prop_B"
spes_table$distance <- 450 - spes_table$cluster_x_coord
spes_table[spes_table$variable_parameter == "cluster_x_coord", "variable_parameter"] <- "distance" 
spes_table$E_volume <- (4/3) * pi * spes_table$E_radius_x * spes_table$E_radius_y * spes_table$E_radius_z
spes_table[spes_table$variable_parameter == "E_radius_z", "variable_parameter"] <- "E_volume" 

# Subset metric_df_lists2D to only include the middle slice
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metric_df_lists2D_subset <- metric_df_lists2D
for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    curr_list <- metric_df_lists2D_subset[[spes_metadata_index]]
    
    for (i in seq(length(curr_list))) {
      curr_df <- curr_list[[i]]
      curr_df <- curr_df[curr_df$slice == 1, ]
      curr_list[[i]] <- curr_df
    }
    metric_df_lists2D_subset[[spes_metadata_index]] <- curr_list
  }
}

# Set up plots metadata
non_gradient_plots_metadata <- list(
  ellipsoid = list(
    arrangement = list(x_aes = "temp_arrangement", y_aes = "metric"),
    bg_prop_A = list(x_aes = "bg_prop_A", y_aes = "metric"),
    bg_prop_B = list(x_aes = "bg_prop_B", y_aes = "metric"),
    E_volume = list(x_aes = "E_volume", y_aes = "metric")
  ),
  network = list(
    arrangement = list(x_aes = "temp_arrangement", y_aes = "metric"),
    bg_prop_A = list(x_aes = "bg_prop_A", y_aes = "metric"),
    bg_prop_B = list(x_aes = "bg_prop_B", y_aes = "metric"),
    N_width = list(x_aes = "N_width", y_aes = "metric")
  )
)

gradient_plots_metadata <- list(
  ellipsoid = list(
    arrangement = list(x_aes = "gradient", y_aes = "metric", color_aes = "temp_arrangement"),
    bg_prop_A = list(x_aes = "gradient", y_aes = "metric", color_aes = "bg_prop_A"),
    bg_prop_B = list(x_aes = "gradient", y_aes = "metric", color_aes = "bg_prop_B"),
    E_volume = list(x_aes = "gradient", y_aes = "metric", color_aes = "E_volume")
  ),
  network = list(
    arrangement = list(x_aes = "gradient", y_aes = "metric", color_aes = "temp_arrangement"),
    bg_prop_A = list(x_aes = "gradient", y_aes = "metric", color_aes = "bg_prop_A"),
    bg_prop_B = list(x_aes = "gradient", y_aes = "metric", color_aes = "bg_prop_B"),
    N_width = list(x_aes = "gradient", y_aes = "metric", color_aes = "N_width")
  )
)

# Generate plots and plots into a list
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics <- c("AMD", "MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "CLR", "COO", "CGR", "PBSAC", "PBP", "PBP_AUC", "EBSAC", "EBP", "EBP_AUC")

background_parameters <- c("bg_prop_A", "bg_prop_B")

arrangement_parameters <- list(mixed = "cluster_prop_B",
                               ringed = "ring_width_factor",
                               separated = "distance")

shape_parameters <- list(ellipsoid = c("E_volume"),
                         network = c("N_width"))


metric_plots2D <- list(mixed_ellipsoid = list(),
                       mixed_network = list(),
                       ringed_ellipsoid = list(),
                       ringed_network = list(),
                       separated_ellipsoid = list(),
                       separated_network = list())

for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    spes_table_subset <- spes_table[spes_table$variable_parameter %in% c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]]), 
                                    c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]], "variable_parameter")]
    
    for (metric in metrics) {
      if (metric %in% c("AMD", "PBSAC", "EBSAC", "PBP_AUC", "EBP_AUC")) {
        metric_plots2D[[spes_metadata_index]][[metric]] <- plot_non_gradient_metric(spes_table_subset, 
                                                                                    metric, 
                                                                                    metric_df_lists2D_subset[[spes_metadata_index]][[metric]], 
                                                                                    arrangement_parameters[[arrangement]], 
                                                                                    non_gradient_plots_metadata[[shape]])
      }
      else if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "CLR", "COO", "CGR", "PBP", "EBP")) {
        metric_plots2D[[spes_metadata_index]][[metric]] <- plot_gradient_metric(spes_table_subset, 
                                                                                metric,
                                                                                metric_df_lists2D_subset[[spes_metadata_index]][[metric]], 
                                                                                arrangement_parameters[[arrangement]], 
                                                                                get_gradient(metric),
                                                                                gradient_plots_metadata[[shape]])
      }
    }
  }
}

# Put plots into a pdf
setwd("~/R/plots/S1")
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")

metrics_set1 <- c("AMD", "ACIN", "CKR", "CLR", "COO", "CGR")
metrics_set2 <- c("MS", "NMS", "ACINP", "AE", "PBSAC", "PBP", "PBP_AUC", "EBSAC", "EBP", "EBP_AUC")

pdf("plots2D_middle_slice.pdf", width = 25, height = 12)

for (metric in metrics_set1) {
  for (shape in shapes) {
    curr_metric_plots <- list()
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[arrangement]] <- metric_plots2D[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
    }
    plot <- plot_grid(plotlist = curr_metric_plots,
                      nrow = 1, 
                      ncol = length(arrangements))
    print(plot)
  }
}

for (metric in metrics_set2) {
  curr_metric_plots <- list()
  for (shape in shapes) {
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[spes_metadata_index]] <- metric_plots2D[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
    }
  }
  plot <- plot_grid(plotlist = curr_metric_plots,
                    nrow = length(shapes), 
                    ncol = length(arrangements))
  print(plot)
}

dev.off()

### Get plots with 2D (middle slice only) on the x-axis and 3D on the y-axis ----------------
setwd("~/R/SPIAT-3D_benchmarking/simulations_and_analysis_S1/S1_data")
spes_table <- read.table("spes_table.csv")
spes_table$cluster_prop_B <- 1 - spes_table$cluster_prop_A
spes_table[spes_table$variable_parameter == "cluster_prop_A", "variable_parameter"] <- "cluster_prop_B" 
spes_table$distance <- 450 - spes_table$cluster_x_coord
spes_table[spes_table$variable_parameter == "cluster_x_coord", "variable_parameter"] <- "distance" 
spes_table$E_volume <- (4/3) * pi * spes_table$E_radius_x * spes_table$E_radius_y * spes_table$E_radius_z
spes_table[spes_table$variable_parameter == "E_radius_z", "variable_parameter"] <- "E_volume" 



# Subset metric_df_lists2D to only include the middle slice
metric_df_lists2D_subset <- metric_df_lists2D
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    curr_list <- metric_df_lists2D_subset[[spes_metadata_index]]
    
    for (i in seq(length(curr_list))) {
      curr_df <- curr_list[[i]]
      curr_df <- curr_df[curr_df$slice == 1, ]
      curr_list[[i]] <- curr_df
    }
    metric_df_lists2D_subset[[spes_metadata_index]] <- curr_list
  }
}




# Set up plots metadata
plots_metadata <- list(
  ellipsoid = list(
    arrangement = list(x_aes = "2D", y_aes = "3D", color_aes = "temp_arrangement"),
    bg_prop_A = list(x_aes = "2D", y_aes = "3D", color_aes = "bg_prop_A"),
    bg_prop_B = list(x_aes = "2D", y_aes = "3D", color_aes = "bg_prop_B"),
    E_volume = list(x_aes = "2D", y_aes = "3D", color_aes = "E_volume")
  ),
  network = list(
    arrangement = list(x_aes = "2D", y_aes = "3D", color_aes = "temp_arrangement"),
    bg_prop_A = list(x_aes = "2D", y_aes = "3D", color_aes = "bg_prop_A"),
    bg_prop_B = list(x_aes = "2D", y_aes = "3D", color_aes = "bg_prop_B"),
    N_width = list(x_aes = "2D", y_aes = "3D", color_aes = "N_width")
  )
)



# Generate plots and plots into a list
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")

background_parameters <- c("bg_prop_A", "bg_prop_B")

arrangement_parameters <- list(mixed = "cluster_prop_B",
                               ringed = "ring_width_factor",
                               separated = "distance")

shape_parameters <- list(ellipsoid = c("E_volume"),
                         network = c("N_width"))

metric_plots2D_vs_3D_middle_slice <- list(mixed_ellipsoid = list(),
                                          mixed_network = list(),
                                          ringed_ellipsoid = list(),
                                          ringed_network = list(),
                                          separated_ellipsoid = list(),
                                          separated_network = list())

for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    spes_table_subset <- spes_table[spes_table$variable_parameter %in% c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]]), 
                                    c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]], "variable_parameter")]
    
    for (metric in metrics) {
      metric_plots2D_vs_3D_middle_slice[[spes_metadata_index]][[metric]] <- plot_3D_vs_2D_metric_one_slice(spes_table_subset, 
                                                                                                           metric, 
                                                                                                           metric_df_lists3D[[spes_metadata_index]][[metric]],
                                                                                                           metric_df_lists2D_subset[[spes_metadata_index]][[metric]], 
                                                                                                           arrangement_parameters[[arrangement]], 
                                                                                                           plots_metadata[[shape]])
    }
  }
}


# Put plots into a pdf
setwd("~/R/plots/S1")
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics_set1 <- c("AMD",  "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC")
metrics_set2 <- c("MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")

pdf("plots2D_vs_3D_middle_slice.pdf", width = 25, height = 10)

for (metric in metrics_set1) {
  for (shape in shapes) {
    curr_metric_plots <- list()
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[arrangement]] <- metric_plots2D_vs_3D_middle_slice[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
    }
    plot <- plot_grid(plotlist = curr_metric_plots,
                      nrow = 1, 
                      ncol = length(arrangements))
    print(plot)
  }
}

for (metric in metrics_set2) {
  curr_metric_plots <- list()
  for (shape in shapes) {
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[spes_metadata_index]] <- metric_plots2D_vs_3D_middle_slice[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
    }
  }
  plot <- plot_grid(plotlist = curr_metric_plots,
                    nrow = length(shapes), 
                    ncol = length(arrangements))
  print(plot)
}

dev.off()






### Get plots with 2D (all slices) on the x-axis and 3D on the y-axis ----------------
setwd("~/R/SPIAT-3D_benchmarking/simulations_and_analysis_S1/S1_data")
spes_table <- read.table("spes_table.csv")
spes_table$cluster_prop_B <- 1 - spes_table$cluster_prop_A
spes_table[spes_table$variable_parameter == "cluster_prop_A", "variable_parameter"] <- "cluster_prop_B" 
spes_table$distance <- 450 - spes_table$cluster_x_coord
spes_table[spes_table$variable_parameter == "cluster_x_coord", "variable_parameter"] <- "distance" 
spes_table$E_volume <- (4/3) * pi * spes_table$E_radius_x * spes_table$E_radius_y * spes_table$E_radius_z
spes_table[spes_table$variable_parameter == "E_radius_z", "variable_parameter"] <- "E_volume" 

# Set up plots metadata
plots_metadata <- list(
  ellipsoid = list(
    arrangement = list(x_aes = "2D", y_aes = "3D", color_aes = "slice", label = "temp_arrangement"),
    bg_prop_A = list(x_aes = "2D", y_aes = "3D", color_aes = "slice", label = "bg_prop_A"),
    bg_prop_B = list(x_aes = "2D", y_aes = "3D", color_aes = "slice", label = "bg_prop_B"),
    E_volume = list(x_aes = "2D", y_aes = "3D", color_aes = "slice", label = "E_volume")
  ),
  network = list(
    arrangement = list(x_aes = "2D", y_aes = "3D", color_aes = "slice", label = "temp_arrangement"),
    bg_prop_A = list(x_aes = "2D", y_aes = "3D", color_aes = "slice", label = "bg_prop_A"),
    bg_prop_B = list(x_aes = "2D", y_aes = "3D", color_aes = "slice", label = "bg_prop_B"),
    N_width = list(x_aes = "2D", y_aes = "3D", color_aes = "slice", label = "N_width")
  )
)



# Generate plots and plots into a list
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")

background_parameters <- c("bg_prop_A", "bg_prop_B")

arrangement_parameters <- list(mixed = "cluster_prop_B",
                               ringed = "ring_width_factor",
                               separated = "distance")

shape_parameters <- list(ellipsoid = c("E_volume"),
                         network = c("N_width"))

metric_plots2D_vs_3D_all_slices <- list(mixed_ellipsoid = list(),
                                        mixed_network = list(),
                                        ringed_ellipsoid = list(),
                                        ringed_network = list(),
                                        separated_ellipsoid = list(),
                                        separated_network = list())

for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    spes_table_subset <- spes_table[spes_table$variable_parameter %in% c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]]), 
                                    c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]], "variable_parameter")]
    
    for (metric in metrics) {
      metric_plots2D_vs_3D_all_slices[[spes_metadata_index]][[metric]] <- plot_3D_vs_2D_metric_all_slices(spes_table_subset, 
                                                                                                          metric, 
                                                                                                          metric_df_lists3D[[spes_metadata_index]][[metric]],
                                                                                                          metric_df_lists2D[[spes_metadata_index]][[metric]], 
                                                                                                          arrangement_parameters[[arrangement]], 
                                                                                                          plots_metadata[[shape]])
    }
  }
}



# Put plots into a pdf
setwd("~/R/plots/S1")
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics_set1 <- c("AMD",  "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC")
metrics_set2 <- c("MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")

pdf("plots2D_vs_3D_all_slices.pdf", width = 25, height = 10)

for (metric in metrics_set1) {
  for (shape in shapes) {
    curr_metric_plots <- list()
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[arrangement]] <- metric_plots2D_vs_3D_all_slices[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
    }
    plot <- plot_grid(plotlist = curr_metric_plots,
                      nrow = 1, 
                      ncol = length(arrangements))
    print(plot)
  }
}

for (metric in metrics_set2) {
  curr_metric_plots <- list()
  for (shape in shapes) {
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[spes_metadata_index]] <- metric_plots2D_vs_3D_all_slices[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
    }
  }
  plot <- plot_grid(plotlist = curr_metric_plots,
                    nrow = length(shapes), 
                    ncol = length(arrangements))
  print(plot)
}

dev.off()






### Get plots for error for non-gradient metrics ----------------------------------
setwd("~/R/SPIAT-3D_benchmarking/simulations_and_analysis_S1/S1_data")
spes_table <- read.table("spes_table.csv")
spes_table$cluster_prop_B <- 1 - spes_table$cluster_prop_A
spes_table[spes_table$variable_parameter == "cluster_prop_A", "variable_parameter"] <- "cluster_prop_B" 
spes_table$distance <- 450 - spes_table$cluster_x_coord
spes_table[spes_table$variable_parameter == "cluster_x_coord", "variable_parameter"] <- "distance" 
spes_table$E_volume <- (4/3) * pi * spes_table$E_radius_x * spes_table$E_radius_y * spes_table$E_radius_z
spes_table[spes_table$variable_parameter == "E_radius_z", "variable_parameter"] <- "E_volume" 


# Set up plots metadata
plots_metadata <- list(
  ellipsoid = list(
    arrangement = list(x_aes = "temp_arrangement", y_aes = "error", color_aes = "slice"),
    bg_prop_A = list(x_aes = "bg_prop_A", y_aes = "error", color_aes = "slice"),
    bg_prop_B = list(x_aes = "bg_prop_B", y_aes = "error", color_aes = "slice"),
    E_volume = list(x_aes = "E_volume", y_aes = "error", color_aes = "slice")
  ),
  network = list(
    arrangement = list(x_aes = "temp_arrangement", y_aes = "error", color_aes = "slice"),
    bg_prop_A = list(x_aes = "bg_prop_A", y_aes = "error", color_aes = "slice"),
    bg_prop_B = list(x_aes = "bg_prop_B", y_aes = "error", color_aes = "slice"),
    N_width = list(x_aes = "N_width", y_aes = "error", color_aes = "slice")
  )
)



# Generate plots and plots into a list
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")

background_parameters <- c("bg_prop_A", "bg_prop_B")

arrangement_parameters <- list(mixed = "cluster_prop_B",
                               ringed = "ring_width_factor",
                               separated = "distance")

shape_parameters <- list(ellipsoid = c("E_volume"),
                         network = c("N_width"))

metric_plots_error_non_gradient <- list(mixed_ellipsoid = list(),
                                        mixed_network = list(),
                                        ringed_ellipsoid = list(),
                                        ringed_network = list(),
                                        separated_ellipsoid = list(),
                                        separated_network = list())

for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    spes_table_subset <- spes_table[spes_table$variable_parameter %in% c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]]), 
                                    c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]], "variable_parameter")]
    
    for (metric in metrics) {
      metric_plots_error_non_gradient[[spes_metadata_index]][[metric]] <- plot_error_non_gradient_metric(spes_table_subset, 
                                                                                                         metric, 
                                                                                                         metric_df_lists3D[[spes_metadata_index]][[metric]],
                                                                                                         metric_df_lists2D[[spes_metadata_index]][[metric]], 
                                                                                                         arrangement_parameters[[arrangement]], 
                                                                                                         plots_metadata[[shape]])
    }
  }
}


# Put plots into a pdf
setwd("~/R/plots/S1")
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics_set1 <- c("AMD",  "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC")
metrics_set2 <- c("MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")

pdf("plots_error_non_gradient_all_slices.pdf", width = 25, height = 10)

for (metric in metrics_set1) {
  for (shape in shapes) {
    curr_metric_plots <- list()
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[arrangement]] <- metric_plots_error_non_gradient[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
    }
    plot <- plot_grid(plotlist = curr_metric_plots,
                      nrow = 1, 
                      ncol = length(arrangements))
    print(plot)
  }
}

for (metric in metrics_set2) {
  curr_metric_plots <- list()
  for (shape in shapes) {
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[spes_metadata_index]] <- metric_plots_error_non_gradient[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
    }
  }
  plot <- plot_grid(plotlist = curr_metric_plots,
                    nrow = length(shapes), 
                    ncol = length(arrangements))
  print(plot)
}

dev.off()





### Get plots for error for gradient metrics one slice ----------------------------------
# Read spes_table
setwd("~/R/SPIAT-3D_benchmarking/simulations_and_analysis_S1/S1_data")
spes_table <- read.table("spes_table.csv")
spes_table$cluster_prop_B <- 1 - spes_table$cluster_prop_A
spes_table[spes_table$variable_parameter == "cluster_prop_A", "variable_parameter"] <- "cluster_prop_B"
spes_table$distance <- 450 - spes_table$cluster_x_coord
spes_table[spes_table$variable_parameter == "cluster_x_coord", "variable_parameter"] <- "distance" 
spes_table$E_volume <- (4/3) * pi * spes_table$E_radius_x * spes_table$E_radius_y * spes_table$E_radius_z
spes_table[spes_table$variable_parameter == "E_radius_z", "variable_parameter"] <- "E_volume" 

# Subset metric_df_lists2D to only include the middle slice
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metric_df_lists2D_subset <- metric_df_lists2D
for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    curr_list <- metric_df_lists2D_subset[[spes_metadata_index]]
    
    for (i in seq(length(curr_list))) {
      curr_df <- curr_list[[i]]
      curr_df <- curr_df[curr_df$slice == 1, ]
      curr_list[[i]] <- curr_df
    }
    metric_df_lists2D_subset[[spes_metadata_index]] <- curr_list
  }
}

gradient_plots_metadata <- list(
  ellipsoid = list(
    arrangement = list(x_aes = "gradient", y_aes = "error", color_aes = "temp_arrangement"),
    bg_prop_A = list(x_aes = "gradient", y_aes = "error", color_aes = "bg_prop_A"),
    bg_prop_B = list(x_aes = "gradient", y_aes = "error", color_aes = "bg_prop_B"),
    E_volume = list(x_aes = "gradient", y_aes = "error", color_aes = "E_volume")
  ),
  network = list(
    arrangement = list(x_aes = "gradient", y_aes = "error", color_aes = "temp_arrangement"),
    bg_prop_A = list(x_aes = "gradient", y_aes = "error", color_aes = "bg_prop_A"),
    bg_prop_B = list(x_aes = "gradient", y_aes = "error", color_aes = "bg_prop_B"),
    N_width = list(x_aes = "gradient", y_aes = "error", color_aes = "N_width")
  )
)


# Generate plots and plots into a list
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics <- c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "CLR", "COO", "CGR", "PBP",  "EBP")

background_parameters <- c("bg_prop_A", "bg_prop_B")

arrangement_parameters <- list(mixed = "cluster_prop_B",
                               ringed = "ring_width_factor",
                               separated = "distance")

shape_parameters <- list(ellipsoid = c("E_volume"),
                         network = c("N_width"))

metric_plots_error_gradient <- list(mixed_ellipsoid = list(),
                                    mixed_network = list(),
                                    ringed_ellipsoid = list(),
                                    ringed_network = list(),
                                    separated_ellipsoid = list(),
                                    separated_network = list())

for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    spes_table_subset <- spes_table[spes_table$variable_parameter %in% c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]]), 
                                    c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]], "variable_parameter")]
    
    for (metric in metrics) {
      metric_plots_error_gradient[[spes_metadata_index]][[metric]] <- plot_error_gradient_metric_one_slice(spes_table_subset, 
                                                                                                           metric, 
                                                                                                           metric_df_lists3D[[spes_metadata_index]][[metric]],
                                                                                                           metric_df_lists2D_subset[[spes_metadata_index]][[metric]], 
                                                                                                           arrangement_parameters[[arrangement]],
                                                                                                           get_gradient(metric),
                                                                                                           gradient_plots_metadata[[shape]])
    }
  }
}


# Put plots into a pdf
setwd("~/R/plots/S1")
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics_set1 <- c("ACIN", "CKR", "CLR", "COO", "CGR")
metrics_set2 <- c("MS", "NMS", "ACINP", "AE", "PBP", "EBP")

pdf("plots_error_gradient_one_slice.pdf", width = 25, height = 10)

for (metric in metrics_set1) {
  for (shape in shapes) {
    curr_metric_plots <- list()
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[arrangement]] <- metric_plots_error_gradient[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
    }
    plot <- plot_grid(plotlist = curr_metric_plots,
                      nrow = 1, 
                      ncol = length(arrangements))
    print(plot)
  }
}

for (metric in metrics_set2) {
  curr_metric_plots <- list()
  for (shape in shapes) {
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[spes_metadata_index]] <- metric_plots_error_gradient[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
    }
  }
  plot <- plot_grid(plotlist = curr_metric_plots,
                    nrow = length(shapes), 
                    ncol = length(arrangements))
  print(plot)
}

dev.off()

### Get plots for 2D metric analysis (all slices with ground truth) -----
setwd("~/R/SPIAT-3D_benchmarking/simulations_and_analysis_S1/S1_data")
spes_table <- read.table("spes_table.csv")
spes_table$cluster_prop_B <- 1 - spes_table$cluster_prop_A
spes_table[spes_table$variable_parameter == "cluster_prop_A", "variable_parameter"] <- "cluster_prop_B" 
spes_table$distance <- 450 - spes_table$cluster_x_coord
spes_table[spes_table$variable_parameter == "cluster_x_coord", "variable_parameter"] <- "distance" 
spes_table$E_volume <- (4/3) * pi * spes_table$E_radius_x * spes_table$E_radius_y * spes_table$E_radius_z
spes_table[spes_table$variable_parameter == "E_radius_z", "variable_parameter"] <- "E_volume" 




# Set up plots metadata
plots_metadata <- list(
  ellipsoid = list(
    arrangement = list(x_aes = "temp_arrangement", y_aes = "2D", color_aes = "slice"),
    bg_prop_A = list(x_aes = "bg_prop_A", y_aes = "2D", color_aes = "slice"),
    bg_prop_B = list(x_aes = "bg_prop_B", y_aes = "2D", color_aes = "slice"),
    E_volume = list(x_aes = "E_volume", y_aes = "2D", color_aes = "slice")
  ),
  network = list(
    arrangement = list(x_aes = "temp_arrangement", y_aes = "2D", color_aes = "slice"),
    bg_prop_A = list(x_aes = "bg_prop_A", y_aes = "2D", color_aes = "slice"),
    bg_prop_B = list(x_aes = "bg_prop_B", y_aes = "2D", color_aes = "slice"),
    N_width = list(x_aes = "N_width", y_aes = "2D", color_aes = "slice")
  )
)



# Generate plots and plots into a list
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")

background_parameters <- c("bg_prop_A", "bg_prop_B")

arrangement_parameters <- list(mixed = "cluster_prop_B",
                               ringed = "ring_width_factor",
                               separated = "distance")

shape_parameters <- list(ellipsoid = c("E_volume"),
                         network = c("N_width"))

metric_plots2D_all_slices_ground_truth <- list(mixed_ellipsoid = list(),
                                               mixed_network = list(),
                                               ringed_ellipsoid = list(),
                                               ringed_network = list(),
                                               separated_ellipsoid = list(),
                                               separated_network = list())

for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    spes_table_subset <- spes_table[spes_table$variable_parameter %in% c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]]), 
                                    c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]], "variable_parameter")]
    
    for (metric in metrics) {
      metric_plots2D_all_slices_ground_truth[[spes_metadata_index]][[metric]] <- plot_non_gradient_metric_all_slices_ground_truth(spes_table_subset, 
                                                                                                                                  metric, 
                                                                                                                                  metric_df_lists3D[[spes_metadata_index]][[metric]], 
                                                                                                                                  metric_df_lists2D[[spes_metadata_index]][[metric]], 
                                                                                                                                  arrangement_parameters[[arrangement]], 
                                                                                                                                  plots_metadata[[shape]])
    }
  }
}


# Put plots into a pdf
setwd("~/R/plots/S1")
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics_set1 <- c("AMD", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC")
metrics_set2 <- c("MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")

pdf("plots2D_all_slices_with_ground_truth.pdf", width = 25, height = 10)

for (metric in metrics_set1) {
  for (shape in shapes) {
    curr_metric_plots <- list()
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[arrangement]] <- metric_plots2D_all_slices_ground_truth[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
    }
    plot <- plot_grid(plotlist = curr_metric_plots,
                      nrow = 1, 
                      ncol = length(arrangements))
    print(plot)
  }
}

for (metric in metrics_set2) {
  curr_metric_plots <- list()
  for (shape in shapes) {
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[spes_metadata_index]] <- metric_plots2D_all_slices_ground_truth[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
    }
  }
  plot <- plot_grid(plotlist = curr_metric_plots,
                    nrow = length(shapes), 
                    ncol = length(arrangements))
  print(plot)
}

dev.off()




### Get plots for violin plots for slice (all slices) -----
setwd("~/R/SPIAT-3D_benchmarking/simulations_and_analysis_S1/S1_data")
spes_table <- read.table("spes_table.csv")
spes_table$cluster_prop_B <- 1 - spes_table$cluster_prop_A
spes_table[spes_table$variable_parameter == "cluster_prop_A", "variable_parameter"] <- "cluster_prop_B" 
spes_table$distance <- 450 - spes_table$cluster_x_coord
spes_table[spes_table$variable_parameter == "cluster_x_coord", "variable_parameter"] <- "distance" 
spes_table$E_volume <- (4/3) * pi * spes_table$E_radius_x * spes_table$E_radius_y * spes_table$E_radius_z
spes_table[spes_table$variable_parameter == "E_radius_z", "variable_parameter"] <- "E_volume" 




# Set up plots metadata
plots_metadata <- list(
  ellipsoid = list(
    arrangement = list(x_aes = "slice", y_aes = "metric", label = "temp_arrangement"),
    bg_prop_A = list(x_aes = "slice", y_aes = "metric", label = "bg_prop_A"),
    bg_prop_B = list(x_aes = "slice", y_aes = "metric", label = "bg_prop_B"),
    E_volume = list(x_aes = "slice", y_aes = "metric", label = "E_volume")
  ),
  network = list(
    arrangement = list(x_aes = "slice", y_aes = "metric", label = "temp_arrangement"),
    bg_prop_A = list(x_aes = "slice", y_aes = "metric", label = "bg_prop_A"),
    bg_prop_B = list(x_aes = "slice", y_aes = "metric", label = "bg_prop_B"),
    N_width = list(x_aes = "slice", y_aes = "metric", label = "N_width")
  )
)



# Generate plots and plots into a list
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")

background_parameters <- c("bg_prop_A", "bg_prop_B")

arrangement_parameters <- list(mixed = "cluster_prop_B",
                               ringed = "ring_width_factor",
                               separated = "distance")

shape_parameters <- list(ellipsoid = c("E_volume"),
                         network = c("N_width"))

metric_plots_violin_all_slices <- list(mixed_ellipsoid = list(),
                                       mixed_network = list(),
                                       ringed_ellipsoid = list(),
                                       ringed_network = list(),
                                       separated_ellipsoid = list(),
                                       separated_network = list())

for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    spes_table_subset <- spes_table[spes_table$variable_parameter %in% c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]]), 
                                    c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]], "variable_parameter")]
    
    for (metric in metrics) {
      metric_plots_violin_all_slices[[spes_metadata_index]][[metric]] <- plot_violin_all_slices(spes_table_subset, 
                                                                                                metric, 
                                                                                                metric_df_lists2D[[spes_metadata_index]][[metric]], 
                                                                                                arrangement_parameters[[arrangement]], 
                                                                                                plots_metadata[[shape]])
    }
  }
}


# Put plots into a pdf
setwd("~/R/plots/S1")
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics_set1 <- c("AMD", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC")
metrics_set2 <- c("MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")

pdf("metric_plots_violin_all_slices.pdf", width = 25, height = 10)

for (metric in metrics_set1) {
  for (shape in shapes) {
    curr_metric_plots <- list()
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[arrangement]] <- metric_plots_violin_all_slices[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
    }
    plot <- plot_grid(plotlist = curr_metric_plots,
                      nrow = 1, 
                      ncol = length(arrangements))
    print(plot)
  }
}

for (metric in metrics_set2) {
  curr_metric_plots <- list()
  for (shape in shapes) {
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[spes_metadata_index]] <- metric_plots_violin_all_slices[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
    }
  }
  plot <- plot_grid(plotlist = curr_metric_plots,
                    nrow = length(shapes), 
                    ncol = length(arrangements))
  print(plot)
}

dev.off()




### Get plots for violin plots for slice (all slices with ground truth) -----
setwd("~/R/SPIAT-3D_benchmarking/simulations_and_analysis_S1/S1_data")
spes_table <- read.table("spes_table.csv")
spes_table$cluster_prop_B <- 1 - spes_table$cluster_prop_A
spes_table[spes_table$variable_parameter == "cluster_prop_A", "variable_parameter"] <- "cluster_prop_B" 
spes_table$distance <- 450 - spes_table$cluster_x_coord
spes_table[spes_table$variable_parameter == "cluster_x_coord", "variable_parameter"] <- "distance" 
spes_table$E_volume <- (4/3) * pi * spes_table$E_radius_x * spes_table$E_radius_y * spes_table$E_radius_z
spes_table[spes_table$variable_parameter == "E_radius_z", "variable_parameter"] <- "E_volume" 




# Set up plots metadata
plots_metadata <- list(
  ellipsoid = list(
    arrangement = list(x_aes = "slice", y_aes = "metric", label = "temp_arrangement"),
    bg_prop_A = list(x_aes = "slice", y_aes = "metric", label = "bg_prop_A"),
    bg_prop_B = list(x_aes = "slice", y_aes = "metric", label = "bg_prop_B"),
    E_volume = list(x_aes = "slice", y_aes = "metric", label = "E_volume")
  ),
  network = list(
    arrangement = list(x_aes = "slice", y_aes = "metric", label = "temp_arrangement"),
    bg_prop_A = list(x_aes = "slice", y_aes = "metric", label = "bg_prop_A"),
    bg_prop_B = list(x_aes = "slice", y_aes = "metric", label = "bg_prop_B"),
    N_width = list(x_aes = "slice", y_aes = "metric", label = "N_width")
  )
)



# Generate plots and plots into a list
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")

background_parameters <- c("bg_prop_A", "bg_prop_B")

arrangement_parameters <- list(mixed = "cluster_prop_B",
                               ringed = "ring_width_factor",
                               separated = "distance")

shape_parameters <- list(ellipsoid = c("E_volume"),
                         network = c("N_width"))

metric_plots_violin_all_slices_ground_truth <- list(mixed_ellipsoid = list(),
                                                    mixed_network = list(),
                                                    ringed_ellipsoid = list(),
                                                    ringed_network = list(),
                                                    separated_ellipsoid = list(),
                                                    separated_network = list())

for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    spes_table_subset <- spes_table[spes_table$variable_parameter %in% c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]]), 
                                    c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]], "variable_parameter")]
    
    for (metric in metrics) {
      metric_plots_violin_all_slices_ground_truth[[spes_metadata_index]][[metric]] <- plot_violin_all_slices_ground_truth(spes_table_subset, 
                                                                                                                          metric, 
                                                                                                                          metric_df_lists3D[[spes_metadata_index]][[metric]], 
                                                                                                                          metric_df_lists2D[[spes_metadata_index]][[metric]], 
                                                                                                                          arrangement_parameters[[arrangement]], 
                                                                                                                          plots_metadata[[shape]])
    }
  }
}


# Put plots into a pdf
setwd("~/R/plots/S1")
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics_set1 <- c("AMD", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC")
metrics_set2 <- c("MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")

pdf("metric_plots_violin_all_slices_with_ground_truth.pdf", width = 25, height = 10)

for (metric in metrics_set1) {
  for (shape in shapes) {
    curr_metric_plots <- list()
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[arrangement]] <- metric_plots_violin_all_slices_ground_truth[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
    }
    plot <- plot_grid(plotlist = curr_metric_plots,
                      nrow = 1, 
                      ncol = length(arrangements))
    print(plot)
  }
}

for (metric in metrics_set2) {
  curr_metric_plots <- list()
  for (shape in shapes) {
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[spes_metadata_index]] <- metric_plots_violin_all_slices_ground_truth[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
    }
  }
  plot <- plot_grid(plotlist = curr_metric_plots,
                    nrow = length(shapes), 
                    ncol = length(arrangements))
  print(plot)
}

dev.off()




