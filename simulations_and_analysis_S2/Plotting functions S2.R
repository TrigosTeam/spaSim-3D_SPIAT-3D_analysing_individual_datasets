library(cowplot)
library(ggplot2)
library(S4Vectors)
library(stringr)
library(dplyr)

### Utility function to get metric cell types -----
get_metric_cell_types <- function(metric) {
  # Get metric_cell_types
  if (metric %in% c("AMD", "ACIN", "CKR", "CLR", "COO", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC")) {
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
  else if (metric %in% c("prop_SAC", "prop_prevalence", "prop_AUC")) {
    metric_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("entropy_SAC", "entropy_prevalence", "entropy_AUC")) {
    metric_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))
  }
  else {
    stop("metric not found")
  }
  return(metric_cell_types)
}

get_metric_cell_types_singular <- function(metric) {
  # Get metric_cell_types
  if (metric %in% c("AMD", "ACIN", "CKR", "CLR", "COO", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "MS", "NMS", "MS_AUC", "NMS_AUC", "ACINP", "ACINP_AUC", "prop_SAC", "prop_prevalence", "prop_AUC")) {
    metric_cell_types <- data.frame(ref = c("A"), tar = c("B"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("AE", "AE_AUC")) {
    metric_cell_types <- data.frame(ref = c("A"), tar = c("A,B"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("entropy_SAC", "entropy_prevalence", "entropy_AUC")) {
    metric_cell_types <- data.frame(cell_types = c("A,B"))
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
  if (metric %in% c("AMD", "ACIN", "CKR", "CLR", "COO", "MS", "NMS", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "MS_AUC", "NMS_AUC", "prop_SAC", "prop_prevalence", "prop_AUC")) {
    metric_df <- metric_df[metric_df$reference == metric_cell_types[index, "ref"] & metric_df$target == metric_cell_types[index, "tar"], ] 
  }
  else if (metric %in% c("ACINP", "AE", "ACINP_AUC", "AE_AUC")) {
    metric_df <- metric_df[metric_df$reference == metric_cell_types[index, "ref"], ] 
  }
  else if (metric %in% c("entropy_SAC", "entropy_prevalence", "entropy_AUC")) {
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
  
  if (metric %in% c("AMD", "ACIN", "ACINP", "AE", "CKR", "CLR", "COO", "MS", "NMS", "ACIN_AUC", "ACINP_AUC", "AE_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "MS_AUC", "NMS_AUC", "prop_SAC", "prop_prevalence", "prop_AUC")) {
    title <- ggdraw() +
      draw_label(paste("Reference:", metric_cell_types[index, "ref"], "Target:", metric_cell_types[index, "tar"]),
                 fontface = 'bold')
  }
  else if (metric %in% c("entropy_SAC", "entropy_prevalence", "entropy_AUC")) {
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
    alpha <- as.numeric(x_sci[ , 1])
    power <- as.integer(x_sci[ , 2])
    paste(alpha, power, sep = "e")
  }
  
  create_plot <- function(data, x_aes, y_aes, title = "") {
    
    size <- 0.5
    
    plot <- ggplot(data, aes_string(x = x_aes, y = y_aes)) +
      labs(title = title, x = x_aes, y = y_aes) +
      theme_bw()
    
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
      title <- plot_def$title
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
    
    breaks <- pretty(c(min(data[[color_aes]]), max(data[[color_aes]])), n = 3)
    
    plot <- ggplot(data, aes_string(x = x_aes, y = y_aes, group = group_aes, color = color_aes)) +
      labs(title = title, x = x_aes, y = y_aes) +
      theme_bw() +
      geom_line()
    
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
    
    plot <- ggplot(data, aes_string(x = x_aes, y = y_aes, color = color_aes)) +
      labs(title = title, x = x_aes, y = y_aes) +
      theme_bw() +
      xlim(min(c(data[[x_aes]], data[[y_aes]])), max(c(data[[x_aes]], data[[y_aes]]))) +
      ylim(min(c(data[[x_aes]], data[[y_aes]])), max(c(data[[x_aes]], data[[y_aes]])))
    
    
    # Use scientific notation for ellipsoid volume
    if (color_aes == "E_volume") {
      plot <- plot + 
        geom_point(size = size) + 
        geom_abline(intercept = 0, slope = 1, color = "black", linetype = "longdash") +
        scale_color_continuous(labels = formatCustomSci)
    }
    else if (typeof(data[[color_aes]]) == "double") {
      breaks <- pretty(c(min(data[[color_aes]]), max(data[[color_aes]])), n = 3)
      plot <- plot + 
        geom_point(size = size) + 
        geom_abline(intercept = 0, slope = 1, color = "black", linetype = "longdash") +
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
      draw_label(label, fontface = 'bold')
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

### Function to compare 3D vs 2D (not annotating for arrangement or shape) all slices ---------------
plot_3D_vs_2D_metric_all_slices_no_annotating <- function(metric, 
                                                          metric_df3D,
                                                          metric_df2D,
                                                          plots_metadata) {
  
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
    
    plot_df <- data.frame(row.names = rownames(metric_df3D_subset))
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

### Function to compare 3D vs 2D (not annotating for arrangement or shape) random slice ---------------
plot_3D_vs_2D_metric_random_slice_no_annotating <- function(metric, 
                                                            metric_df3D,
                                                            metric_df2D,
                                                            plots_metadata) {
  
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
    
    wilcox_test  <- wilcox.test(data[[x_aes]], data[[y_aes]], paired = TRUE)
    p_value <- wilcox_test$p.value
    if (p_value == 0) p_value <- 2.2e-308
    if (0 < p_value && p_value < 1e-3)  {
      p_value <- formatCustomSci(p_value)
    }
    else {
      p_value  <- round(p_value, 3)
    }
    title <- paste("p =", p_value)
    
    plot <- ggplot(data, aes_string(x = x_aes, y = y_aes)) +
      labs(x = x_aes, y = y_aes) +
      theme_bw() +
      ggtitle(title) +
      theme(plot.title = element_text(size = 9)) +
      xlim(min(c(data[[x_aes]], data[[y_aes]]), na.rm = T), max(c(data[[x_aes]], data[[y_aes]]), na.rm = T)) +
      ylim(min(c(data[[x_aes]], data[[y_aes]]), na.rm = T), max(c(data[[x_aes]], data[[y_aes]]), na.rm = T)) +
      geom_point(size = 0.5) +
      geom_abline(intercept = 0, slope = 1, color = "red", linetype = "longdash")
    
    return(plot)
  }
  
  # Put plots into an organised list
  plots_list <- list()
  
  for (i in seq(nrow(metric_cell_types))) {
    
    # Subset metric_df for chosen pair/cell types
    metric_df3D_subset <- subset_metric_df(metric, metric_df3D, metric_cell_types, i)
    metric_df2D_subset <- subset_metric_df(metric, metric_df2D, metric_cell_types, i)
    
    plot_df <- data.frame(row.names = rownames(metric_df3D_subset))
    plot_df[[paste(metric, "3D", sep = "_")]] <- metric_df3D_subset[[metric]]

    # Choose a random slice from metric_df2D_subset
    n_slices <- length(unique(metric_df2D[["slice"]]))
    
    metric_df2D_subset$key <- paste(metric_df2D_subset[["spe"]], metric_df2D_subset[["slice"]], sep = "_")
    plot_df[[paste(metric, "2D", sep = "_")]] <- metric_df2D_subset[metric_df2D_subset$key %in% paste(unique(metric_df2D_subset[["spe"]]), sample(seq(n_slices), nrow(metric_df3D_subset), replace = TRUE), sep = "_"), 
                                                                    metric]
    
    # Factor
    if (!is.null(plot_df$shape)) plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    if (!is.null(plot_df$slice)) plot_df$slice <- as.character(plot_df$slice)
    
    # Generate plots based on plots_metadata, use final column of metric_cell_types
    plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]] <- lapply(plots_metadata, function(plot_def) {
      x_aes <- paste(metric, "3D", sep = "_")
      y_aes <- paste(metric, "2D", sep = "_")
      title <- plot_def$title
      plot <- create_plot(data = plot_df, x_aes = x_aes, y_aes = y_aes, title = title)
      return(plot)
    })
  }
  
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
  
  return(combined_plots)
}

### Function to compare ERROR vs 2D (not annotating for arrangement or shape) random slice ---------------
plot_error_vs_2D_metric_random_slice_no_annotating <- function(metric, 
                                                               metric_df3D,
                                                               metric_df2D,
                                                               plots_metadata) {
  
  # Get metric_cell_types
  metric_cell_types <- get_metric_cell_types(metric)
  
  # Define plotting function
  formatCustomSci <- function(x) {
    x_sci <- str_split_fixed(formatC(x, format = "e"), "e", 2)
    alpha <- as.numeric(x_sci[ , 1])
    power <- as.integer(x_sci[ , 2])
    paste(alpha, power, sep = "e")
  }
  
  create_plot <- function(data, x_aes, y_aes, title = "") {
    
    plot <- ggplot(data, aes_string(x = x_aes, y = y_aes)) +
      labs(title = title, x = x_aes, y = y_aes) +
      theme_bw() +
      ylim(min(c(0, data[[y_aes]]), na.rm = T), max(c(0, data[[y_aes]]), na.rm = T)) +
      geom_point(size = 0.5) +
      geom_abline(intercept = 0, slope = 0, color = "red", linetype = "longdash")
    
    return(plot)
  }
  
  # Put plots into an organised list
  plots_list <- list()
  
  for (i in seq(nrow(metric_cell_types))) {
    
    # Subset metric_df for chosen pair/cell types
    metric_df3D_subset <- subset_metric_df(metric, metric_df3D, metric_cell_types, i)
    metric_df2D_subset <- subset_metric_df(metric, metric_df2D, metric_cell_types, i)
    
    plot_df <- data.frame(row.names = rownames(metric_df3D_subset))
    plot_df[[paste(metric, "3D", sep = "_")]] <- metric_df3D_subset[[metric]]
    
    # Choose a random slice from metric_df2D_subset
    n_slices <- length(unique(metric_df2D[["slice"]]))
    
    metric_df2D_subset$key <- paste(metric_df2D_subset[["spe"]], metric_df2D_subset[["slice"]], sep = "_")
    plot_df[[paste(metric, "2D", sep = "_")]] <- metric_df2D_subset[metric_df2D_subset$key %in% paste(unique(metric_df2D_subset[["spe"]]), sample(seq(n_slices), nrow(metric_df3D_subset), replace = TRUE), sep = "_"), 
                                                                    metric]
    
    # Add error
    plot_df[[paste(metric, "error", sep = "_")]] <- 100 * (plot_df[[paste(metric, "2D", sep = "_")]] - plot_df[[paste(metric, "3D", sep = "_")]]) / (plot_df[[paste(metric, "3D", sep = "_")]])
    
    # Factor
    if (!is.null(plot_df$shape)) plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    if (!is.null(plot_df$slice)) plot_df$slice <- as.character(plot_df$slice)
    
    # Generate plots based on plots_metadata, use final column of metric_cell_types
    plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]] <- lapply(plots_metadata, function(plot_def) {
      x_aes <- paste(metric, "3D", sep = "_")
      y_aes <- paste(metric, "error", sep = "_")
      title <- plot_def$title
      plot <- create_plot(data = plot_df, x_aes = x_aes, y_aes = y_aes, title = title)
      return(plot)
    })
  }
  
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
  
  return(combined_plots)
}

### Function to compare 3D vs 2D (not annotating for arrangement or shape) random slice (one cell type pair A/B or equivalent) ---------------
plot_3D_vs_2D_metric_random_slice_no_annotating_one_cell_type_pair <- function(metric, 
                                                                               metric_df3D,
                                                                               metric_df2D,
                                                                               plots_metadata) {
  
  # Get metric_cell_types
  metric_cell_types <- get_metric_cell_types_singular(metric)
  
  # Define plotting function
  formatCustomSci <- function(x) {
    x_sci <- str_split_fixed(formatC(x, format = "e"), "e", 2)
    alpha <- round(as.numeric(x_sci[ , 1]), 1)
    power <- as.integer(x_sci[ , 2])
    paste(alpha, power, sep = "e")
  }
  
  create_plot <- function(data, x_aes, y_aes, title = "") {
    
    wilcox_test  <- wilcox.test(data[[x_aes]], data[[y_aes]], paired = TRUE)
    p_value <- wilcox_test$p.value
    if (p_value == 0) p_value <- 2.2e-308
    if (0 < p_value && p_value < 1e-3)  {
      p_value <- formatCustomSci(p_value)
    }
    else {
      p_value  <- round(p_value, 3)
    }
    title <- paste("p =", p_value)
    
    plot <- ggplot(data, aes_string(x = x_aes, y = y_aes)) +
      labs(x = x_aes, y = y_aes) +
      theme_bw() +
      ggtitle(title) +
      theme(plot.title = element_text(size = 9)) +
      xlim(min(c(data[[x_aes]], data[[y_aes]]), na.rm = T), max(c(data[[x_aes]], data[[y_aes]]), na.rm = T)) +
      ylim(min(c(data[[x_aes]], data[[y_aes]]), na.rm = T), max(c(data[[x_aes]], data[[y_aes]]), na.rm = T)) +
      geom_point(size = 0.5) +
      geom_abline(intercept = 0, slope = 1, color = "red", linetype = "longdash")
    
    return(plot)
  }
  
  # Put plots into an organised list
  plots_list <- list()
  
  for (i in seq(nrow(metric_cell_types))) {
    
    # Subset metric_df for chosen pair/cell types
    metric_df3D_subset <- subset_metric_df(metric, metric_df3D, metric_cell_types, i)
    metric_df2D_subset <- subset_metric_df(metric, metric_df2D, metric_cell_types, i)
    
    plot_df <- data.frame(row.names = rownames(metric_df3D_subset))
    plot_df[[paste(metric, "3D", sep = "_")]] <- metric_df3D_subset[[metric]]
    
    # Choose a random slice from metric_df2D_subset
    n_slices <- length(unique(metric_df2D[["slice"]]))
    
    metric_df2D_subset$key <- paste(metric_df2D_subset[["spe"]], metric_df2D_subset[["slice"]], sep = "_")
    plot_df[[paste(metric, "2D", sep = "_")]] <- metric_df2D_subset[metric_df2D_subset$key %in% paste(unique(metric_df2D_subset[["spe"]]), sample(seq(n_slices), nrow(metric_df3D_subset), replace = TRUE), sep = "_"), 
                                                                    metric]
    
    # Factor
    if (!is.null(plot_df$shape)) plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    if (!is.null(plot_df$slice)) plot_df$slice <- as.character(plot_df$slice)
    
    # Generate plots based on plots_metadata, use final column of metric_cell_types
    plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]] <- lapply(plots_metadata, function(plot_def) {
      x_aes <- paste(metric, "3D", sep = "_")
      y_aes <- paste(metric, "2D", sep = "_")
      title <- plot_def$title
      plot <- create_plot(data = plot_df, x_aes = x_aes, y_aes = y_aes, title = title)
      return(plot)
    })
  }
  
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
  
  return(combined_plots)
}

### Function to compare ERROR vs 2D (not annotating for arrangement or shape) random slice (one cell type pair A/B or equivalent) ---------------
plot_error_vs_2D_metric_random_slice_no_annotating_one_cell_type_pair <- function(metric, 
                                                                                  metric_df3D,
                                                                                  metric_df2D,
                                                                                  plots_metadata) {
  
  # Get metric_cell_types
  metric_cell_types <- get_metric_cell_types_singular(metric)
  
  # Define plotting function
  formatCustomSci <- function(x) {
    x_sci <- str_split_fixed(formatC(x, format = "e"), "e", 2)
    alpha <- as.numeric(x_sci[ , 1])
    power <- as.integer(x_sci[ , 2])
    paste(alpha, power, sep = "e")
  }
  
  create_plot <- function(data, x_aes, y_aes, title = "") {
    
    plot <- ggplot(data, aes_string(x = x_aes, y = y_aes)) +
      labs(title = title, x = x_aes, y = y_aes) +
      theme_bw() +
      ylim(min(c(0, data[[y_aes]]), na.rm = T), max(c(0, data[[y_aes]]), na.rm = T)) +
      geom_point(size = 0.5) +
      geom_abline(intercept = 0, slope = 0, color = "red", linetype = "longdash")
    
    return(plot)
  }
  
  # Put plots into an organised list
  plots_list <- list()
  
  for (i in seq(nrow(metric_cell_types))) {
    
    # Subset metric_df for chosen pair/cell types
    metric_df3D_subset <- subset_metric_df(metric, metric_df3D, metric_cell_types, i)
    metric_df2D_subset <- subset_metric_df(metric, metric_df2D, metric_cell_types, i)
    
    plot_df <- data.frame(row.names = rownames(metric_df3D_subset))
    plot_df[[paste(metric, "3D", sep = "_")]] <- metric_df3D_subset[[metric]]
    
    # Choose a random slice from metric_df2D_subset
    n_slices <- length(unique(metric_df2D[["slice"]]))
    
    metric_df2D_subset$key <- paste(metric_df2D_subset[["spe"]], metric_df2D_subset[["slice"]], sep = "_")
    plot_df[[paste(metric, "2D", sep = "_")]] <- metric_df2D_subset[metric_df2D_subset$key %in% paste(unique(metric_df2D_subset[["spe"]]), sample(seq(n_slices), nrow(metric_df3D_subset), replace = TRUE), sep = "_"), 
                                                                    metric]
    
    # Add error
    plot_df[[paste(metric, "error", sep = "_")]] <- 100 * (plot_df[[paste(metric, "2D", sep = "_")]] - plot_df[[paste(metric, "3D", sep = "_")]]) / (plot_df[[paste(metric, "3D", sep = "_")]])
    
    # Factor
    if (!is.null(plot_df$shape)) plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    if (!is.null(plot_df$slice)) plot_df$slice <- as.character(plot_df$slice)
    
    # Generate plots based on plots_metadata, use final column of metric_cell_types
    plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]] <- lapply(plots_metadata, function(plot_def) {
      x_aes <- paste(metric, "3D", sep = "_")
      y_aes <- paste(metric, "error", sep = "_")
      title <- plot_def$title
      plot <- create_plot(data = plot_df, x_aes = x_aes, y_aes = y_aes, title = title)
      return(plot)
    })
  }
  
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
  
  return(combined_plots)
}




### Function to ERROR box plot for each metric (not annotating for arrangement or shape) random slice (one cell type pair A/B or equivalent) ---------------
plot_error_box_plot_metric_random_slice_no_annotating_one_cell_type_pair <- function(metrics,
                                                                                     metric_df3D,
                                                                                     metric_df2D) {
  
  # Put all data into this data frame
  combined_plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_cell_types
    metric_cell_types <- get_metric_cell_types_singular(metric)
    
    for (i in seq(nrow(metric_cell_types))) {
      
      metric_df3D_subset <- metric_df3D[[metric]]
      metric_df2D_subset <- metric_df2D[[metric]]
      
      # Subset metric_df for chosen pair/cell types
      metric_df3D_subset <- subset_metric_df(metric, metric_df3D_subset, metric_cell_types, i)
      metric_df2D_subset <- subset_metric_df(metric, metric_df2D_subset, metric_cell_types, i)
      
      plot_df <- data.frame(row.names = rownames(metric_df3D_subset))
      plot_df[[paste(metric, "3D", sep = "_")]] <- metric_df3D_subset[[metric]]
      
      # Choose a random slice from metric_df2D_subset
      n_slices <- length(unique(metric_df2D_subset[["slice"]]))
      
      metric_df2D_subset$key <- paste(metric_df2D_subset[["spe"]], metric_df2D_subset[["slice"]], sep = "_")
      plot_df[[paste(metric, "2D", sep = "_")]] <- metric_df2D_subset[metric_df2D_subset$key %in% paste(unique(metric_df2D_subset[["spe"]]), sample(seq(n_slices), nrow(metric_df3D_subset), replace = TRUE), sep = "_"), 
                                                                      metric]
      
      # Add error
      plot_df[[paste(metric, "error", sep = "_")]] <- 100 * (plot_df[[paste(metric, "2D", sep = "_")]] - plot_df[[paste(metric, "3D", sep = "_")]]) / (plot_df[[paste(metric, "3D", sep = "_")]])
      
      # Factor
      if (!is.null(plot_df$shape)) plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
      if (!is.null(plot_df$slice)) plot_df$slice <- as.character(plot_df$slice)
      
      combined_plot_df <- rbind(combined_plot_df, data.frame(error = plot_df[[paste(metric, "error", sep = "_")]],
                                                             metric = metric))
    }    
  }
  
  # Remove outliers using IQR method for each metric group
  df_clean <- combined_plot_df %>%
    group_by(metric) %>%
    filter(
      error > quantile(error, 0.25, na.rm = TRUE) - 20 * IQR(error, na.rm = TRUE) &
        error < quantile(error, 0.75, na.rm = TRUE) + 20 * IQR(error, na.rm = TRUE)
    )

  # Plot cleaned data
  fig <- ggplot(df_clean, aes(x = metric, y = error)) +
    geom_boxplot() +
    labs(title = "Error Distribution by Metric",
         x = "Metric",
         y = "Error (%)") +
    theme_minimal()
  
  return(fig)
}



