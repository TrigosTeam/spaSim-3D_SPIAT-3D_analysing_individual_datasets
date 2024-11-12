library(cowplot)
library(ggplot2)
library(S4Vectors)
library(stringr)
library(dplyr)

### Utility functions -------------------------
get_gradient <- function(metric) {
  if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR")) return("radius")
  return("threshold")
}

### Read metric_df_lists --------------------------------------------------------------
setwd("~/R/spaSim-3D/scripts/simulations and analysis S1/S1 data")
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

gradient_radii_metrics <- c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR")


for (arrangement in arrangements) {
  
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    for (metric in gradient_radii_metrics) {
      metric_AUC_name <- paste(metric, "AUC", sep = "_")
      
      if (metric %in% c("MS", "NMS", "ACIN", "CKR")) {
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
    
    # prop_AUC 3D
    prop_prevalence_df <- metric_df_lists3D[[spes_metadata_index]][["prop_prevalence"]]
    prop_prevalence_df$prop_AUC <- apply(prop_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
    prop_AUC_df <- prop_prevalence_df[ , c("spe", "reference", "target", "prop_AUC")]
    metric_df_lists3D[[spes_metadata_index]][["prop_AUC"]] <- prop_AUC_df
    
    # entropy_AUC 3D
    entropy_prevalence_df <- metric_df_lists3D[[spes_metadata_index]][["entropy_prevalence"]]
    entropy_prevalence_df$entropy_AUC <- apply(entropy_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
    entropy_AUC_df <- entropy_prevalence_df[ , c("spe", "cell_types", "entropy_AUC")]
    metric_df_lists3D[[spes_metadata_index]][["entropy_AUC"]] <- entropy_AUC_df
    
    # prop_AUC 2D
    prop_prevalence_df <- metric_df_lists2D[[spes_metadata_index]][["prop_prevalence"]]
    prop_prevalence_df$prop_AUC <- apply(prop_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
    prop_AUC_df <- prop_prevalence_df[ , c("spe", "slice", "reference", "target", "prop_AUC")]
    metric_df_lists2D[[spes_metadata_index]][["prop_AUC"]] <- prop_AUC_df
    
    # entropy_AUC 2D
    entropy_prevalence_df <- metric_df_lists2D[[spes_metadata_index]][["entropy_prevalence"]]
    entropy_prevalence_df$entropy_AUC <- apply(entropy_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
    entropy_AUC_df <- entropy_prevalence_df[ , c("spe", "slice", "cell_types", "entropy_AUC")]
    metric_df_lists2D[[spes_metadata_index]][["entropy_AUC"]] <- entropy_AUC_df
  }
}


### Get plots for 3D metric analysis (non-gradient) ------------------------------------------
# Read spes_table
setwd("~/R/spaSim-3D/scripts/simulations and analysis S1/S1 data")
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
metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")

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

metrics_set1 <- c("AMD",  "ACIN_AUC", "CKR_AUC")
metrics_set2 <- c("MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")

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
setwd("~/R/spaSim-3D/scripts/simulations and analysis S1/S1 data")
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
metrics <- c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_prevalence", "entropy_prevalence")

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

metrics_set1 <- c("ACIN", "CKR")
metrics_set2 <- c("MS", "NMS", "ACINP", "AE", "prop_prevalence", "entropy_prevalence")

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
setwd("~/R/spaSim-3D/scripts/simulations and analysis S1/S1 data")
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
metrics <- c("AMD", "MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "prop_prevalence", "prop_AUC", "entropy_SAC", "entropy_prevalence", "entropy_AUC")

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
      if (metric %in% c("AMD", "prop_SAC", "entropy_SAC", "prop_AUC", "entropy_AUC")) {
        metric_plots2D[[spes_metadata_index]][[metric]] <- plot_non_gradient_metric(spes_table_subset, 
                                                                                    metric, 
                                                                                    metric_df_lists2D_subset[[spes_metadata_index]][[metric]], 
                                                                                    arrangement_parameters[[arrangement]], 
                                                                                    non_gradient_plots_metadata[[shape]])
      }
      else if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_prevalence", "entropy_prevalence")) {
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
# metrics <- c("AMD", "MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "prop_prevalence", "prop_AUC", "entropy_SAC", "entropy_prevalence", "entropy_AUC")

metrics_set1 <- c("AMD", "ACIN", "CKR")
metrics_set2 <- c("MS", "NMS", "ACINP", "AE", "prop_SAC", "prop_prevalence", "prop_AUC", "entropy_SAC", "entropy_prevalence", "entropy_AUC")

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
setwd("~/R/spaSim-3D/scripts/simulations and analysis S1/S1 data")
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
metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")

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
metrics_set1 <- c("AMD",  "ACIN_AUC", "CKR_AUC")
metrics_set2 <- c("MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")

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
setwd("~/R/spaSim-3D/scripts/simulations and analysis S1/S1 data")
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
metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")

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
metrics_set1 <- c("AMD",  "ACIN_AUC", "CKR_AUC")
metrics_set2 <- c("MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")

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
setwd("~/R/spaSim-3D/scripts/simulations and analysis S1/S1 data")
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
metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")

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
metrics_set1 <- c("AMD",  "ACIN_AUC", "CKR_AUC")
metrics_set2 <- c("MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")

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
setwd("~/R/spaSim-3D/scripts/simulations and analysis S1/S1 data")
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
metrics <- c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_prevalence",  "entropy_prevalence")

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
metrics_set1 <- c("ACIN", "CKR")
metrics_set2 <- c("MS", "NMS", "ACINP", "AE", "prop_prevalence", "entropy_prevalence")

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
setwd("~/R/spaSim-3D/scripts/simulations and analysis S1/S1 data")
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
metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")

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
metrics_set1 <- c("AMD", "ACIN_AUC", "CKR_AUC")
metrics_set2 <- c("MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")

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
setwd("~/R/spaSim-3D/scripts/simulations and analysis S1/S1 data")
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
metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")

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
metrics_set1 <- c("AMD", "ACIN_AUC", "CKR_AUC")
metrics_set2 <- c("MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")

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
setwd("~/R/spaSim-3D/scripts/simulations and analysis S1/S1 data")
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
metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")

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
metrics_set1 <- c("AMD", "ACIN_AUC", "CKR_AUC")
metrics_set2 <- c("MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")

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



