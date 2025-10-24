library(cowplot)
library(ggplot2)

### Utility functions -------------------------
get_gradient <- function(metric) {
  if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "CLR", "COO", "CGR")) return("radius")
  return("threshold")
}

### Read metric_df_lists --------------------------------------------------------------
setwd("~/R/SPIAT-3D_benchmarking/simulations_and_analysis_S2/S2_data")
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

gradient_radii_metrics <- c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "CLR", "COO", "CGR")


for (arrangement in arrangements) {
  
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    for (metric in gradient_radii_metrics) {
      metric_AUC_name <- paste(metric, "AUC", sep = "_")
      
      if (metric %in% c("MS", "NMS", "ACIN", "CKR", "CLR", "COO", "CGR")) {
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


### Get plots for 3D metric analysis ------------------------------------------
# Read spes_table
setwd("~/R/SPIAT-3D_benchmarking/simulations_and_analysis_S2/S2_data")
spes_table <- read.table("spes_table.csv")
spes_table$cluster_prop_B <- 1 - spes_table$cluster_prop_A
spes_table$distance <- 450 - spes_table$cluster_x_coord
spes_table$E_volume <- (4/3) * pi * spes_table$E_radius_x * spes_table$E_radius_y * spes_table$E_radius_z

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


metric_plots3D <- list(mixed_ellipsoid = list(),
                       mixed_network = list(),
                       ringed_ellipsoid = list(),
                       ringed_network = list(),
                       separated_ellipsoid = list(),
                       separated_network = list())

for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    for (metric in metrics) {
      if (metric %in% c("AMD", "PBSAC", "EBSAC", "PBP_AUC", "EBP_AUC")) {
        metric_plots3D[[spes_metadata_index]][[metric]] <- plot_non_gradient_metric(spes_table, 
                                                                                    metric, 
                                                                                    metric_df_lists3D[[spes_metadata_index]][[metric]], 
                                                                                    arrangement_parameters[[arrangement]], 
                                                                                    non_gradient_plots_metadata[[shape]])
      }
      else if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "CLR", "COO", "CGR", "PBP", "EBP")) {
        metric_plots3D[[spes_metadata_index]][[metric]] <- plot_gradient_metric(spes_table, 
                                                                                metric,
                                                                                metric_df_lists3D[[spes_metadata_index]][[metric]], 
                                                                                arrangement_parameters[[arrangement]], 
                                                                                get_gradient(metric),
                                                                                gradient_plots_metadata[[shape]])
      }
    }
  }
}

# Put plots into a pdf
setwd("~/R/plots/S2")
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")

metrics_set1 <- c("AMD", "ACIN", "CKR", "CLR", "COO", "CGR")
metrics_set2 <- c("MS", "NMS", "ACINP", "AE", "PBSAC", "PBP", "PBP_AUC", "EBSAC", "EBP", "EBP_AUC")

pdf("plots3D.pdf", width = 25, height = 12)

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
setwd("~/R/SPIAT-3D_benchmarking/simulations_and_analysis_S2/S2_data")
spes_table <- read.table("spes_table.csv")
spes_table$cluster_prop_B <- 1 - spes_table$cluster_prop_A
spes_table$distance <- 450 - spes_table$cluster_x_coord
spes_table$E_volume <- (4/3) * pi * spes_table$E_radius_x * spes_table$E_radius_y * spes_table$E_radius_z

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
    
    for (metric in metrics) {
      if (metric %in% c("AMD", "PBSAC", "EBSAC", "PBP_AUC", "EBP_AUC")) {
        metric_plots2D[[spes_metadata_index]][[metric]] <- plot_non_gradient_metric(spes_table, 
                                                                                    metric, 
                                                                                    metric_df_lists2D_subset[[spes_metadata_index]][[metric]], 
                                                                                    arrangement_parameters[[arrangement]], 
                                                                                    non_gradient_plots_metadata[[shape]])
      }
      else if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "CLR", "COO", "CGR", "PBP", "EBP")) {
        metric_plots2D[[spes_metadata_index]][[metric]] <- plot_gradient_metric(spes_table, 
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
setwd("~/R/plots/S2")
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
setwd("~/R/SPIAT-3D_benchmarking/simulations_and_analysis_S2/S2_data")
spes_table <- read.table("spes_table.csv")
spes_table$cluster_prop_B <- 1 - spes_table$cluster_prop_A
spes_table$distance <- 450 - spes_table$cluster_x_coord
spes_table$E_volume <- (4/3) * pi * spes_table$E_radius_x * spes_table$E_radius_y * spes_table$E_radius_z


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

metric_plots2D_vs_3D <- list(mixed_ellipsoid = list(),
                             mixed_network = list(),
                             ringed_ellipsoid = list(),
                             ringed_network = list(),
                             separated_ellipsoid = list(),
                             separated_network = list())

for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    for (metric in metrics) {
      metric_plots2D_vs_3D[[spes_metadata_index]][[metric]] <- plot_3D_vs_2D_metric_one_slice(spes_table, 
                                                                                              metric, 
                                                                                              metric_df_lists3D[[spes_metadata_index]][[metric]],
                                                                                              metric_df_lists2D_subset[[spes_metadata_index]][[metric]], 
                                                                                              arrangement_parameters[[arrangement]], 
                                                                                              plots_metadata[[shape]])
    }
  }
}


# Put plots into a pdf
setwd("~/R/plots/S2")
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics_set1 <- c("AMD",  "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC")
metrics_set2 <- c("MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")

pdf("plots2D_vs_3D_middle_slice.pdf", width = 25, height = 10)

for (metric in metrics_set1) {
  for (shape in shapes) {
    curr_metric_plots <- list()
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[arrangement]] <- metric_plots2D_vs_3D[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
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
      curr_metric_plots[[spes_metadata_index]] <- metric_plots2D_vs_3D[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
    }
  }
  plot <- plot_grid(plotlist = curr_metric_plots,
                    nrow = length(shapes), 
                    ncol = length(arrangements))
  print(plot)
}

dev.off()






### Get plots with 2D (all slices) on the x-axis and 3D on the y-axis ----------------
setwd("~/R/SPIAT-3D_benchmarking/simulations_and_analysis_S2/S2_data")
spes_table <- read.table("spes_table.csv")
spes_table$cluster_prop_B <- 1 - spes_table$cluster_prop_A
spes_table$distance <- 450 - spes_table$cluster_x_coord
spes_table$E_volume <- (4/3) * pi * spes_table$E_radius_x * spes_table$E_radius_y * spes_table$E_radius_z

# Subset metric_df_lists2D to include 1, 2, 3 slice (I know it says all above, deal with it)
metric_df_lists2D_subset <- metric_df_lists2D
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    curr_list <- metric_df_lists2D_subset[[spes_metadata_index]]
    
    for (i in seq(length(curr_list))) {
      curr_df <- curr_list[[i]]
      curr_df <- curr_df[curr_df$slice %in% c(1, 2, 3), ]
      curr_list[[i]] <- curr_df
    }
    metric_df_lists2D_subset[[spes_metadata_index]] <- curr_list
  }
}

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
    
    for (metric in metrics) {
      metric_plots2D_vs_3D_all_slices[[spes_metadata_index]][[metric]] <- plot_3D_vs_2D_metric_all_slices(spes_table, 
                                                                                                          metric, 
                                                                                                          metric_df_lists3D[[spes_metadata_index]][[metric]],
                                                                                                          metric_df_lists2D_subset[[spes_metadata_index]][[metric]], 
                                                                                                          arrangement_parameters[[arrangement]], 
                                                                                                          plots_metadata[[shape]])
    }
  }
}



# Put plots into a pdf
setwd("~/R/plots/S2")
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics_set1 <- c("AMD",  "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC")
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






### Get plots with 2D (all slices) on the x-axis and 3D on the y-axis (not annotating for arrangement or shape) ----------------

# Set up plots metadata
plots_metadata <- list(
  temp <- list(x_aes = "2D", y_aes = "3D", color_aes = "slice")
)

# Generate plots and plots into a list
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")


# Merge lists in metric_lists
metric_df_lists3D_merged <- list()
metric_df_lists2D_merged <- list()

i <- 1
for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    for (metric in metrics) {
      if (i == 1)  {
        metric_df_lists3D_merged[[metric]] <- data.frame()
        metric_df_lists2D_merged[[metric]] <- data.frame()
      }
      if (i > 1) {
        temp <- nrow(metric_df_lists3D[[spes_metadata_index]][[metric]])
        n_slices <- length(unique(metric_df_lists2D[[spes_metadata_index]][[metric]][["slice"]]))
        if (metric %in% c("AMD", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC")) {
          metric_df_lists3D[[spes_metadata_index]][[metric]][["spe"]] <- 
            paste("spe", rep(seq((temp/4) * (i - 1) + 1, (temp/4) * (i - 1) + (temp/4)), each = 4), sep = "_")
          metric_df_lists2D[[spes_metadata_index]][[metric]][["spe"]] <-
            paste("spe", rep(seq((temp/4) * (i - 1) + 1, (temp)/4 * (i - 1) + (temp/4)), each = 4 * n_slices), sep = "_")
        }
        else {
          metric_df_lists3D[[spes_metadata_index]][[metric]][["spe"]] <- 
            paste("spe", rep(seq((temp/2) * (i - 1) + 1, (temp/2) * (i - 1) + (temp/2)), each = 2), sep = "_")
          metric_df_lists2D[[spes_metadata_index]][[metric]][["spe"]] <-
            paste("spe", rep(seq((temp/2) * (i - 1) + 1, (temp/2) * (i - 1) + (temp/2)), each = 2 * n_slices), sep = "_")
        }
      }
      metric_df_lists3D_merged[[metric]] <- rbind(metric_df_lists3D_merged[[metric]], metric_df_lists3D[[spes_metadata_index]][[metric]])
      metric_df_lists2D_merged[[metric]] <- rbind(metric_df_lists2D_merged[[metric]], metric_df_lists2D[[spes_metadata_index]][[metric]])
    }
    
    i <- i + 1
  }
}

metric_plots_3D_vs_2D_all_slices_no_annotating <- list()

for (metric in metrics) {
  metric_plots_3D_vs_2D_all_slices_no_annotating[[metric]] <- plot_3D_vs_2D_metric_all_slices_no_annotating(metric, 
                                                                                                            metric_df_lists3D_merged[[metric]],
                                                                                                            metric_df_lists2D_merged[[metric]], 
                                                                                                            plots_metadata)
}



# Put plots into a pdf
setwd("~/R/plots/S2")
metrics_set1 <- c("AMD",  "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC")
metrics_set2 <- c("MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC")
metrics_set3 <- c("PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")

pdf("plot_3D_vs_2D_metric_all_slices_no_annotating.pdf", width = 12, height = 12)

curr_metric_plots1 <- list()
for (metric in metrics_set1) {
    curr_metric_plots1[[metric]] <- metric_plots_3D_vs_2D_all_slices_no_annotating[[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
}
plot1 <- plot_grid(plotlist = curr_metric_plots1,
                   nrow = 1, 
                   ncol = length(metrics_set1))
print(plot1)

curr_metric_plots2 <- list()
for (metric in metrics_set2) {
  curr_metric_plots2[[metric]] <- metric_plots_3D_vs_2D_all_slices_no_annotating[[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
}
plot2 <- plot_grid(plotlist = curr_metric_plots2,
                   nrow = 2, 
                   ncol = 2)

curr_metric_plots3 <- list()
for (metric in metrics_set3) {
  curr_metric_plots3[[metric]] <- metric_plots_3D_vs_2D_all_slices_no_annotating[[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
}
plot3 <- plot_grid(plotlist = curr_metric_plots3,
                   nrow = 2, 
                   ncol = 2)

plot4 <- plot_grid(plot2, plot3,
                   nrow = 1,
                   ncol = 2)
print(plot4)
dev.off()






### Get plots with 2D (all slices) on the x-axis and 3D on the y-axis (not annotating for arrangement or shape and choosing random slice) ----------------

# Set up plots metadata
plots_metadata <- list(
  temp <- list(x_aes = "3D", y_aes = "2D")
)

# Generate plots and plots into a list
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")


# Merge lists in metric_lists
metric_df_lists3D_merged <- list()
metric_df_lists2D_merged <- list()

i <- 1
for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    for (metric in metrics) {
      if (i == 1)  {
        metric_df_lists3D_merged[[metric]] <- data.frame()
        metric_df_lists2D_merged[[metric]] <- data.frame()
      }
      if (i > 1) {
        temp <- nrow(metric_df_lists3D[[spes_metadata_index]][[metric]])
        n_slices <- length(unique(metric_df_lists2D[[spes_metadata_index]][[metric]][["slice"]]))
        if (metric %in% c("AMD", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC")) {
          metric_df_lists3D[[spes_metadata_index]][[metric]][["spe"]] <- 
            paste("spe", rep(seq((temp/4) * (i - 1) + 1, (temp/4) * (i - 1) + (temp/4)), each = 4), sep = "_")
          metric_df_lists2D[[spes_metadata_index]][[metric]][["spe"]] <-
            paste("spe", rep(seq((temp/4) * (i - 1) + 1, (temp)/4 * (i - 1) + (temp/4)), each = 4 * n_slices), sep = "_")
        }
        else {
          metric_df_lists3D[[spes_metadata_index]][[metric]][["spe"]] <- 
            paste("spe", rep(seq((temp/2) * (i - 1) + 1, (temp/2) * (i - 1) + (temp/2)), each = 2), sep = "_")
          metric_df_lists2D[[spes_metadata_index]][[metric]][["spe"]] <-
            paste("spe", rep(seq((temp/2) * (i - 1) + 1, (temp/2) * (i - 1) + (temp/2)), each = 2 * n_slices), sep = "_")
        }
      }
      metric_df_lists3D_merged[[metric]] <- rbind(metric_df_lists3D_merged[[metric]], metric_df_lists3D[[spes_metadata_index]][[metric]])
      metric_df_lists2D_merged[[metric]] <- rbind(metric_df_lists2D_merged[[metric]], metric_df_lists2D[[spes_metadata_index]][[metric]])
    }
    
    i <- i + 1
  }
}

metric_plots_3D_vs_2D_random_slice <- list()

for (metric in metrics) {
  metric_plots_3D_vs_2D_random_slice[[metric]] <- plot_3D_vs_2D_metric_random_slice_no_annotating(metric, 
                                                                                                  metric_df_lists3D_merged[[metric]],
                                                                                                  metric_df_lists2D_merged[[metric]], 
                                                                                                  plots_metadata)
}



# Put plots into a pdf
setwd("~/R/plots/S2")
metrics_set1 <- c("AMD",  "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC")
metrics_set2 <- c("MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC")
metrics_set3 <- c("PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")

pdf("plots2D_vs_3D_random_slice_all_cell_pairs.pdf", width = 12, height = 12)

curr_metric_plots1 <- list()
for (metric in metrics_set1) {
  curr_metric_plots1[[metric]] <- metric_plots_3D_vs_2D_random_slice[[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
}
plot1 <- plot_grid(plotlist = curr_metric_plots1,
                   nrow = 1, 
                   ncol = length(metrics_set1))
print(plot1)

curr_metric_plots2 <- list()
for (metric in metrics_set2) {
  curr_metric_plots2[[metric]] <- metric_plots_3D_vs_2D_random_slice[[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
}
plot2 <- plot_grid(plotlist = curr_metric_plots2,
                   nrow = 2, 
                   ncol = 2)

curr_metric_plots3 <- list()
for (metric in metrics_set3) {
  curr_metric_plots3[[metric]] <- metric_plots_3D_vs_2D_random_slice[[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
}
plot3 <- plot_grid(plotlist = curr_metric_plots3,
                   nrow = 2, 
                   ncol = 2)

plot4 <- plot_grid(plot2, plot3,
                   nrow = 1,
                   ncol = 2)
print(plot4)
dev.off()






### Get plots with 2D (all slices) on the x-axis and ERROR on the y-axis (not annotating for arrangement or shape and choosing random slice) ----------------

# Set up plots metadata
plots_metadata <- list(
  temp <- list(x_aes = "3D", y_aes = "error")
)

# Generate plots and plots into a list
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")


# Merge lists in metric_lists
metric_df_lists3D_merged <- list()
metric_df_lists2D_merged <- list()

i <- 1
for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    for (metric in metrics) {
      if (i == 1)  {
        metric_df_lists3D_merged[[metric]] <- data.frame()
        metric_df_lists2D_merged[[metric]] <- data.frame()
      }
      if (i > 1) {
        temp <- nrow(metric_df_lists3D[[spes_metadata_index]][[metric]])
        n_slices <- length(unique(metric_df_lists2D[[spes_metadata_index]][[metric]][["slice"]]))
        if (metric %in% c("AMD", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC")) {
          metric_df_lists3D[[spes_metadata_index]][[metric]][["spe"]] <- 
            paste("spe", rep(seq((temp/4) * (i - 1) + 1, (temp/4) * (i - 1) + (temp/4)), each = 4), sep = "_")
          metric_df_lists2D[[spes_metadata_index]][[metric]][["spe"]] <-
            paste("spe", rep(seq((temp/4) * (i - 1) + 1, (temp)/4 * (i - 1) + (temp/4)), each = 4 * n_slices), sep = "_")
        }
        else {
          metric_df_lists3D[[spes_metadata_index]][[metric]][["spe"]] <- 
            paste("spe", rep(seq((temp/2) * (i - 1) + 1, (temp/2) * (i - 1) + (temp/2)), each = 2), sep = "_")
          metric_df_lists2D[[spes_metadata_index]][[metric]][["spe"]] <-
            paste("spe", rep(seq((temp/2) * (i - 1) + 1, (temp/2) * (i - 1) + (temp/2)), each = 2 * n_slices), sep = "_")
        }
      }
      metric_df_lists3D_merged[[metric]] <- rbind(metric_df_lists3D_merged[[metric]], metric_df_lists3D[[spes_metadata_index]][[metric]])
      metric_df_lists2D_merged[[metric]] <- rbind(metric_df_lists2D_merged[[metric]], metric_df_lists2D[[spes_metadata_index]][[metric]])
    }
    
    i <- i + 1
  }
}

metric_plots_error_vs_2D_random_slice <- list()

for (metric in metrics) {
  metric_plots_error_vs_2D_random_slice[[metric]] <- plot_error_vs_2D_metric_random_slice_no_annotating(metric, 
                                                                                                        metric_df_lists3D_merged[[metric]],
                                                                                                        metric_df_lists2D_merged[[metric]], 
                                                                                                        plots_metadata)
}



# Put plots into a pdf
setwd("~/R/plots/S2")
metrics_set1 <- c("AMD",  "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC")
metrics_set2 <- c("MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC")
metrics_set3 <- c("PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")

pdf("plots_error_vs_3D_random_slice_all_cell_pairs.pdf", width = 12, height = 12)

curr_metric_plots1 <- list()
for (metric in metrics_set1) {
  curr_metric_plots1[[metric]] <- metric_plots_error_vs_2D_random_slice[[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
}
plot1 <- plot_grid(plotlist = curr_metric_plots1,
                   nrow = 1, 
                   ncol = length(metrics_set1))
print(plot1)

curr_metric_plots2 <- list()
for (metric in metrics_set2) {
  curr_metric_plots2[[metric]] <- metric_plots_error_vs_2D_random_slice[[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
}
plot2 <- plot_grid(plotlist = curr_metric_plots2,
                   nrow = 2, 
                   ncol = 2)

curr_metric_plots3 <- list()
for (metric in metrics_set3) {
  curr_metric_plots3[[metric]] <- metric_plots_error_vs_2D_random_slice[[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
}
plot3 <- plot_grid(plotlist = curr_metric_plots3,
                   nrow = 2, 
                   ncol = 2)

plot4 <- plot_grid(plot2, plot3,
                   nrow = 1,
                   ncol = 2)
print(plot4)
dev.off()






### Get plots with 2D (all slices) on the x-axis and 3D on the y-axis (not annotating for arrangement or shape and choosing random slice) (one cell type pair A/B or equivalent) ----------------

# Set up plots metadata
plots_metadata <- list(
  temp <- list(x_aes = "3D", y_aes = "2D")
)

# Generate plots and plots into a list
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")


# Merge lists in metric_lists
metric_df_lists3D_merged <- list()
metric_df_lists2D_merged <- list()

i <- 1
for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    for (metric in metrics) {
      if (i == 1)  {
        metric_df_lists3D_merged[[metric]] <- data.frame()
        metric_df_lists2D_merged[[metric]] <- data.frame()
      }
      if (i > 1) {
        temp <- nrow(metric_df_lists3D[[spes_metadata_index]][[metric]])
        n_slices <- length(unique(metric_df_lists2D[[spes_metadata_index]][[metric]][["slice"]]))
        if (metric %in% c("AMD", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC")) {
          metric_df_lists3D[[spes_metadata_index]][[metric]][["spe"]] <- 
            paste("spe", rep(seq((temp/4) * (i - 1) + 1, (temp/4) * (i - 1) + (temp/4)), each = 4), sep = "_")
          metric_df_lists2D[[spes_metadata_index]][[metric]][["spe"]] <-
            paste("spe", rep(seq((temp/4) * (i - 1) + 1, (temp)/4 * (i - 1) + (temp/4)), each = 4 * n_slices), sep = "_")
        }
        else {
          metric_df_lists3D[[spes_metadata_index]][[metric]][["spe"]] <- 
            paste("spe", rep(seq((temp/2) * (i - 1) + 1, (temp/2) * (i - 1) + (temp/2)), each = 2), sep = "_")
          metric_df_lists2D[[spes_metadata_index]][[metric]][["spe"]] <-
            paste("spe", rep(seq((temp/2) * (i - 1) + 1, (temp/2) * (i - 1) + (temp/2)), each = 2 * n_slices), sep = "_")
        }
      }
      metric_df_lists3D_merged[[metric]] <- rbind(metric_df_lists3D_merged[[metric]], metric_df_lists3D[[spes_metadata_index]][[metric]])
      metric_df_lists2D_merged[[metric]] <- rbind(metric_df_lists2D_merged[[metric]], metric_df_lists2D[[spes_metadata_index]][[metric]])
    }
    
    i <- i + 1
  }
}

metric_plots_3D_vs_2D_random_slice <- list()

for (metric in metrics) {
  metric_plots_3D_vs_2D_random_slice[[metric]] <- plot_3D_vs_2D_metric_random_slice_no_annotating_one_cell_type_pair(metric, 
                                                                                                                     metric_df_lists3D_merged[[metric]],
                                                                                                                     metric_df_lists2D_merged[[metric]], 
                                                                                                                     plots_metadata)
}



# Put plots into a pdf
setwd("~/R/plots/S2")

pdf("plots2D_vs_3D_random_slice_one_cell_pair.pdf", width = 12, height = 12)

metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")
curr_metric_plots <- list()
for (metric in metrics) {
  curr_metric_plots[[metric]] <- metric_plots_3D_vs_2D_random_slice[[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
}

plot <- plot_grid(plotlist = curr_metric_plots,
                  nrow = 4, 
                  ncol = 4)

print(plot)

dev.off()






### Get plots with 2D (all slices) on the x-axis and ERROR on the y-axis (not annotating for arrangement or shape and choosing random slice) (one cell type pair A/B or equivalent) ----------------

# Set up plots metadata
plots_metadata <- list(
  temp <- list(x_aes = "3D", y_aes = "error")
)

# Generate plots and plots into a list
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")


# Merge lists in metric_lists
metric_df_lists3D_merged <- list()
metric_df_lists2D_merged <- list()

i <- 1
for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    for (metric in metrics) {
      if (i == 1)  {
        metric_df_lists3D_merged[[metric]] <- data.frame()
        metric_df_lists2D_merged[[metric]] <- data.frame()
      }
      if (i > 1) {
        temp <- nrow(metric_df_lists3D[[spes_metadata_index]][[metric]])
        n_slices <- length(unique(metric_df_lists2D[[spes_metadata_index]][[metric]][["slice"]]))
        if (metric %in% c("AMD", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC")) {
          metric_df_lists3D[[spes_metadata_index]][[metric]][["spe"]] <- 
            paste("spe", rep(seq((temp/4) * (i - 1) + 1, (temp/4) * (i - 1) + (temp/4)), each = 4), sep = "_")
          metric_df_lists2D[[spes_metadata_index]][[metric]][["spe"]] <-
            paste("spe", rep(seq((temp/4) * (i - 1) + 1, (temp)/4 * (i - 1) + (temp/4)), each = 4 * n_slices), sep = "_")
        }
        else {
          metric_df_lists3D[[spes_metadata_index]][[metric]][["spe"]] <- 
            paste("spe", rep(seq((temp/2) * (i - 1) + 1, (temp/2) * (i - 1) + (temp/2)), each = 2), sep = "_")
          metric_df_lists2D[[spes_metadata_index]][[metric]][["spe"]] <-
            paste("spe", rep(seq((temp/2) * (i - 1) + 1, (temp/2) * (i - 1) + (temp/2)), each = 2 * n_slices), sep = "_")
        }
      }
      metric_df_lists3D_merged[[metric]] <- rbind(metric_df_lists3D_merged[[metric]], metric_df_lists3D[[spes_metadata_index]][[metric]])
      metric_df_lists2D_merged[[metric]] <- rbind(metric_df_lists2D_merged[[metric]], metric_df_lists2D[[spes_metadata_index]][[metric]])
    }
    
    i <- i + 1
  }
}

metric_plots_error_vs_2D_random_slice <- list()

for (metric in metrics) {
  metric_plots_error_vs_2D_random_slice[[metric]] <- plot_error_vs_2D_metric_random_slice_no_annotating_one_cell_type_pair(metric, 
                                                                                                                           metric_df_lists3D_merged[[metric]],
                                                                                                                           metric_df_lists2D_merged[[metric]], 
                                                                                                                           plots_metadata)
}



# Put plots into a pdf
setwd("~/R/plots/S2")

pdf("plots_error_vs_3D_random_slice_one_cell_pair.pdf", width = 12, height = 12)

metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")
curr_metric_plots <- list()
for (metric in metrics) {
  curr_metric_plots[[metric]] <- metric_plots_error_vs_2D_random_slice[[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
}

plot <- plot_grid(plotlist = curr_metric_plots,
                  nrow = 4, 
                  ncol = 4)

print(plot)

dev.off()





### Get plots with ERROR box plots for each metric (not annotating for arrangement or shape and choosing random slice) (one cell type pair A/B or equivalent) ----------------

# Generate plots and plots into a list
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC", "PBSAC", "PBP_AUC", "EBSAC", "EBP_AUC")


# Merge lists in metric_lists
metric_df_lists3D_merged <- list()
metric_df_lists2D_merged <- list()

i <- 1
for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    for (metric in metrics) {
      if (i == 1)  {
        metric_df_lists3D_merged[[metric]] <- data.frame()
        metric_df_lists2D_merged[[metric]] <- data.frame()
      }
      if (i > 1) {
        temp <- nrow(metric_df_lists3D[[spes_metadata_index]][[metric]])
        n_slices <- length(unique(metric_df_lists2D[[spes_metadata_index]][[metric]][["slice"]]))
        if (metric %in% c("AMD", "ACIN_AUC", "CKR_AUC", "CLR_AUC", "COO_AUC", "CGR_AUC")) {
          metric_df_lists3D[[spes_metadata_index]][[metric]][["spe"]] <- 
            paste("spe", rep(seq((temp/4) * (i - 1) + 1, (temp/4) * (i - 1) + (temp/4)), each = 4), sep = "_")
          metric_df_lists2D[[spes_metadata_index]][[metric]][["spe"]] <-
            paste("spe", rep(seq((temp/4) * (i - 1) + 1, (temp)/4 * (i - 1) + (temp/4)), each = 4 * n_slices), sep = "_")
        }
        else {
          metric_df_lists3D[[spes_metadata_index]][[metric]][["spe"]] <- 
            paste("spe", rep(seq((temp/2) * (i - 1) + 1, (temp/2) * (i - 1) + (temp/2)), each = 2), sep = "_")
          metric_df_lists2D[[spes_metadata_index]][[metric]][["spe"]] <-
            paste("spe", rep(seq((temp/2) * (i - 1) + 1, (temp/2) * (i - 1) + (temp/2)), each = 2 * n_slices), sep = "_")
        }
      }
      metric_df_lists3D_merged[[metric]] <- rbind(metric_df_lists3D_merged[[metric]], metric_df_lists3D[[spes_metadata_index]][[metric]])
      metric_df_lists2D_merged[[metric]] <- rbind(metric_df_lists2D_merged[[metric]], metric_df_lists2D[[spes_metadata_index]][[metric]])
    }
    
    i <- i + 1
  }
}

fig <- plot_error_box_plot_metric_random_slice_no_annotating_one_cell_type_pair(metrics, 
                                                                                metric_df_lists3D_merged, 
                                                                                metric_df_lists2D_merged)


# Put plots into a pdf
setwd("~/R/plots/S2")

pdf("plot_error_box_plot_metric_random_slice_no_annotating_one_cell_type_pair.pdf", width = 12, height = 6)

print(fig)

dev.off()


