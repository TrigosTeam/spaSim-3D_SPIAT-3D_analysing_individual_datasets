### 1. Set up all fixed parameters for tables ------------------
## Background
bg_prop_A_fixed <- 0
bg_prop_B_fixed <- 0

## Ellipsoid radii
E_radius_x_fixed <- 75
E_radius_y_fixed <- 100
E_radius_z_fixed <- 125

## Network width
N_width_fixed <- 30

## Mixed clusters parameter
cluster_prop_A_fixed <- 0.7

## Ringed clusters parameter
ring_width_factor_fixed <- 0.15  # The width of the ring is proportional to the size (radius/width) of the cluster, for an ellipsoid, use the average of the three radii

## Separated clusters parameter
cluster_x_coord_fixed <- 150


### 2. Set up all variable parameters for tables-----------------------------------------------

## Background
bg_prop_A_range <- c("min" = 0, "max" = 0.10)   # bg_prop_A ranges from 0 to 0.1
bg_prop_B_range <- c("min" = 0, "max" = 0.10)   # bg_prop_B ranges from 0 to 0.1

## Ellipsoid radii
E_radius_x_range <- c("min" = 75, "max" = 125)
E_radius_y_range <- c("min" = 75, "max" = 125)
E_radius_z_range <- c("min" = 75, "max" = 125)

## Network width
N_width_range <- c("min" = 25, "max" = 35)


## Mixed clusters parameter
cluster_prop_A_range <- c("min" = 0.5, "max" = 0.9)   # Proportion of A cells in cluster will range from 0.5 to 0.9, proportion of B cells in cluster will be 1 - prop(B)


## Ringed clusters parameter
ring_width_factor_range <- c("min" = 0.1, "max" = 0.2)  # The width of the ring is proportional to the size (radius/width) of the cluster, for an ellipsoid, use the average of the three radii


## Separated clusters parameter
cluster_x_coord_range <- c("min" = 125, "max" = 175)


### 3. Generate spes_table----------------------------------------------
n_simulations <- 200
n_parameters <- 7 # The three ellipsoid radii count as 1 parameter
spes_table_colnames <- c("bg_prop_A", "bg_prop_B", "E_radius_x", "E_radius_y", "E_radius_z", "N_width", "cluster_prop_A", "ring_width_factor", "cluster_x_coord")
spes_table <- data.frame(matrix(nrow = n_simulations * n_parameters, ncol = length(spes_table_colnames))) 
colnames(spes_table) <- spes_table_colnames

# Fill in spes_table with fixed parameter values
spes_table$bg_prop_A <- bg_prop_A_fixed
spes_table$bg_prop_B <- bg_prop_B_fixed
spes_table$E_radius_x <- E_radius_x_fixed
spes_table$E_radius_y <- E_radius_y_fixed
spes_table$E_radius_z <- E_radius_z_fixed
spes_table$N_width <- N_width_fixed
spes_table$cluster_prop_A <- cluster_prop_A_fixed
spes_table$ring_width_factor <- ring_width_factor_fixed
spes_table$cluster_x_coord <- cluster_x_coord_fixed
spes_table$variable_parameter <- ""

# For every n_simulations rows, change one of the parameters to have a continuous range
parameter_continuous_data <- list(
  bg_prop_A = c(0 ,runif(n_simulations - 1, bg_prop_A_range["min"], bg_prop_A_range["max"])), # Make first bg_prop_A and bg_prop_B equal to 0
  bg_prop_B = c(0, runif(n_simulations - 1, bg_prop_B_range["min"], bg_prop_B_range["max"])),
  E_radius_x = runif(n_simulations, E_radius_x_range["min"], E_radius_x_range["max"]),
  E_radius_y = runif(n_simulations, E_radius_y_range["min"], E_radius_y_range["max"]),
  E_radius_z = runif(n_simulations, E_radius_z_range["min"], E_radius_z_range["max"]),
  N_width = runif(n_simulations, N_width_range["min"], N_width_range["max"]),
  cluster_prop_A = runif(n_simulations, cluster_prop_A_range["min"], cluster_prop_A_range["max"]),
  ring_width_factor = runif(n_simulations, ring_width_factor_range["min"], ring_width_factor_range["max"]),
  cluster_x_coord = runif(n_simulations, cluster_x_coord_range["min"], cluster_x_coord_range["max"])
)

i <- 1
for (parameter in spes_table_colnames) {
  if (parameter %in% c("E_radius_y", "E_radius_z")) i <- i - 1
  spes_table[((i - 1) * n_simulations + 1):(i * n_simulations), parameter] <- parameter_continuous_data[[parameter]]
  
  spes_table[((i - 1) * n_simulations + 1):(i * n_simulations), "variable_parameter"] <- parameter
  i <- i + 1
}

setwd("~/R/SPIAT-3D_benchmarking/simulations_and_analysis_S1/S1_data")
write.table(spes_table, "spes_table.csv")


### 4. Generate sps metadata ----
# constant metadata for simulations
bg_metadata <- spe_metadata_background_template("random")
bg_metadata$background$n_cells <- 30000
bg_metadata$background$length <- 600
bg_metadata$background$width <- 600
bg_metadata$background$height <- 300
bg_metadata$background$minimum_distance_between_cells <- 10
bg_metadata$background$cell_types <- c("A", "B", "O") # Cell proportions will change later

# Network
N_radius <- 125
N_n_edges <- 20

# Mixed clusters
mixed_cluster_cell_types <- c("A", "B")
mixed_cluster_centre_loc <- c(300, 300, 150)

# Ringed clusters
ringed_cluster_cell_type <- "A"
ringed_cluster_cell_prop <- 1
ringed_ring_cell_type <- "B"
ringed_ring_cell_prop <- 1
ringed_cluster_centre_loc <- c(300, 300, 150)

# Separated clusters
separated_cluster_cell_type <- "A"
separated_cluster_cell_prop <- 1




# Get table for simulations
setwd("~/R/SPIAT-3D_benchmarking/simulations_and_analysis_S1/S1_data")
spes_table <- read.table("spes_table.csv")

# Setup
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")

background_parameters <- c("bg_prop_A", "bg_prop_B")

arrangement_parameters <- list(mixed = "cluster_prop_A",
                               ringed = "ring_width_factor",
                               separated = "cluster_x_coord")

shape_parameters <- list(ellipsoid = c("E_radius_x", "E_radius_y", "E_radius_z"),
                         network = c("N_width"))

# Set up metadata list
spes_metadata <- list(mixed_ellipsoid = list(),
                      mixed_network = list(),
                      ringed_ellipsoid = list(),
                      ringed_network = list(),
                      separated_ellipsoid = list(),
                      separated_network = list())

for (arrangement in arrangements) {
  
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    spes_table_subset <- spes_table[spes_table$variable_parameter %in% c(background_parameters, arrangement_parameters[[arrangement]], shape_parameters[[shape]]), 
                                     c(background_parameters, arrangement_parameters[[arrangement]], shape_parameters[[shape]])]
    
    for (i in seq(nrow(spes_table_subset))) {
      if (arrangement %in% c("mixed", "separated")) curr_metadata <- spe_metadata_cluster_template("regular", shape, bg_metadata)
      if (arrangement == "ringed") curr_metadata <- spe_metadata_cluster_template("ring", shape, bg_metadata)
      
      curr_metadata$background$cell_proportions <- c(spes_table_subset$bg_prop_A[i], 
                                                     spes_table_subset$bg_prop_B[i],
                                                     1 - spes_table_subset$bg_prop_A[i] - spes_table_subset$bg_prop_B[i]) # prop(O) = 1 - prop(A) - prop(B)
      
      
      if (shape == "ellipsoid") {
        curr_metadata$cluster_1$radii <- c(spes_table_subset$E_radius_x[i], spes_table_subset$E_radius_y[i], spes_table_subset$E_radius_z[i])
        curr_metadata$cluster_1$axes_rotation <- c(runif(1, 0, 180), runif(1, 0, 180), runif(1, 0, 180))
      }
      
      else if (shape == "network") {
        curr_metadata$cluster_1$n_edges <- N_n_edges
        curr_metadata$cluster_1$width <- spes_table_subset$N_width[i]
        curr_metadata$cluster_1$radius <- N_radius
      }
      
      
      if (arrangement == "mixed") {
        curr_metadata$cluster_1$cluster_cell_types <- mixed_cluster_cell_types
        curr_metadata$cluster_1$cluster_cell_proportions <- c(spes_table_subset$cluster_prop_A[i], 1 - spes_table_subset$cluster_prop_A[i])
        curr_metadata$cluster_1$centre_loc <- mixed_cluster_centre_loc
      }
      else if (arrangement == "ringed") {
        curr_metadata$cluster_1$cluster_cell_types <- ringed_cluster_cell_type
        curr_metadata$cluster_1$cluster_cell_proportions <- ringed_cluster_cell_prop
        curr_metadata$cluster_1$centre_loc <- ringed_cluster_centre_loc
        curr_metadata$cluster_1$ring_cell_types <- ringed_ring_cell_type
        curr_metadata$cluster_1$ring_cell_proportions <- ringed_ring_cell_prop
        curr_metadata$cluster_1$ring_width <- spes_table_subset$ring_width_factor[i] * ifelse(shape == "network", 
                                                                                       spes_table_subset$N_width[i],
                                                                                       apply(spes_table_subset[, c("E_radius_x", "E_radius_y", "E_radius_z")], 1, mean))
      }
      else if (arrangement == "separated") {
        curr_metadata$cluster_1$cluster_cell_types <- separated_cluster_cell_type
        curr_metadata$cluster_1$cluster_cell_proportions <- separated_cluster_cell_prop
        curr_metadata$cluster_1$centre_loc <- c(spes_table_subset$cluster_x_coord[i], 300, 150)
        
        curr_metadata <- spe_metadata_cluster_template("regular", "sphere", curr_metadata)
        curr_metadata$cluster_2$cluster_cell_types <- "B"
        curr_metadata$cluster_2$cluster_cell_proportions <- 1
        curr_metadata$cluster_2$centre_loc <- c(450, 300, 150)
        curr_metadata$cluster_2$radius <- 100
      }
      
      spes_metadata[[spes_metadata_index]][[i]] <- curr_metadata
    }
  }
}

setwd("~/R/SPIAT-3D_benchmarking/simulations_and_analysis_S1/S1_data")
saveRDS(spes_metadata, "spes_metadata.RDS")

