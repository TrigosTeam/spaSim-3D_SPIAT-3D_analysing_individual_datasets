### 1. Set up all variable parameters for tables-----------------------------------------------

## Background
bg_prop_A_range <- c("min" = 0, "max" = 0.10)   # bg_prop_A ranges from 0 to 0.1
bg_prop_A_equals_zero_prob <- 1/2               # Give higher probability of bg_prop_A equals 0
bg_prop_B_range <- c("min" = 0, "max" = 0.10)   # bg_prop_B ranges from 0 to 0.1
bg_prop_B_equals_zero_prob <- 1/2               # Give higher probability of bg_prop_B equals 0

# Function to get bg props using bg_prop_range and probability prop equals 0
get_bg_props <- function(n_simulations, bg_prop_range, bg_prop_equals_zero_prob) {
  probs <- runif(n_simulations)
  props <- ifelse(probs <= bg_prop_equals_zero_prob, 0, round(runif(n_simulations, min = bg_prop_range["min"], max = bg_prop_range["max"]), 3))
  return(props)
}

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


### 2. Generate spes_table----------------------------------------------
n_simulations <- 1000
n_parameters <- 7 # The three ellipsoid radii count as 1 parameter
spes_table_colnames <- c("bg_prop_A", "bg_prop_B", "E_radius_x", "E_radius_y", "E_radius_z", "N_width", "cluster_prop_A", "ring_width_factor", "cluster_x_coord")
spes_table <- data.frame(matrix(nrow = n_simulations, ncol = length(spes_table_colnames))) 
colnames(spes_table) <- spes_table_colnames

# Fill in spes_table with variable parameter values
parameter_continuous_data <- list(
  bg_prop_A = get_bg_props(n_simulations, bg_prop_A_range, bg_prop_A_equals_zero_prob), # Make first bg_prop_A and bg_prop_B equal to 0
  bg_prop_B = get_bg_props(n_simulations, bg_prop_B_range, bg_prop_B_equals_zero_prob),
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
  spes_table[, parameter] <- parameter_continuous_data[[parameter]]
  i <- i + 1
}

setwd("~/R/spaSim-3D/scripts/simulations and analysis S2/S2 data")
write.table(spes_table, "spes_table.csv")


### 3. Generate sps metadata ----
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
setwd("~/R/spaSim-3D/scripts/simulations and analysis S2/S2 data")
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
    
    for (i in seq(nrow(spes_table))) {
      if (arrangement %in% c("mixed", "separated")) curr_metadata <- spe_metadata_cluster_template("regular", shape, bg_metadata)
      if (arrangement == "ringed") curr_metadata <- spe_metadata_cluster_template("ring", shape, bg_metadata)
      
      curr_metadata$background$cell_proportions <- c(spes_table$bg_prop_A[i], 
                                                     spes_table$bg_prop_B[i],
                                                     1 - spes_table$bg_prop_A[i] - spes_table$bg_prop_B[i]) # prop(O) = 1 - prop(A) - prop(B)
      
      
      if (shape == "ellipsoid") {
        curr_metadata$cluster_1$x_radius <- spes_table$E_radius_x[i]
        curr_metadata$cluster_1$y_radius <- spes_table$E_radius_y[i]
        curr_metadata$cluster_1$z_radius <- spes_table$E_radius_z[i]
        curr_metadata$cluster_1$x_y_rotation <- runif(1, 0, 180)
        curr_metadata$cluster_1$x_z_rotation <- runif(1, 0, 180)
        curr_metadata$cluster_1$y_z_rotation <- runif(1, 0, 180)
      }
      
      else if (shape == "network") {
        curr_metadata$cluster_1$n_edges <- N_n_edges
        curr_metadata$cluster_1$width <- spes_table$N_width[i]
        curr_metadata$cluster_1$radius <- N_radius
      }
      
      
      if (arrangement == "mixed") {
        curr_metadata$cluster_1$cluster_cell_types <- mixed_cluster_cell_types
        curr_metadata$cluster_1$cluster_cell_proportions <- c(spes_table$cluster_prop_A[i], 1 - spes_table$cluster_prop_A[i])
        curr_metadata$cluster_1$centre_loc <- mixed_cluster_centre_loc
      }
      else if (arrangement == "ringed") {
        curr_metadata$cluster_1$cluster_cell_types <- ringed_cluster_cell_type
        curr_metadata$cluster_1$cluster_cell_proportions <- ringed_cluster_cell_prop
        curr_metadata$cluster_1$centre_loc <- ringed_cluster_centre_loc
        curr_metadata$cluster_1$ring_cell_types <- ringed_ring_cell_type
        curr_metadata$cluster_1$ring_cell_proportions <- ringed_ring_cell_prop
        curr_metadata$cluster_1$ring_width <- spes_table$ring_width_factor[i] * ifelse(shape == "network", 
                                                                                       spes_table$N_width[i],
                                                                                       apply(spes_table[, c("E_radius_x", "E_radius_y", "E_radius_z")], 1, mean))
      }
      else if (arrangement == "separated") {
        curr_metadata$cluster_1$cluster_cell_types <- separated_cluster_cell_type
        curr_metadata$cluster_1$cluster_cell_proportions <- separated_cluster_cell_prop
        curr_metadata$cluster_1$centre_loc <- c(spes_table$cluster_x_coord[i], 300, 150)
        
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

setwd("~/R/spaSim-3D/scripts/simulations and analysis S2/S2 data")
saveRDS(spes_metadata, "spes_metadata.RDS")

