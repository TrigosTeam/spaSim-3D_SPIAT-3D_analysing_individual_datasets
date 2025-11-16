# Generate data frame containing values of important parameters

generate_random_parameters <- function(
  n_simulations,
  bg_prop_A_range,
  bg_prop_B_range,
  E_radius_x_range,
  E_radius_y_range,
  E_radius_z_range,
  N_width_range,
  cluster_prop_A_range,
  ring_width_factor_range,
  cluster_x_coord_range
  ) {
  
  parameters_df <- data.frame(
    bg_prop_A = runif(n_simulations, bg_prop_A_range["min"], bg_prop_A_range["max"]),
    bg_prop_B = runif(n_simulations, bg_prop_B_range["min"], bg_prop_B_range["max"]),
    E_radius_x = runif(n_simulations, E_radius_x_range["min"], E_radius_x_range["max"]),
    E_radius_y = runif(n_simulations, E_radius_y_range["min"], E_radius_y_range["max"]),
    E_radius_z = runif(n_simulations, E_radius_z_range["min"], E_radius_z_range["max"]),
    N_width = runif(n_simulations, N_width_range["min"], N_width_range["max"]),
    cluster_prop_A = runif(n_simulations, cluster_prop_A_range["min"], cluster_prop_A_range["max"]),
    ring_width_factor = runif(n_simulations, ring_width_factor_range["min"], ring_width_factor_range["max"]),
    cluster_x_coord = runif(n_simulations, cluster_x_coord_range["min"], cluster_x_coord_range["max"])
  )
  
  return(parameters_df)
}

parameters_df <- generate_random_parameters(
  n_simulations = 10000,
  bg_prop_A_range = c("min" = 0, "max" = 0.10),
  bg_prop_B_range = c("min" = 0, "max" = 0.10),
  E_radius_x_range = c("min" = 75, "max" = 125),
  E_radius_y_range = c("min" = 75, "max" = 125),
  E_radius_z_range = c("min" = 75, "max" = 125),
  N_width_range = c("min" = 25, "max" = 35),
  cluster_prop_A_range = c("min" = 0.5, "max" = 0.9),
  ring_width_factor_range = c("min" = 0.1, "max" = 0.2) ,
  cluster_x_coord_range = c("min" = 125, "max" = 175)
)
