### Read updated data --------------------------------------------------------

setwd("~/Lin et al - human colorectal cancer/CRC1_data_updated/raw_slice_data")

slice_nums <- c("002", "007", "014", "020", "025", "029",
                "034", "039", "044", "049", "050", "051",
                "052", "054", "059", "064", "069", "074",
                "078", "084", "086", "091", "097")
file_names <- paste("CRC01-", slice_nums, ".csv", sep = "")

# Get z coord of each slice in micrometers using the slice number
get_slice_z_coord <- function(slice_num) {
  
  slice_num <- as.integer(slice_num)
  
  if (slice_num <= 85) {
    return(slice_num * 5)
  }
  
  else {
    return(85 * 5 + (slice_num - 85) * 4)
  }
}

# Store combined result in data frame
df <- data.frame()

# Loop through each slice
for (i in seq(length(slice_nums))) {
  
  file_name <- file_names[i]
  
  CRC1_df <- read.csv(file_name, sep = ",")
  CRC1_df <- CRC1_df[ , c("NewType", "Xtt", "Ytt")]
  CRC1_df$Ztt <- get_slice_z_coord(slice_nums[i])
  
  df <- rbind(df, CRC1_df)

}

# Current data frame has a number to represent each cell type
# Use the cell type dictionary to determine the cell type from the number

setwd("~/Lin et al - human colorectal cancer/Other_data")

cell_type_dictionary <- read.csv("cell_type_dictionary.csv")

specific_cell_types <- unique(cell_type_dictionary$Type_Name) # E.g. Tumour/Epi, Ki67+ Tumour/Epi...
df$Cell.Type.Specific <- specific_cell_types[df$NewType]

generic_cell_types <- unique(cell_type_dictionary$Category) # E.g. Tumour, Stroma, Immune
df$Cell.Type.Generic <- generic_cell_types[ifelse(df$NewType <= 3, 1, ifelse(df$NewType <= 6, 2, 3))]

df$NewType <- NULL

# Rename coordinate columns
colnames(df)[1:3] <- c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")

# Save
# setwd("~/Lin et al - human colorectal cancer/CRC1_data_updated")
# saveRDS(df, "CRC1_df.rds")




### Plot data ----------------------------------------------------------------

setwd("~/Lin et al - human colorectal cancer/CRC1_data_updated")
df <- readRDS("CRC1_df.rds")

library(dplyr)
library(rgl)
plot_cells_rgl_3D <- function(df, n_cells, plot_cell_types = NULL, feature_colname) {
  
  ## Get all unique cell types
  cell_types <- c("Tumor/Epi.", "Ki67+ Tumor/Epi.", "PDL1+ Tumor/Epi.", 
                  
                  "Endothelial", "Muscle/Fibroblast",      
                  
                  "Macrophage(I)", "Macrophage(II)", "Macrophage(III)", "Macrophage(IV)", "PDL1+ Macrophage", 
                  "PDL1+ lymphocyte",  "DN Lymphocyte", "DP Lymphocyte", "Lymphocyte(III)",    
                  "T helper", "PD1+ T helper", "Tc cell", "PD1+ Tc", "Treg", 
                  "B cells",
                  
                  "Other")
  
  ## Assign colour to each cell type
  cell_colours <- c("orange", "orange2", "orange3",       # Tumour
                    
                    "brown1", "brown",             # Endothelial & fibroblast
                    
                    "green1", "green", "green3", "green4", "darkgreen", # Macrophages
                    "purple1",  "purple2", "purple3", "purple4",    # Lymphocytes
                    "turquoise1", "turquoise2", "steelblue1", "steelblue2", "steelblue3", # T cells
                    "orchid",                 # B cells
                    
                    "lightgray")
  
  names(cell_colours) <- cell_types
  
  if (!is.null(plot_cell_types)) df <- df[df[[feature_colname]] %in% plot_cell_types, ]
  
  df_plot <- sample_n(df, n_cells)  
  df_plot$Cell.Colour <- cell_colours[df_plot[[feature_colname]]]
  
  options(rgl.printRglwidget = T)
  
  open3d()
  plot3d(df_plot$Cell.X.Position,
         df_plot$Cell.Y.Position,
         df_plot$Cell.Z.Position,
         col = df_plot$Cell.Colour,
         size = 4,
         xlab = 'x',
         ylab = 'y',
         zlab = 'z',
         xlim = NULL,
         ylim = NULL,
         zlim = c(0, 500),
         forceClipregion = TRUE)
  aspect3d(5, 5, 1)
  highlevel()
  
}

# Tumour cells
tumour_cell_types <- c("Tumor/Epi.", "Ki67+ Tumor/Epi.", "PDL1+ Tumor/Epi.")

# Stromal cells
stromal_cell_types <- c("Endothelial", "Muscle/Fibroblast")

# Immune cells
immune_cell_types <- c("Macrophage(I)", "Macrophage(II)", "Macrophage(III)", "Macrophage(IV)", "PDL1+ Macrophage",
                       "PDL1+ lymphocyte",  "DN Lymphocyte", "DP Lymphocyte", "Lymphocyte(III)",
                       "T helper", "PD1+ T helper", "Tc cell", "PD1+ Tc", "Treg",
                       "B cells")


plot_cells_rgl_3D(df, 50000, NULL, "Cell.Type.Specific")
plot_cells_rgl_3D(df, 50000, tumour_cell_types, "Cell.Type.Specific")
plot_cells_rgl_3D(df, 50000, stromal_cell_types, "Cell.Type.Specific")
plot_cells_rgl_3D(df, 50000, immune_cell_types, "Cell.Type.Specific")


plot_cells_rgl_3D(df, 50000, "Tumor/Epi.", "Cell.Type.Specific")


### calculate_minimum distance_df functions ------------------------------

calculate_minimum_distances_between_cell_types3D <- function(df,
                                                             cell_types_of_interest = NULL,
                                                             feature_colname = "Cell.Type",
                                                             show_summary = TRUE,
                                                             plot_image = TRUE) {
  
  if (is.null(df[["Cell.ID"]])) {
    warning("Temporarily adding Cell.Id column to your df")
    df$Cell.ID <- paste("Cell", seq(nrow(df)), sep = "_")
  }
  
  # Select all rows in data frame which only contains the cells of interest
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spe object
    unknown_cell_types <- setdiff(cell_types_of_interest, df[[feature_colname]])
    if (length(unknown_cell_types) != 0) {
      warning(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                    paste(unknown_cell_types, collapse = ", ")))
      return(data.frame())
    }
    
    df <- df[df[[feature_colname]] %in% cell_types_of_interest, ]
  }
  
  # If there are no cells, give error
  if (nrow(df) <= 1) {
    warning("There must be at least 2 cells of type 'cell_types_of_interest' in spe")
    return(data.frame())
  }
  
  # Create a list of the number of cell types with their corresponding cell ID's
  cell_types <- list()
  for (eachType in unique(df[[feature_colname]])) {
    cell_types[[eachType]] <- as.character(df$Cell.ID[df[[feature_colname]] == eachType])
  }
  
  # Get different possible cell type combinations
  # Each row represents a combination
  # If a row is [1 , 2], then we are comparing cell type 1 and cell type 2
  unique_cells <- unique(df[[feature_colname]]) # unique cell types
  permu <- gtools::permutations(length(unique_cells), 2, repeats.allowed = TRUE)
  
  result <- vector()
  
  for (i in seq(nrow(permu))) {
    name1 <- unique_cells[permu[i, 1]]
    name2 <- unique_cells[permu[i, 2]]
    
    # Get x,y,z coords for all cells of cell_type1 and cell_type2
    all_cell_type1_coord <- df[df[, feature_colname] == name1, 
                               c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
    
    all_cell_type2_coord <- df[df[, feature_colname] == name2, 
                               c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
    
    # Find all of closest points
    # For each cell of cell_type1, find the closest cell of cell_type2
    if (name1 != name2) {
      all_closest <- RANN::nn2(data = all_cell_type2_coord, 
                               query = all_cell_type1_coord, 
                               k = 1)  
    }
    # If we are comparing the same cell_type, and there is only one of this cell type, move on
    else if (nrow(all_cell_type1_coord) == 1) {
      warning("There is only 1 '", name1, "' cell in your data. It has no nearest neighbour of the same cell type.", sep = "")
      next
    }
    else {
      # If we are comparing the same cell_type, use the second closest neighbour
      all_closest <- RANN::nn2(data = all_cell_type2_coord, 
                               query = all_cell_type1_coord, 
                               k = 2)
      all_closest[['nn.idx']] <- all_closest[['nn.idx']][, 2]
      all_closest[['nn.dists']] <- all_closest[['nn.dists']][, 2]
    }
    
    # Create the data frame containing the chosen cells and their ids, as well as
    # the nearest cell to them and their ids, and the distance between
    cell_type2_cell_IDs <- df[df[ , feature_colname] == name2, "Cell.ID"]
    
    local_dist_mins <- data.frame(
      ref_cell_id = cell_types[[name1]],
      ref_cell_type = name1,
      nearest_cell_id = cell_type2_cell_IDs[as.vector(all_closest$nn.idx)],
      nearest_cell_type = name2,
      distance = all_closest$nn.dists
    )
    result <- rbind(result, local_dist_mins)
  }
  
  result$pair <- paste(result$ref_cell_type, result$nearest_cell_type,sep = "/")
  
  # Plot
  if (plot_image) {
    fig <- plot_cell_distances_violin3D(result)
    methods::show(fig)
  }
  
  # Print summary
  if (show_summary) {
    print(summarise_distances_between_cell_types3D(result))  
  }
  
  return(result)
}

summarise_distances_between_cell_types3D <- function(df) {
  
  pair <- distance <- NULL
  
  # summarise the results
  summarised_dists <- df %>% 
    dplyr::group_by(pair) %>%
    dplyr::summarise(mean(distance, na.rm = TRUE), 
                     min(distance, na.rm = TRUE), 
                     max(distance, na.rm = TRUE),
                     stats::median(distance, na.rm = TRUE), 
                     stats::sd(distance, na.rm = TRUE))
  
  summarised_dists <- data.frame(summarised_dists)
  
  colnames(summarised_dists) <- c("pair", 
                                  "mean", 
                                  "min", 
                                  "max", 
                                  "median", 
                                  "std_dev")
  
  for (i in seq(nrow(summarised_dists))) {
    # Get cell_types for each pair
    cell_types <- strsplit(summarised_dists[i,"pair"], "/")[[1]]
    
    summarised_dists[i, "reference"] <- cell_types[1]
    summarised_dists[i, "target"] <- cell_types[2]
  }
  
  return(summarised_dists)
}


## For scales parameter, use "free_x" or "free". "free_y" looks silly
plot_cell_distances_violin3D <- function(cell_to_cell_dist, scales = "free_x") {
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  pair <- distance <- NULL
  
  fig <- ggplot(cell_to_cell_dist, aes(x = pair, y = distance)) + 
    geom_violin() +
    facet_wrap(~pair, scales=scales, strip.position="bottom") +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
    labs(title="Cell distances", x = "Reference/Target pair", y = "Distance") +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult= 1), colour = "red")
  
  message("Plots show mean ± sd")
  
  return(fig)
}

### Get the minimum distance dfs between cell types within each slice -------------

# Read data
setwd("~/Lin et al - human colorectal cancer/CRC1_data_updated")
df <- readRDS("CRC1_df.rds")

# Get z-coords for each slice
slice_z_coords <- unique(df$Cell.Z.Position)

# Get minimum distances between each cell type for each slice
cell_types <- c("Tumor/Epi.", "Ki67+ Tumor/Epi.", "PDL1+ Tumor/Epi.", 
                "Endothelial", "Muscle/Fibroblast", "Macrophage(I)", 
                "Macrophage(II)", "Macrophage(III)", "Macrophage(IV)", "PDL1+ Macrophage",
                "PDL1+ lymphocyte",  "DN Lymphocyte", "DP Lymphocyte", "Lymphocyte(III)",
                "T helper", "PD1+ T helper", "Tc cell", "PD1+ Tc", "Treg",
                "B cells") # Excludes "Other cell type

minimum_distances_slice_df <- data.frame()
for (slice_z_coord in slice_z_coords) {
 
  # Get cells which have the current slice_z_coord
  slice_df <- df[df$Cell.Z.Position == slice_z_coord, ]
  
  for (cell_type in cell_types) {
    # Get minimum_distances_df for current slice and cell
    curr_minimum_distances_slice_df <- calculate_minimum_distances_between_cell_types3D(slice_df,
                                                                                        cell_types_of_interest = cell_type,
                                                                                        feature_colname = "Cell.Type.Specific",
                                                                                        show_summary = FALSE,
                                                                                        plot_image = FALSE) 
    
    if (nrow(curr_minimum_distances_slice_df) != 0) curr_minimum_distances_slice_df$slice_z_coord <- slice_z_coord
    minimum_distances_slice_df <- rbind(minimum_distances_slice_df, curr_minimum_distances_slice_df)
  }
}

# setwd("~/Lin et al - human colorectal cancer/CRC1_data_updated/minimum_distances_data")
# saveRDS(minimum_distances_slice_df, "minimum_distances_within_slices_df.rds")


### Get the minimum distance df between cell types after combining slices ---------------
# Read data
setwd("~/Lin et al - human colorectal cancer/CRC1_data_updated")
df <- readRDS("CRC1_df.rds")

# Get minimum distances between each cell type for each slice
cell_types <- c("Tumor/Epi.", "Ki67+ Tumor/Epi.", "PDL1+ Tumor/Epi.", 
                "Endothelial", "Muscle/Fibroblast", "Macrophage(I)", 
                "Macrophage(II)", "Macrophage(III)", "Macrophage(IV)", "PDL1+ Macrophage",
                "PDL1+ lymphocyte",  "DN Lymphocyte", "DP Lymphocyte", "Lymphocyte(III)",
                "T helper", "PD1+ T helper", "Tc cell", "PD1+ Tc", "Treg",
                "B cells") # Excludes "Other cell type

minimum_distances_slices_combined_df <- data.frame()


for (cell_type in cell_types) {
  # Get minimum_distances_df for current slice and cell
  curr_minimum_distances_slice_df <- calculate_minimum_distances_between_cell_types3D(df,
                                                                                      cell_types_of_interest = cell_type,
                                                                                      feature_colname = "Cell.Type.Specific",
                                                                                      show_summary = FALSE,
                                                                                      plot_image = FALSE) 
  
  minimum_distances_slices_combined_df <- rbind(minimum_distances_slices_combined_df, curr_minimum_distances_slice_df)
}
setwd("~/Lin et al - human colorectal cancer/objects/minimum_distances_data")
saveRDS(minimum_distances_slices_combined_df, "minimum_distance_all_slices_specific_df.rds")

### Get 'average_minimum_distance', 'lower_quantile_minimum_distance', 'average_shortest_5_percent_minimum_distance' from minimum distance dfs ----------------

# Read data
setwd("~/Lin et al - human colorectal cancer/CRC1_data_updated/minimum_distances_data")
df <- readRDS("minimum_distances_within_slices_df.rds")

# Get z-coords for each slice
slice_z_coords <- unique(df$slice_z_coord)

# Get minimum distances between each cell type for each slice
cell_types <- c("Tumor/Epi.", "Ki67+ Tumor/Epi.", "PDL1+ Tumor/Epi.", 
                "Endothelial", "Muscle/Fibroblast", "Macrophage(I)", 
                "Macrophage(II)", "Macrophage(III)", "Macrophage(IV)", "PDL1+ Macrophage",
                "PDL1+ lymphocyte",  "DN Lymphocyte", "DP Lymphocyte", "Lymphocyte(III)",
                "T helper", "PD1+ T helper", "Tc cell", "PD1+ Tc", "Treg",
                "B cells") # Excludes "Other cell type

### Define output data frames
average_minimum_distance_within_slices_df <- data.frame(matrix(nrow = length(slice_z_coords), ncol = length(cell_types)))
rownames(average_minimum_distance_within_slices_df) <- slice_z_coords
colnames(average_minimum_distance_within_slices_df) <- cell_types

lower_quantile_minimum_distance_within_slices_df <- data.frame(matrix(nrow = length(slice_z_coords), ncol = length(cell_types)))
rownames(lower_quantile_minimum_distance_within_slices_df) <- slice_z_coords
colnames(lower_quantile_minimum_distance_within_slices_df) <- cell_types

average_shortest_5_percent_minimum_distance_within_slices_df <- data.frame(matrix(nrow = length(slice_z_coords), ncol = length(cell_types)))
rownames(average_shortest_5_percent_minimum_distance_within_slices_df) <- slice_z_coords
colnames(average_shortest_5_percent_minimum_distance_within_slices_df) <- cell_types

for (slice_z_coord in slice_z_coords) {
  
  curr_slice_df <- df[df$slice_z_coord == slice_z_coord, ]
  for (cell_type in cell_types) {
    curr_slice_cell_type_df <- curr_slice_df[curr_slice_df$ref_cell_type == cell_type, ]
    
    # Fill in each output data frame
    average_minimum_distance_within_slices_df[as.character(slice_z_coord), cell_type] <- mean(curr_slice_cell_type_df$distance)
    lower_quantile_minimum_distance_within_slices_df[as.character(slice_z_coord), cell_type] <- quantile(curr_slice_cell_type_df$distance, 0.25)
    average_shortest_5_percent_minimum_distance_within_slices_df[as.character(slice_z_coord), cell_type] <- mean(curr_slice_cell_type_df$distance[curr_slice_cell_type_df$distance <= quantile(curr_slice_cell_type_df$distance, 0.05)])
  }  
}

setwd("~/Lin et al - human colorectal cancer/CRC1_data_updated/minimum_distances_data")
saveRDS(average_minimum_distance_within_slices_df, "average_minimum_distance_within_slices_df.rds")
saveRDS(lower_quantile_minimum_distance_within_slices_df, "lower_quantile_minimum_distance_within_slices_df.rds")
saveRDS(average_shortest_5_percent_minimum_distance_within_slices_df, "average_shortest_5_percent_minimum_distance_within_slices_df.rds")




### Plot 'average_minimum_distance', 'lower_quantile_minimum_distance', 'average_5_percent_shortest_distance' dfs ----------------------

# Read data
setwd("~/Lin et al - human colorectal cancer/CRC1_data_updated/minimum_distances_data")
average_minimum_distance_within_slices_df <- readRDS("average_minimum_distance_within_slices_df.rds")
lower_quantile_minimum_distance_within_slices_df <- readRDS("lower_quantile_minimum_distance_within_slices_df.rds")
average_shortest_5_percent_minimum_distance_within_slices_df <- readRDS("average_shortest_5_percent_minimum_distance_within_slices_df.rds")

library(ggplot2)
# function for number of observations 
median.n <- function(x) {
  return(c(y = 5, label = round(median(x), 2)))
  # experiment with the multiplier to find the perfect position
}

# function for mean labels
mean.n <- function(x) {
  return(c(y = 4.75, label = round(mean(x), 2))) 
  # experiment with the multiplier to find the perfect position
}

plot_minimum_distances_within_slices <- function(minimum_distances_within_slices_df, minimum_distance_metric) {
  
  df <- reshape2::melt(minimum_distances_within_slices_df)
  fig <- ggplot(df, aes(variable, value, color = variable)) +
    geom_boxplot(outliers = FALSE) +
    # facet_wrap(~variable, scale = "free_x") +
    stat_summary(fun.data = median.n, geom = "text", fun.y = median, colour = "black") +
    stat_summary(fun.data = mean.n, geom = "text", fun.y = mean, colour = "red") +
    labs(x = "cell type", y = minimum_distance_metric) +
    scale_color_discrete(name = "cell type") +
    theme_bw()
  methods::show(fig)
}

plot_minimum_distances_within_slices(average_minimum_distance_within_slices_df, "average minimum distance")
plot_minimum_distances_within_slices(lower_quantile_minimum_distance_within_slices_df, "lower quantile minimum distance")
plot_minimum_distances_within_slices(average_shortest_5_percent_minimum_distance_within_slices_df, "average shortest 5 percent minimum distance")


## Plotting just the averages
cell_types <- c("Tumor/Epi.", "Ki67+ Tumor/Epi.", "PDL1+ Tumor/Epi.", 
                "Endothelial", "Muscle/Fibroblast", "Macrophage(I)", 
                "Macrophage(II)", "Macrophage(III)", "Macrophage(IV)", "PDL1+ Macrophage",
                "PDL1+ lymphocyte",  "DN Lymphocyte", "DP Lymphocyte", "Lymphocyte(III)",
                "T helper", "PD1+ T helper", "Tc cell", "PD1+ Tc", "Treg",
                "B cells") # Excludes "Other cell type

minimum_distance_metrics_slice_averages <- data.frame(AMD = apply(average_minimum_distance_within_slices_df, 2, mean, na.rm = TRUE),
                                                      LQMD = apply(lower_quantile_minimum_distance_within_slices_df, 2, mean, na.rm = TRUE),
                                                      AS_5_PMD = apply(average_shortest_5_percent_minimum_distance_within_slices_df, 2, mean, na.rm = TRUE))
setwd("~/Lin et al - human colorectal cancer/CRC1_data_updated")
df <- readRDS("CRC1_df.rds")
minimum_distance_metrics_slice_averages$log_n_cells <- log(as.numeric(table(df$Cell.Type.Specific)[cell_types]), 2)
minimum_distance_metrics_slice_averages$Cell.Type <- rownames(minimum_distance_metrics_slice_averages)

plot_df <- reshape2::melt(minimum_distance_metrics_slice_averages, "Cell.Type")

plot_df$Cell.Type <- factor(plot_df$Cell.Type, cell_types)
ggplot(plot_df, aes(Cell.Type, value, color = Cell.Type)) + 
  geom_point(size = 5) +
  facet_wrap(~variable, nrow = 4, scales = "free_y") +
  labs(x = "cell type", y = "minimum distance value") +
  scale_color_discrete(name = "cell type") +
  theme_bw()


### Plot violin plots for minimum distances between cell types within each slice ----------------------

# Read minimum distance df for each slice (altogether)
setwd("~/Lin et al - human colorectal cancer/objects/minimum_distances_data")
minimum_distances_within_slices_specific_df <- readRDS("minimum_distances_within_slices_specific_df.rds")
minimum_distances_within_slices_specific_df <- minimum_distances_within_slices_specific_df[ , c("ref_cell_type", "distance", "slice_z_coord")]

# Get each cell type
cell_types <- c("Tumor/Epi.", "Ki67+ Tumor/Epi.", "PDL1+ Tumor/Epi.", 
                "Endothelial", "Muscle/Fibroblast", "Macrophage(I)", 
                "Macrophage(II)", "Macrophage(III)", "Macrophage(IV)", "PDL1+ Macrophage",
                "PDL1+ lymphocyte",  "DN Lymphocyte", "DP Lymphocyte", "Lymphocyte(III)",
                "T helper", "PD1+ T helper", "Tc cell", "PD1+ Tc", "Treg",
                "B cells") # Excludes "Other cell type

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Assign color to each cell type
colors <- gg_color_hue(length(cell_types))
names(colors) <- cell_types

# Factor df by cell types
minimum_distances_within_slices_specific_df$ref_cell_type <- factor(minimum_distances_within_slices_specific_df$ref_cell_type, cell_types)

# Get z-coords for each slice
slice_z_coords <- unique(minimum_distances_within_slices_specific_df$slice_z_coord)

library(ggplot2)


#### Save plots into pdf

# Violin plot
setwd("~/Lin et al - human colorectal cancer/plots")
# pdf("minimum_distances_within_slices_specific_violin_plots.pdf")
for (slice_z_coord in slice_z_coords) {
  slice_df <- minimum_distances_within_slices_specific_df[minimum_distances_within_slices_specific_df$slice_z_coord == slice_z_coord, ]
    
  fig <- ggplot(slice_df, aes(ref_cell_type, distance, col = ref_cell_type, fill = ref_cell_type)) +
    geom_violin() + 
    labs(title = paste("slice z-coord:", slice_z_coord, "microns"), x = "") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
    facet_wrap(.~ref_cell_type, scales = "free")
  print(fig)
}
dev.off()


# Density plot
setwd("~/Lin et al - human colorectal cancer/plots")
pdf("minimum_distances_within_slices_specific_density_plots.pdf")
for (slice_z_coord in slice_z_coords) {
  slice_df <- minimum_distances_within_slices_specific_df[minimum_distances_within_slices_specific_df$slice_z_coord == slice_z_coord, ]
  
  fig <- ggplot(slice_df, aes(x = distance, col = ref_cell_type, fill = ref_cell_type)) +
    geom_density() + 
    labs(title = paste("slice z-coord:", slice_z_coord, "microns"), x = "") +
    scale_x_continuous(n.breaks = 3) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
    facet_wrap(.~ref_cell_type, scales = "free")
  print(fig)
}
dev.off()




### Plot violin plots for minimum distances between cell types for all slices -----------------------
setwd("~/Lin et al - human colorectal cancer/objects/minimum_distances_data")

minimum_distances_all_slices_specific <- readRDS("minimum_distance_all_slices_specific_df.rds")

# Get each cell type
cell_types <- c("Tumor/Epi.", "Ki67+ Tumor/Epi.", "PDL1+ Tumor/Epi.", 
                "Endothelial", "Muscle/Fibroblast", "Macrophage(I)", 
                "Macrophage(II)", "Macrophage(III)", "Macrophage(IV)", "PDL1+ Macrophage",
                "PDL1+ lymphocyte",  "DN Lymphocyte", "DP Lymphocyte", "Lymphocyte(III)",
                "T helper", "PD1+ T helper", "Tc cell", "PD1+ Tc", "Treg",
                "B cells") # Excludes "Other cell type

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Assign color to each cell type
colors <- gg_color_hue(length(cell_types))
names(colors) <- cell_types

# Factor df by cell types
minimum_distances_all_slices_specific$ref_cell_type <- factor(minimum_distances_all_slices_specific$ref_cell_type, cell_types)


library(ggplot2)


#### Save plots into pdf

# Violin plot
setwd("~/Lin et al - human colorectal cancer/plots")
pdf("minimum_distances_all_slices_specific_violin_plots.pdf")

fig <- ggplot(minimum_distances_all_slices_specific, aes(ref_cell_type, distance, col = ref_cell_type, fill = ref_cell_type)) +
  geom_violin() + 
  labs(x = "") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  facet_wrap(.~ref_cell_type, scales = "free")
print(fig)

dev.off()


# Density plot
setwd("~/Lin et al - human colorectal cancer/plots")
pdf("minimum_distances_all_slices_specific_density_plots.pdf")

fig <- ggplot(minimum_distances_all_slices_specific, aes(x = distance, col = ref_cell_type, fill = ref_cell_type)) +
  geom_density() + 
  labs(x = "") +
  scale_x_continuous(n.breaks = 3) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  facet_wrap(.~ref_cell_type, scales = "free")
print(fig)

dev.off()
