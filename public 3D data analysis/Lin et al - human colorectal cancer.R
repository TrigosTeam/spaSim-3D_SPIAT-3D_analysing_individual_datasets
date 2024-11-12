### Get cell type dictionary to help annotate cell types ----------------------
setwd("~/Lin et al - human colorectal cancer/Other data")

cell_type_dict <- read.csv("cell_type_dictionary.csv")

# format column names of dictionary
all_cell_markers <- c("Keratin", "Ki67", "CD3", "CD20", 
                      "CD45RO", "CD4", "CD8a", "CD68", 
                      "CD163", "FOXP3", "PD1", "PDL1", 
                      "CD31", "aSMA", "Desmin", "CD45")
colnames(cell_type_dict)[1:length(all_cell_markers)] <- all_cell_markers 

# Get phenotype for each row
get_phenotype <- function(num_vec) {

  temp_cell_markers <- c("Keratin", "Ki67", "CD3", "CD20", 
                        "CD45RO", "CD4", "CD8a", "CD68", 
                        "CD163", "FOXP3", "PD1", "PDL1", 
                        "CD31", "aSMA", "Desmin", "CD45")
  
  temp_cell_markers <- temp_cell_markers[num_vec == 1]
  phenotype <- paste(temp_cell_markers, collapse = ',')
  
  return(phenotype)
}

cell_type_dict$phenotype <- apply(cell_type_dict[, 1:length(all_cell_markers)], 1, get_phenotype)

### 1a. Combine each slice, then determine cell types ---------

library(SPIAT)

### Get all slice data
setwd("~/Lin et al - human colorectal cancer/CRC1 data")

## Read data
slice_nums <- c("002", "007", "014", "020", "025", "029",
                "034", "039", "044", "049", "050", "051",
                "052", "054", "059", "064", "069", "074",
                "078", "084", "086", "091", "097", "102", "106")
file_names <- paste("WD-76845-", slice_nums, ".csv", sep = "")

dataset <- vector(mode = "list", length = length(file_names))

for (i in seq(length(file_names))) {
  dataset[[i]] <- read.csv(file_names[i])
}



## Loop through each slice and combine all slice data first (also add the z coord)
combined_data <- data.frame()

for (i in seq(length(dataset))) {
  data <- dataset[[i]]
  data$Z <- as.integer(slice_nums[i]) * 5
  
  combined_data <- rbind(combined_data, data)
}

# Get chosen cell markers
all_cell_markers <- c("Keratin", "Ki67", "CD3", "CD20", 
                      "CD45RO", "CD4", "CD8a", "CD68", 
                      "CD163", "FOXP3", "PD1", "PDL1", 
                      "CD31", "aSMA", "Desmin", "CD45")

## Get intensity matrix
intensity_matrix <- t(combined_data[, all_cell_markers])
colnames(intensity_matrix) <- paste("Cell_", seq(ncol(intensity_matrix)), sep="")

## Get x, y and z coords
coord_x <- combined_data[ , "X"]
coord_y <- combined_data[ , "Y"]
coord_z <- combined_data[ , "Z"]

# Format into spe object
general_format_image <- format_image_to_spe(format = "general",
                                            intensity_matrix = intensity_matrix,
                                            phenotypes = "unknown",
                                            coord_x = coord_x,
                                            coord_y = coord_y)

# Predict phenotypes
predicted_image <- predict_phenotypes(spe_object = general_format_image,
                                      thresholds = NULL,
                                      tumour_marker = "Keratin",
                                      baseline_markers = all_cell_markers[all_cell_markers != "Keratin"],
                                      reference_phenotype = FALSE)

all_phenotypes <- predicted_image$Phenotype
phenotypes <- unique(all_phenotypes)

# Get cell types from phenotypes
cell_types1 <- ifelse(phenotypes %in% cell_type_dict$phenotype, cell_type_dict$Category, "Other") ## Tumour, immune, stromal
cell_types2 <- ifelse(phenotypes %in% cell_type_dict$phenotype, cell_type_dict$Type_Name, "Other") ## All cell types

all_cell_types1 <- cell_types1[match(all_phenotypes, phenotypes)]
all_cell_types2 <- cell_types2[match(all_phenotypes, phenotypes)]

# Get final data frame
df <- data.frame(Cell.X.Position = coord_x,
                 Cell.Y.Position = coord_y,
                 Cell.Z.Position = coord_z,
                 Cell.Type1 = all_cell_types1,
                 Cell.Type2 = all_cell_types2)

### 1b. Plot (tumour, immune, stromal) --------------------

setwd("~/Lin et al - human colorectal cancer/Other data")
df <- readRDS("lin_et_al_3D_spatial_data_combine_slices.rds")

## Get all unique cell types
cell_types <- c("Tumor", "Immune", "Stroma", "Other")

## Assign colour to each cell type
cell_colours <- c("orange", "skyblue2", "brown", "lightgray")

names(cell_colours) <- cell_types
df$Cell.Colour <- cell_colours[df$Cell.Type1]


## Plot
library(dplyr)
library(rgl)
options(rgl.printRglwidget = T)

df_plot <- sample_n(df, 100000)
plot3d(df_plot$Cell.X.Position,
       df_plot$Cell.Y.Position,
       df_plot$Cell.Z.Position,
       col = df_plot$Cell.Colour,
       size = 8,
       xlab = 'x',
       ylab = 'y',
       zlab = 'z')

# legend3d("topright", legend = cell_types, pch = 16, col = cell_colours, cex=0.5, inset=0.02)



### 1c. Plot (all cell types) --------------------

setwd("~/Lin et al - human colorectal cancer/Other data")
df <- readRDS("lin_et_al_3D_spatial_data_combine_slices.rds")

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
df$Cell.Colour <- cell_colours[df$Cell.Type2]


## Plot
library(dplyr)
library(rgl)
options(rgl.printRglwidget = T)

df_plot <- sample_n(df, 100000)
plot3d(df_plot$Cell.X.Position,
       df_plot$Cell.Y.Position,
       df_plot$Cell.Z.Position,
       col = df_plot$Cell.Colour,
       size = 8,
       xlab = 'x',
       ylab = 'y',
       zlab = 'z')



## Plot for a particular cell type
library(dplyr)
library(rgl)
options(rgl.printRglwidget = T)

# cells_chosen <- c("Tumor/Epi.", "Ki67+ Tumor/Epi.", "PDL1+ Tumor/Epi.")
# cells_chosen <- c("Endothelial", "Muscle/Fibroblast")
# cells_chosen <- c("Macrophage(I)", "Macrophage(II)", "Macrophage(III)", "Macrophage(IV)", "PDL1+ Macrophage",
#                   "PDL1+ lymphocyte",  "DN Lymphocyte", "DP Lymphocyte", "Lymphocyte(III)",
#                   "T helper", "PD1+ T helper", "Tc cell", "PD1+ Tc", "Treg",
#                   "B cells")

df_cells_chosen <- df[df$Cell.Type2 %in% cells_chosen, ]
df_plot <- sample_n(df_cells_chosen, 100000)
plot3d(df_plot$Cell.X.Position,
       df_plot$Cell.Y.Position,
       df_plot$Cell.Z.Position,
       col = df_plot$Cell.Colour,
       size = 8,
       xlab = 'x',
       ylab = 'y',
       zlab = 'z')


# legend3d("topright", legend = cell_types, pch = 16, col = cell_colours, cex=0.5, inset=0.02)








### 3a. Examine each slice separately to determine cell types


### 2a. Examine each slice separately and determine cell types for each slice ---------

library(SPIAT)

### Get all slice data
setwd("~/Lin et al - human colorectal cancer/CRC1 data")

## Read data
slice_nums <- c("002", "007", "014", "020", "025", "029",
                "034", "039", "044", "049", "050", "051",
                "052", "054", "059", "064", "069", "074",
                "078", "084", "086", "091", "097", "102", "106")
file_names <- paste("WD-76845-", slice_nums, ".csv", sep = "")

dataset <- vector(mode = "list", length = length(file_names))

for (i in seq(length(file_names))) {
  dataset[[i]] <- read.csv(file_names[i])
}



# Get chosen cell markers
all_cell_markers <- c("Keratin", "Ki67", "CD3", "CD20", 
                      "CD45RO", "CD4", "CD8a", "CD68", 
                      "CD163", "FOXP3", "PD1", "PDL1", 
                      "CD31", "aSMA", "Desmin", "CD45")


## Define result data.frame
df <- data.frame()


## Loop through each slice
i <- 1
for (data in dataset) {
  
  ## Get intensity matrix
  intensity_matrix <- t(data[, all_cell_markers])
  colnames(intensity_matrix) <- paste("Cell_", seq(ncol(intensity_matrix)), sep="")
  
  ## Get x, y and z coords
  coord_x <- data[ , "X"]
  coord_y <- data[ , "Y"]
  coord_z <- rep(as.integer(slice_nums[i]) * 5, nrow(data))
  
  # Format into spe object
  general_format_image <- format_image_to_spe(format = "general",
                                              intensity_matrix = intensity_matrix,
                                              phenotypes = "unknown",
                                              coord_x = coord_x,
                                              coord_y = coord_y)
  
  # Predict phenotypes
  predicted_image <- predict_phenotypes(spe_object = general_format_image,
                                        thresholds = NULL,
                                        tumour_marker = "Keratin",
                                        baseline_markers = all_cell_markers[all_cell_markers != "Keratin"],
                                        reference_phenotype = FALSE)
  
  all_phenotypes <- predicted_image$Phenotype
  phenotypes <- unique(all_phenotypes)
  
  # Get cell types from phenotypes
  cell_types1 <- ifelse(phenotypes %in% cell_type_dict$phenotype, cell_type_dict$Category, "Other") ## Tumour, immune, stromal
  cell_types2 <- ifelse(phenotypes %in% cell_type_dict$phenotype, cell_type_dict$Type_Name, "Other") ## All cell types
  
  all_cell_types1 <- cell_types1[match(all_phenotypes, phenotypes)]
  all_cell_types2 <- cell_types2[match(all_phenotypes, phenotypes)]
  
  df_temp <- data.frame(Cell.X.Position = coord_x,
                        Cell.Y.Position = coord_y,
                        Cell.Z.Position = coord_z,
                        Cell.Type1 = all_cell_types1,
                        Cell.Type2 = all_cell_types2)
  
  df <- rbind(df, df_temp)
  
  i <- i + 1
  
  break
}



### 2b. Plot (tumour, immune, stromal) --------------------

setwd("~/Lin et al - human colorectal cancer/Other data")
df <- readRDS("lin_et_al_3D_spatial_data_separate_slices.rds")

## Get all unique cell types
cell_types <- c("Tumor", "Immune", "Stroma", "Other")

## Assign colour to each cell type
cell_colours <- c("orange", "skyblue2", "brown", "lightgray")

names(cell_colours) <- cell_types
df$Cell.Colour <- cell_colours[df$Cell.Type1]


## Plot
library(dplyr)
library(rgl)
options(rgl.printRglwidget = T)

df_plot <- sample_n(df, 100000)
plot3d(df_plot$Cell.X.Position,
       df_plot$Cell.Y.Position,
       df_plot$Cell.Z.Position,
       col = df_plot$Cell.Colour,
       size = 8,
       xlab = 'x',
       ylab = 'y',
       zlab = 'z')

# legend3d("topright", legend = cell_types, pch = 16, col = cell_colours, cex=0.5, inset=0.02)



### 2c. Plot (all cell types) --------------------

setwd("~/Lin et al - human colorectal cancer/Other data")
df <- readRDS("lin_et_al_3D_spatial_data_separate_slices.rds")

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
df$Cell.Colour <- cell_colours[df$Cell.Type2]


## Plot
library(dplyr)
library(rgl)
options(rgl.printRglwidget = T)

df_plot <- sample_n(df, 100000)
plot3d(df_plot$Cell.X.Position,
       df_plot$Cell.Y.Position,
       df_plot$Cell.Z.Position,
       col = df_plot$Cell.Colour,
       size = 8,
       xlab = 'x',
       ylab = 'y',
       zlab = 'z')


## Plot for a particular cell type
library(dplyr)
library(rgl)
options(rgl.printRglwidget = T)

# cells_chosen <- c("Tumor/Epi.", "Ki67+ Tumor/Epi.", "PDL1+ Tumor/Epi.")
# cells_chosen <- c("Endothelial", "Muscle/Fibroblast")
cells_chosen <- c("Macrophage(I)", "Macrophage(II)", "Macrophage(III)", "Macrophage(IV)", "PDL1+ Macrophage",
                  "PDL1+ lymphocyte",  "DN Lymphocyte", "DP Lymphocyte", "Lymphocyte(III)",
                  "T helper", "PD1+ T helper", "Tc cell", "PD1+ Tc", "Treg",
                  "B cells")

df_cells_chosen <- df[df$Cell.Type2 %in% cells_chosen, ]
df_plot <- sample_n(df_cells_chosen, 100000)
 
df_plot$Cell.Type2 <- factor(df_plot$Cell.Type2, rev(cells_chosen))

plot3d(df_plot$Cell.X.Position,
       df_plot$Cell.Y.Position,
       df_plot$Cell.Z.Position,
       col = df_plot$Cell.Colour,
       size = 4,
       xlab = 'x',
       ylab = 'y',
       zlab = 'z',
       alpha = 0.5)

# legend3d("topright", legend = cell_types, pch = 16, col = cell_colours, cex=0.5, inset=0.02)



### 3a. Combine each slice and get thresholds ---------------------------------
library(SPIAT)

### Get all slice data
setwd("~/Lin et al - human colorectal cancer/CRC1 data")

## Read data
slice_nums <- c("002", "007", "014", "020", "025", "029",
                "034", "039", "044", "049", "050", "051",
                "052", "054", "059", "064", "069", "074",
                "078", "084", "086", "091", "097", "102", "106")
file_names <- paste("WD-76845-", slice_nums, ".csv", sep = "")

dataset <- vector(mode = "list", length = length(file_names))

for (i in seq(length(file_names))) {
  dataset[[i]] <- read.csv(file_names[i])
}



## Loop through each slice and combine all slice data first (also add the z coord)
combined_data <- data.frame()

for (i in seq(length(dataset))) {
  data <- dataset[[i]]
  
  combined_data <- rbind(combined_data, data)
}

# Get chosen cell markers
all_cell_markers <- c("Keratin", "Ki67", "CD3", "CD20", 
                      "CD45RO", "CD4", "CD8a", "CD68", 
                      "CD163", "FOXP3", "PD1", "PDL1", 
                      "CD31", "aSMA", "Desmin", "CD45")

## Get intensity matrix
intensity_matrix <- t(combined_data[, all_cell_markers])
colnames(intensity_matrix) <- paste("Cell_", seq(ncol(intensity_matrix)), sep="")

## Get x, y coords
coord_x <- combined_data[ , "X"]
coord_y <- combined_data[ , "Y"]

# Format into spe object
general_format_image <- format_image_to_spe(format = "general",
                                            intensity_matrix = intensity_matrix,
                                            phenotypes = "unknown",
                                            coord_x = coord_x,
                                            coord_y = coord_y)

# Get thresholds
thresholds_df <- get_thresholds(spe_object = general_format_image,
                                thresholds = NULL,
                                tumour_marker = "Keratin",
                                baseline_markers = all_cell_markers[all_cell_markers != "Keratin"],
                                reference_phenotype = FALSE)

thresholds_df <- data.frame(thresholds_df)
rownames(thresholds_df) <- 1


setwd("~/Lin et al - human colorectal cancer/Other data")
write.table(thresholds_df, "thresholds_combine_slices_df.csv")




### 3b. Examine each slice separately and get thresholds for each slice ------
library(SPIAT)

### Get all slice data
setwd("~/Lin et al - human colorectal cancer/CRC1 data")

## Read data
slice_nums <- c("002", "007", "014", "020", "025", "029",
                "034", "039", "044", "049", "050", "051",
                "052", "054", "059", "064", "069", "074",
                "078", "084", "086", "091", "097", "102", "106")
file_names <- paste("WD-76845-", slice_nums, ".csv", sep = "")

dataset <- vector(mode = "list", length = length(file_names))

for (i in seq(length(file_names))) {
  dataset[[i]] <- read.csv(file_names[i])
}



# Get chosen cell markers
all_cell_markers <- c("Keratin", "Ki67", "CD3", "CD20", 
                      "CD45RO", "CD4", "CD8a", "CD68", 
                      "CD163", "FOXP3", "PD1", "PDL1", 
                      "CD31", "aSMA", "Desmin", "CD45")


## Define threshold data frame
thresholds_df <- data.frame()


## Loop through each slice
i <- 1
for (data in dataset) {
  
  ## Get intensity matrix
  intensity_matrix <- t(data[, all_cell_markers])
  colnames(intensity_matrix) <- paste("Cell_", seq(ncol(intensity_matrix)), sep="")
  
  ## Get x
  coord_x <- data[ , "X"]
  coord_y <- data[ , "Y"]
  
  # Format into spe object
  general_format_image <- format_image_to_spe(format = "general",
                                              intensity_matrix = intensity_matrix,
                                              phenotypes = "unknown",
                                              coord_x = coord_x,
                                              coord_y = coord_y)
  
  # Get thresholds
  thresholds <- get_thresholds(spe_object = general_format_image,
                               thresholds = NULL,
                               tumour_marker = "Keratin",
                               baseline_markers = all_cell_markers[all_cell_markers != "Keratin"],
                               reference_phenotype = FALSE)
  
  # Add to thresholds data.frame
  thresholds_df <- rbind(thresholds_df, data.frame(thresholds))
  
  i <- i + 1
}
rownames(thresholds_df) <- seq(nrow(thresholds_df))

# setwd("~/Lin et al - human colorectal cancer/Other data")
# write.table(thresholds_df, "thresholds_separate_slices_df.csv")


### 3c. Plot thresholds (only separate slices) -------------------------------------------------------
thresholds_df <- read.table("thresholds_separate_slices_df.csv")
thresholds_df$slice <- seq(nrow(thresholds_df))
thresholds_plot_df <- reshape2::melt(thresholds_df, "slice")

fig <- ggplot(thresholds_plot_df, aes(slice, value, color = variable)) + 
  geom_line() +
  facet_wrap(~variable, nrow = 4, ncol = 4, scales = "free") +
  theme_bw()
methods::show(fig)


### 3d. Plot thresholds (combine slices and separate slices) --------------------------
thresholds_combine_slices_df <- read.table("thresholds_combine_slices_df.csv")
thresholds_separate_slices_df <- read.table("thresholds_separate_slices_df.csv")

thresholds_df <- rbind(thresholds_separate_slices_df, thresholds_combine_slices_df)

thresholds_df$slice <- seq(nrow(thresholds_df))
thresholds_df[nrow(thresholds_df), "slice"] <- -1
thresholds_df$key <- ifelse(thresholds_df$slice == -1, "combine", "separate")


thresholds_plot_df <- reshape2::melt(thresholds_df, c("slice", "key"))

fig <- ggplot(data = thresholds_plot_df[thresholds_plot_df$key == "separate", ],
              aes(slice, value, color = variable)) +
  geom_line() +
  geom_point(data = thresholds_plot_df[thresholds_plot_df$key == "combine", ],
             aes(slice, value, color = variable),
             shape = 8) +
  facet_wrap(~variable, nrow = 4, ncol = 4, scales = "free") +
  theme_bw()

methods::show(fig)



### 3e. Get thresholds function --------------------------------------------

get_thresholds <- function(spe_object, thresholds = NULL, tumour_marker,
                           baseline_markers, nuclear_marker = NULL,
                           reference_phenotypes = FALSE, markers_to_phenotype = NULL,
                           plot_distribution=TRUE){
  
  Marker_level <- NULL
  
  formatted_data <- SummarizedExperiment::colData(spe_object)
  intensity_matrix <- SummarizedExperiment::assay(spe_object)
  
  if(is.null(markers_to_phenotype)){
    markers <- rownames(intensity_matrix)
  }else{
    markers <- markers_to_phenotype
    if (sum(markers %in% rownames(intensity_matrix)) != length(markers)) {
      missing <- markers[!(markers %in% rownames(intensity_matrix))]
      stop(sprintf("There are no intensity values for %s in your input dataset. Please check markers\n", 
                   missing))
    }
    intensity_matrix <- intensity_matrix[markers,]
  }
  
  #CHECK
  if (is.element(tumour_marker, markers) == FALSE) {
    stop("Tumour marker not found")
  }
  
  cell_ids <- colnames(intensity_matrix)
  
  rownames(intensity_matrix) <- NULL
  colnames(intensity_matrix) <- NULL
  intensity_matrix_t <- t(intensity_matrix)
  intensity_df <- data.frame(intensity_matrix_t)
  colnames(intensity_df) <- markers
  
  formatted_data <- cbind(formatted_data, intensity_df)
  
  #add actual intensity boolean value to formatted_data
  markers_no_tumour <- markers[markers != tumour_marker]
  if(!is.null(nuclear_marker)){
    markers_no_tumour <- markers_no_tumour[markers_no_tumour != nuclear_marker]
  }
  
  selected_valley_xcord <- list()
  
  ##Add actual marker levels
  for (marker in markers){
    if ((!is.null(nuclear_marker)) && (marker == nuclear_marker)){
      actual_phenotype <- data.frame(rep(1, nrow(formatted_data)))
      colnames(actual_phenotype) <- paste(marker, "_actual_phenotype", sep="")
      formatted_data <- cbind(formatted_data, actual_phenotype)
    } else{
      cell_IDs <- data.frame(formatted_data$Cell.ID)
      colnames(cell_IDs) <- "Cell.ID"
      actual_phenotype <- data.frame(rep(0, nrow(formatted_data)))
      colnames(actual_phenotype) <- "Phenotype_status"
      actual_phenotype <- cbind(actual_phenotype, cell_IDs)
      
      rows <- formatted_data[grepl(marker, formatted_data$Phenotype), ]
      if(nrow(rows) > 0){
        actual_phenotype[actual_phenotype$Cell.ID %in% rows$Cell.ID,]$Phenotype_status <- 1
        
        actual_phenotype <- data.frame(actual_phenotype[,1])
        colnames(actual_phenotype) <- paste(marker, "_actual_phenotype", sep="")
        formatted_data <- cbind(formatted_data, actual_phenotype)  
      }
    }
  }
  
  for (marker in markers_no_tumour) {
    
    #extract the marker intensity column
    marker_specific_level <- formatted_data[,marker]
    
    #calculate the predictions
    if (!is.null(thresholds) && !is.na(thresholds[match(marker,markers)])) {
      #there is a threshold value specified for the marker, use the threshold
      marker_threshold <- thresholds[match(marker,markers)]
      methods::show( paste0(marker, " has threshold specified: ",as.character(marker_threshold)))
      selected_valley_xcord[[marker]] <- marker_threshold
      
      #get the threshold predictions
      predictions_by_threshold <- data.frame(mmand::threshold(marker_specific_level, level = marker_threshold))
      
    } else {
      #calculate the valleys
      intensity_density <- stats::density(marker_specific_level, na.rm=TRUE)
      valleys <- pracma::findpeaks(-(intensity_density)$y)
      valley_ycords <- valleys[,1] * -1
      index <- match(valley_ycords, intensity_density$y)
      valley_xcords <- intensity_density$x[index]
      
      #create a df for the valley coordinates
      valley_df <- data.frame(cbind(valley_xcords, valley_ycords))
      
      #select the first valley that's greater than the maximum density and below 25% density
      ycord_max_density <- max(intensity_density$y)
      xcord_max_density_index <- match(ycord_max_density, intensity_density$y)
      xcord_max_density <- intensity_density$x[xcord_max_density_index]
      
      density_threshold_for_valley <- 0.25 * ycord_max_density
      
      valley_df <- valley_df[valley_df$valley_xcords >= xcord_max_density, ]
      valley_df <- valley_df[valley_df$valley_ycords <= density_threshold_for_valley, ]
      
      selected_valley_xcord[[marker]] <- valley_df$valley_xcords[1]
      #using the selected valley as the threshold
      predictions_by_threshold <- data.frame(mmand::threshold(marker_specific_level, level = selected_valley_xcord[[marker]]))
    }
    colnames(predictions_by_threshold) <- paste(marker, "_predicted_phenotype", sep="")
    formatted_data <- cbind(formatted_data, predictions_by_threshold)
  }
  
  
  ###Prediction for tumour marker
  #Select cells that are positive for a marker that tumor cells don't have
  
  baseline_cells <- vector()
  for(marker in baseline_markers){
    temp <- formatted_data[, colnames(formatted_data) %in% c("Cell.ID", paste0(marker, "_predicted_phenotype"))]
    temp <- temp[temp[,2] == 1,]
    baseline_cells <- c(baseline_cells, temp$Cell.ID)
  }
  baseline_cells <- unique(baseline_cells)
  
  #Tumor marker levels in these cells
  formatted_data_baseline <- formatted_data[formatted_data$Cell.ID %in% baseline_cells,tumour_marker]
  if(length(unique(formatted_data_baseline)) != 2){
    cutoff_for_tumour <- stats::quantile(formatted_data_baseline, 0.95)
  }else{
    cutoff_for_tumour <- max(formatted_data_baseline)-min(formatted_data_baseline)/2
  }
  
  #extract the marker intensity column
  tumour_specific_level <- formatted_data[,tumour_marker]
  
  #calculate the predictions
  if (!is.null(thresholds)) {
    #there is a threshold value specified for the marker, use the threshold
    marker_threshold <- thresholds[match(tumour_marker,markers)]
    methods::show(paste0(tumour_marker, " has threshold specified: ", as.character(marker_threshold)))
    selected_valley_xcord[[tumour_marker]] <- marker_threshold
    
    #get the threshold predictions
    predictions_by_threshold <- data.frame(mmand::threshold(tumour_specific_level, level = marker_threshold))
    
  } else {
    #calculate the valleys
    intensity_density <- stats::density(tumour_specific_level, na.rm=TRUE)
    valleys <- pracma::findpeaks(-(intensity_density)$y)
    valley_ycords <- valleys[,1] * -1
    index <- match(valley_ycords, intensity_density$y)
    valley_xcords <- intensity_density$x[index]
    
    #create a df for the valley coordinates
    valley_df <- data.frame(cbind(valley_xcords, valley_ycords))
    selected_valley_xcord[[tumour_marker]] <- valley_df$valley_xcords[1]
    
    #using the selected valley as the threshold if it is lower than the
    #level of intensity of the tumour marker in non-tumour cells
    
    final_threshold <- ifelse(selected_valley_xcord[[tumour_marker]] < cutoff_for_tumour,
                              selected_valley_xcord[[tumour_marker]], cutoff_for_tumour)
    selected_valley_xcord[[tumour_marker]] <- final_threshold
    predictions_by_threshold <- data.frame(mmand::threshold(tumour_specific_level, level = final_threshold))
  }
  
  return(selected_valley_xcord)
}

### 4a. Plot gene expression values for each slice ---------------------------
library(ggplot2)
library(dplyr)
library(cowplot)

setwd("~/Lin et al - human colorectal cancer/CRC1_data")
protein_expression_data <- read.table("protein_expression_data.csv")

## Melt
plot_df <- reshape2::melt(protein_expression_data,  c("Z", "slice"))

plot_sample_df <- sample_n(plot_df, 10000000)

# Plot with no outliers
fig <- ggplot(plot_sample_df, aes(Z, value, color = variable, group = Z)) +
  geom_boxplot(outlier.size = 0.01) +
  stat_summary(fun = median, geom = "line", aes(group = 1), color = "grey") +
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "black") +
  facet_wrap(~variable, nrow = 4, ncol = 4, scales = "free") +
  theme_bw()
methods::show(fig)



# Plot without outliers
all_cell_markers <- c("Keratin", "Ki67", "CD3", "CD20", 
                      "CD45RO", "CD4", "CD8a", "CD68", 
                      "CD163", "FOXP3", "PD1", "PDL1", 
                      "CD31", "aSMA", "Desmin", "CD45")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

colors <- gg_color_hue(length(all_cell_markers))
names(colors) <- all_cell_markers

plot_sample_df <- sample_n(plot_df, 10000000)
fig_list <- list()

for (cell_marker in all_cell_markers) {
  
  plot_sample_marker_df <- plot_sample_df[plot_sample_df$variable == cell_marker, ]
  
  y_max <- boxplot.stats(plot_sample_marker_df$value)[["stats"]][5]
  
  # y_max <- 0
  # for (slice in seq(25)) {
  #   plot_sample_marker_slice_df <- plot_sample_marker_df[plot_sample_marker_df$slice == slice, ]
  #   curr_y_max <- boxplot.stats(plot_sample_marker_slice_df$value)[["stats"]][5]
  #   if (curr_y_max > y_max) y_max <- curr_y_max
  #   
  # 
  #   View(plot_sample_marker_slice_df)
  #   break
  # }
  # 
  
  fig <- ggplot(plot_sample_marker_df, aes(Z, value, group = Z)) +
    geom_boxplot(outlier.shape = NA, color = colors[cell_marker]) +
    stat_summary(fun = median, geom = "line", aes(group = 1), color = "grey") +
    stat_summary(fun = mean, geom = "line", aes(group = 1), color = "black") +
    scale_y_continuous(limits = c(0, y_max)) +
    theme_bw() +
    labs(x = "z-coord", y = "", title = cell_marker) +
    theme(plot.title = element_text(hjust = 0.5))

  fig_list[[cell_marker]] <- fig    
  
}

plot_grid(fig_list$Keratin, fig_list$Ki67, fig_list$CD3, fig_list$CD20,
          fig_list$CD45RO, fig_list$CD4, fig_list$CD8a, fig_list$CD68,
          fig_list$CD163, fig_list$FOXP3,fig_list$PD1, fig_list$PDL1,
          fig_list$CD31, fig_list$aSMA, fig_list$Desmin, fig_list$CD45,
          nrow = 4, ncol = 4)




### 5a. SPIAT analysis -------------------------------------------------------
setwd("~/Lin et al - human colorectal cancer/Other data")
df <- readRDS("lin_et_al_3D_spatial_data_separate_slices.rds")

df$Cell.ID <- paste("Cell", seq(nrow(df)), sep = "_")

# minimum_distance_df <- calculate_minimum_distances_between_cell_types_df3D(df,
#                                                                            cell_types_of_interest = NULL,
#                                                                            feature_colname = "Cell.Type1",
#                                                                            show_summary = FALSE,
#                                                                            plot_image = FALSE)


setwd("~/Lin et al - human colorectal cancer/Other data")
minimum_distances_df <- readRDS("minimum_distance_df.rds")

distances <- minimum_distances_df$distance
sorted_distances <- distances[order(distances)]
unique(sorted_distances)[1:200]


### 5b. Calculate minimum distance function using a data frame (rather than a spe object) ----------------
calculate_minimum_distances_between_cell_types_df3D <- function(df,
                                                                cell_types_of_interest = NULL,
                                                                feature_colname = "Cell.Type",
                                                                show_summary = TRUE,
                                                                plot_image = TRUE) {
  
  # If there are no cells, give error
  if (nrow(df) == 0) stop("There are no cells in spe")
  
  # Select all rows in data frame which only contains the cells of interest
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spe object
    unknown_cell_types <- setdiff(cell_types_of_interest, df$Cell.Type)
    if (length(unknown_cell_types) != 0) {
      stop(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                 paste(unknown_cell_types, collapse = ", ")))
    }
    
    df <- df[df[[feature_colname]] %in% cell_types_of_interest, ]
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


# Get absolute minimum distance between cell types (which is not 0)
min(minimum_distance[minimum_distance$distance != 0, "distance"]) # Result: 0.001 microns

setwd("~/Lin et al - human colorectal cancer/Other data")
# write.table(minimum_distance, "minimum_distance_df.csv")

### 6a. FIXING THE DATA -----------------------------------------------------
setwd("~/Lin et al - human colorectal cancer/CRC1_data_updated")
df <- readRDS("CRC1_df.rds")

# Fix df generated from first combining slices
setwd("~/Lin et al - human colorectal cancer/Other_data")
combine_slices_df <- readRDS("lin_et_al_3D_spatial_data_combine_slices.rds")

combine_slices_df <- combine_slices_df[seq(nrow(df)), ]
combine_slices_df[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")] <- df[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]


# Fix df generated from examining separate slices
setwd("~/Lin et al - human colorectal cancer/Other_data")
separate_slices_df <- readRDS("lin_et_al_3D_spatial_data_separate_slices.rds")

separate_slices_df <- separate_slices_df[seq(nrow(df)), ]
separate_slices_df[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")] <- df[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]

# Save
# setwd("~/Lin et al - human colorectal cancer/CRC1_data")
# saveRDS(combine_slices_df, "combine_slices_df.rds")
# saveRDS(separate_slices_df, "separate_slices_df.rds")

### 6b. Plot fixed data -------------------------------------------------------
setwd("~/Lin et al - human colorectal cancer/CRC1_data")

combine_slices_df <- readRDS("combine_slices_df.rds")
separate_slices_df <- readRDS("separate_slices_df.rds")

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

# All cells
plot_cells_rgl_3D(df, 10000, NULL, "Cell.Type.Specific")


# Tumour cells
tumour_cell_types <- c("Tumor/Epi.", "Ki67+ Tumor/Epi.", "PDL1+ Tumor/Epi.")

# Stromal cells
stromal_cell_types <- c("Endothelial", "Muscle/Fibroblast")

# Immune cells
immune_cell_types <- c("Macrophage(I)", "Macrophage(II)", "Macrophage(III)", "Macrophage(IV)", "PDL1+ Macrophage",
                       "PDL1+ lymphocyte",  "DN Lymphocyte", "DP Lymphocyte", "Lymphocyte(III)",
                       "T helper", "PD1+ T helper", "Tc cell", "PD1+ Tc", "Treg",
                       "B cells")


plot_cells_rgl_3D(combine_slices_df, 50000, NULL, "Cell.Type2")
plot_cells_rgl_3D(combine_slices_df, 50000, tumour_cell_types, "Cell.Type2")
plot_cells_rgl_3D(combine_slices_df, 50000, stromal_cell_types, "Cell.Type2")
plot_cells_rgl_3D(combine_slices_df, 50000, immune_cell_types, "Cell.Type2")

plot_cells_rgl_3D(separate_slices_df, 50000, NULL, "Cell.Type2")
plot_cells_rgl_3D(separate_slices_df, 50000, tumour_cell_types, "Cell.Type2")
plot_cells_rgl_3D(separate_slices_df, 50000, stromal_cell_types, "Cell.Type2")
plot_cells_rgl_3D(separate_slices_df, 50000, immune_cell_types, "Cell.Type2")

