setwd("C:/Users/Me/OneDrive - The University of Melbourne/PeterMac/Honours 2024/Publicly available 3D spatial data/Rheumatoid Arthritis")

## Read textfiles
datasets <- list()
for (i in 1:5) {
  file_path = paste("RA4_", i, "_annotations.txt" ,sep = "")
  datasets[[i]] <- read.delim(file = file_path,
                              header = TRUE) 
}


## Combine x, y, z coords and cell_type info into one data frame
x <- c()
y <- c()
z <- c()
cell_type <- c()
for (i in 1:5) {
  data <- datasets[[i]]
  x <- append(x, data$x)
  y <- append(y, data$y)
  z <- append(z, rep(i, nrow(data))) # z coord is the current i value
  cell_type <- append(cell_type, data$value)
}

df <- data.frame(x = x,
                 y = y,
                 z = z,
                 cell_type = cell_type)


## Plot in 3D
library(rgl)
library("RColorBrewer")

colors <- brewer.pal(n = 11, name = 'Spectral')

plot3d(x = df$x,
       y = df$y,
       z = df$z,
       col = get_colors(df$cell_type, colors),
       size = 4)

aspect3d(5, 5, 1)

cells <- as.factor(df$cell_type)
legend3d("topright", legend = levels(cells) , pch = 19, col = colors[seq(length(cells))], inset = c(0.02))




get_colors <- function(groups, group.col = palette()) {
  groups <- as.factor(groups)
  ngrps <- length(levels(groups))
  
  if (ngrps > length(group.col)) { 
    group.col <- rep(group.col, ngrps)
  }
  color <- group.col[as.numeric(groups)]
  names(color) <- as.vector(groups)
  
  return(color)
}
