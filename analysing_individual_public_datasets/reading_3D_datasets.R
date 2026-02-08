# Code to read and pre-process in each 3D public dataset.
# CyCIF colorectal cancer ------
setwd("~/R/data3D/CyCIF_colorectal_cancer")
data3D <- read.csv("colorectal_cancer_df.csv")
data3D <- data3D[, -1]
colnames(data3D)[1:3] <- c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")
colnames(data3D)[4:5] <- c("Cell.Type.Specific", "Cell.Type")
data3D$Cell.Type[data3D$Cell.Type == "Tumor/Epi"] <- "Tumour"

# pick and choose the cell types...
cell_types <- c(
  "Tumour",
  "Immune"
)


# Merfish cortex ------
setwd("~/R/data3D/merfish_mouse_brain")
data3D <- read.csv("mouse_cortex_100um_df.csv")
colnames(data3D) <- c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position", "Cell.Type")
# Bin Cell.Z.Position so that any z-coord between 0 and 9 is labelled 0, 10 and 19 is labelled 10 and so on.
data3D$Cell.Z.Position <- (data3D$Cell.Z.Position %/% 10) * 10

data3D$Cell.Type[data3D$Cell.Type == "Astro"] <- "ASC"
data3D$Cell.Type[data3D$Cell.Type == "Oligo"] <- "OGC"
data3D$Cell.Type[data3D$Cell.Type == "Micro"] <- "MGC"

cell_types <- c("INC", "EXC", "OGC", "OPC", "Endo", "MGC", "ASC")
# Merfish hypothalamus -----
setwd("~/R/data3D/merfish_mouse_brain")
data3D <- read.csv("mouse_hypothalamus_200um_df.csv")
colnames(data3D) <- c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position", "Cell.Type")
# Bin Cell.Z.Position so that any z-coord between 0 and 9 is labelled 0, 10 and 19 is labelled 10 and so on.
data3D$Cell.Z.Position <- (data3D$Cell.Z.Position %/% 10) * 10

cell_types <- c("INC", "EXC", "OGC", "OPC", "Endo", "MGC", "ASC")

# openST human metastatic lymph node -----
setwd("~/R/data3D/openST_human_metastatic_lymph_node")
data3D <- read.csv("human_metastatic_lymph_node_df.csv")

cell_types <- unique(data3D$Cell.Type)
cell_types <- cell_types[cell_types != "unknown"]

# spateo mouse embryo -----
# spateo
setwd("~/R/data3D/spateo")
data3D <- read.csv("mouse_E11.5_embryo.csv")
data3D$Cell.Type[data3D$Cell.Type == ""] <- "Empty"

# Pick 10 most abundant cell types
cell_types <- names(sort(table(data3D$Cell.Type), decreasing = TRUE))[1:10]

# STOMICS fly datasets -------
setwd("~/R/data3D/stomics_fly")
data3D <- read.csv("E14-16h_df.csv")
cell_types <- c("CNS", "epidermis", "carcass", "muscle", "fat body")

setwd("~/R/data3D/stomics_fly")
data3D <- read.csv("E16-18h_df.csv")
cell_types <- c("CNS", "epidermis", "carcass", "muscle", "fat body")

setwd("~/R/data3D/stomics_fly")
data3D <- read.csv("L1_df.csv")
cell_types <- c("CNS", "epidermis", "carcass", "muscle", "fat body")

setwd("~/R/data3D/stomics_fly")
data3D <- read.csv("L2_df.csv")
cell_types <- c("CNS", "epidermis", "carcass", "muscle", "fat body")

setwd("~/R/data3D/stomics_fly")
data3D <- read.csv("L3_df.csv")
cell_types <- c("CNS", "epidermis", "carcass", "muscle", "fat body")


