# install.packages("ggplot2")
library(ggplot2)
library(Seurat)
library(dplyr)

# Load and merge the two slides
# UpdateSeuratObject() is required to migrate older RDS files to the current
# SeuratObject class definition (e.g. FOV slots added in recent versions)
larc_a_path <- "/larc_datasets/seuratObject_NLG_SS_LARC_A_6K_260226.RDS"
larc_b_path <- "/larc_datasets/seuratObject_NLG_SS_LARC_B_6K_260226.RDS"
larc_a <- UpdateSeuratObject(readRDS(larc_a_path))
larc_b <- UpdateSeuratObject(readRDS(larc_b_path))

# larc_obj <- merge(larc_a, larc_b) # Incorrect code
larc_seurat_meta_data <- rbind(larc_a@meta.data, larc_b@meta.data)
# print column names of larc_seurat_meta_data
print(colnames(larc_seurat_meta_data))
# Export to CSV
write.csv(larc_seurat_meta_data, "larc_datasets/larc_merged_meta_data.csv", row.names = TRUE)

