library(Seurat)

larc_seurat <- readRDS("larc_datasets/larc_merged_6k.rds")

rna_counts <- larc_seurat@assays$RNA@counts

# export as MEX format (sparse-friendly, R/Python compatible)
mex_dir <- "rna_mex"
dir.create(mex_dir, recursive = TRUE, showWarnings = FALSE)
Matrix::writeMM(rna_counts, file = paste0(mex_dir, "/matrix.mtx"))
write.csv(rownames(rna_counts), paste0(mex_dir, "/genes.csv"), row.names = FALSE)
write.csv(colnames(rna_counts), paste0(mex_dir, "/barcodes.csv"), row.names = FALSE)