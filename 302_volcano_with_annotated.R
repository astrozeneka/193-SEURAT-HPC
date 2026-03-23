library(ggplot2)
library(Seurat)
library(dplyr)
library(Matrix)
library(future.apply)

larc_obj <- readRDS("larc_datasets/larc_merged_6k.rds")

# ── Filter to annotated cells ────────────────────────────────────────────────
annotated_cells <- do.call(c, lapply(
  list.files("../191-SPLIT-COSMIX-SAMPLES/splitted/roi_annotated", pattern = "\\.csv$", full.names = TRUE),
  function(f) { df <- read.csv(f, row.names = 1); rownames(df)[df$selection_mask == 1] }
))
larc_obj <- subset(larc_obj, cells = annotated_cells)
