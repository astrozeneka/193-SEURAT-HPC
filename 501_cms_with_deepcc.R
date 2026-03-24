library(ggplot2)
library(Seurat)
library(dplyr)
library(Matrix)
library(future.apply)

larc_obj <- readRDS("larc_datasets/larc_merged_6k.rds")

# ── Filter to annotated cells & tag source sample ─────────────────────────────
annotation_files <- list.files(
  "../191-SPLIT-COSMIX-SAMPLES/splitted/roi_annotated",
  pattern = "\\.csv$", full.names = TRUE
)

cell_sample_list <- lapply(annotation_files, function(f) {
  sname <- basename(f)
  sname <- sub("^LARC_[AB]_", "", sname)
  sname <- sub("_annotated\\.csv$", "", sname)
  df    <- read.csv(f, row.names = 1)
  data.frame(
    cell      = rownames(df)[df$selection_mask == 1],
    sample_id = sname,
    stringsAsFactors = FALSE
  )
})

cell_sample_df <- do.call(rbind, cell_sample_list)
larc_obj <- subset(larc_obj, cells = cell_sample_df$cell)

# Add sample_id to each cell's metadata
sample_vec <- setNames(cell_sample_df$sample_id, cell_sample_df$cell)
larc_obj   <- AddMetaData(larc_obj, metadata = sample_vec[Cells(larc_obj)], col.name = "sample_id")

# ══════════════════════════════════════════════════════════════════════════════
# PSEUDO-BULK + DESeq2
# ══════════════════════════════════════════════════════════════════════════════
library(DESeq2)

# ── Aggregate raw counts per gene per sample ──────────────────────────────────
# Uses the counts slot (integer, un-normalised) — DESeq2 requires raw counts
raw_counts <- GetAssayData(larc_obj, assay = "RNA", layer = "counts")
bulk_meta  <- larc_obj@meta.data

pseudo_bulk <- sapply(unique(bulk_meta$sample_id), function(s) {
  cells <- rownames(bulk_meta)[bulk_meta$sample_id == s]
  rowSums(raw_counts[, cells, drop = FALSE])
})

# Save RDS for DESeq2 branch (for alter use)
saveRDS(pseudo_bulk, file = "rdata/pseudo_bulk_counts.rds")

#   lines 7-48  → load, filter annotated cells, aggregate → pseudo_bulk
#                                 │
#                 ┌───────────────┴─────────────────┐
#                 ▼                                   ▼
#          [DESeq2 branch]                    [CMS branch] ← NEW 501
#       lines 51-95 (existing)         raw counts → TPM → log2 → Entrez
#                                      → getFunctionalSpectra → DeepCC label

