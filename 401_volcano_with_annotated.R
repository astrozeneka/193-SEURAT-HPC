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
# Result: genes × samples integer matrix

# ── Join response label from clinical data ────────────────────────────────────
clinical    <- read.csv(
  "../125-WHOLE-SLIDE-ANALYSIS/density_with_cd/Clinical_data.csv",
  stringsAsFactors = FALSE
)
response_map <- setNames(clinical$pCR, clinical$SampleId)

sample_info <- data.frame(
  sample_id = colnames(pseudo_bulk),
  response  = response_map[colnames(pseudo_bulk)],
  row.names = colnames(pseudo_bulk),
  stringsAsFactors = FALSE
)

# Warn about any samples that didn't match clinical data
unmatched <- sample_info$sample_id[is.na(sample_info$response)]
if (length(unmatched) > 0)
  message("Samples with no clinical match (excluded): ", paste(unmatched, collapse = ", "))

# Keep only samples with a known response
keep        <- rownames(sample_info)[!is.na(sample_info$response)]
pseudo_bulk <- pseudo_bulk[, keep, drop = FALSE]
sample_info <- sample_info[keep, , drop = FALSE]

# ── DESeq2 ───────────────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = pseudo_bulk,
  colData   = sample_info,
  design    = ~ response
)
dds$response <- relevel(factor(dds$response), ref = "Non-responder")

dds <- DESeq(dds)
res <- results(dds, contrast = c("response", "Responder", "Non-responder"))

res_df       <- as.data.frame(res)
res_df$gene  <- rownames(res_df)
res_df       <- res_df[order(res_df$padj), ]

write.csv(
  res_df,
  file      = "stats/de_results_pseudobulk_deseq2.csv",
  row.names = FALSE,
  quote     = FALSE
)
