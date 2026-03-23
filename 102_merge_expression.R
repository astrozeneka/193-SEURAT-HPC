library(ggplot2)
library(Seurat)
library(dplyr)
library(Matrix)
library(future.apply)

# Load and merge the two slides
# UpdateSeuratObject() is required to migrate older RDS files to the current
# SeuratObject class definition (e.g. FOV slots added in recent versions)
larc_a_path <- "/app/larc_datasets/seuratObject_NLG_SS_LARC_A_6K_260226.RDS"
larc_b_path <- "/app/larc_datasets/seuratObject_NLG_SS_LARC_B_6K_260226.RDS"
larc_a <- UpdateSeuratObject(readRDS(larc_a_path))
larc_b <- UpdateSeuratObject(readRDS(larc_b_path))


# ── Merge ──────────────────────────────────────────────────────────────────
larc_obj <- merge(
  x = larc_a,
  y = larc_b,
  add.cell.ids = NULL,   # NULL because your IDs are already unique
  merge.data  = TRUE,    # carries over any existing normalized layers too
  project     = "LARC_6K"
)
# save rds
saveRDS(larc_obj, file = "larc_datasets/larc_merged_6k.rds")
# Load using
# larc_obj <- readRDS("larc_datasets/larc_merged_6k.rds")

# Merge QC
# Total cell count should equal sum of both objects
ncol(larc_a) + ncol(larc_b) == ncol(larc_obj)   # should be TRUE

# Confirm gene set is identical (should be 6000 for both)
nrow(larc_a) == nrow(larc_b)   # TRUE → safe to merge
nrow(larc_obj)                 # should still be 6000


# ── Tag responder / non-responder ──────────────────────────────────────────
# Cell lists were built by 101_cell_ids_per_treatment.py from clinical pCR data
responder_cells     <- readLines("cell_list/Responder.txt")
non_responder_cells <- readLines("cell_list/Non-responder.txt")

response_vec <- setNames(
  rep(NA_character_, ncol(larc_obj)),
  Cells(larc_obj)
)
response_vec[intersect(Cells(larc_obj), responder_cells)]     <- "Responder"
response_vec[intersect(Cells(larc_obj), non_responder_cells)] <- "Non-responder"

larc_obj <- AddMetaData(larc_obj, metadata = response_vec, col.name = "response")

# Sanity check
table(larc_obj$response, useNA = "ifany")


# -── Normalize ───────────────────────────────────────────────────────────────
larc_obj <- NormalizeData( # Need Server to Compute
  larc_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  assay = "RNA"   # CosMX objects use "Nanostring", not "RNA"
)
# remove cells with unassigned
larc_obj_de <- subset(larc_obj, subset = response %in% c("Responder", "Non-responder"))


# -─ Extract Matrix, Compute Log2FC + Wilcoxon p-value ──────────────────────────
# Extract sparse normalized matrix (genes × cells)
norm_matrix <- GetAssayData(larc_obj_de, assay = "RNA", layer = "data")

# Get cell indices per group
responder_cells    <- WhichCells(larc_obj_de, expression = response == "Responder")
non_responder_cells <- WhichCells(larc_obj_de, expression = response == "Non-responder")

# Subset matrix by group (still sparse)
mat_resp    <- norm_matrix[, responder_cells]
mat_nonresp <- norm_matrix[, non_responder_cells]

# Row means per group (dense vector, one value per gene)
mean_resp    <- rowMeans(mat_resp)
mean_nonresp <- rowMeans(mat_nonresp)

# Log2FC: Responder vs Non-responder
# Already log-normalized so subtraction = log2FC (approx)
log2fc <- mean_resp - mean_nonresp


# --- Wilcoxon rank-sum test per gene (parallelized) ---
# Convert to dense matrix once — avoids repeated sparse slicing overhead
mat_resp_dense    <- as.matrix(mat_resp)
mat_nonresp_dense <- as.matrix(mat_nonresp)
genes <- rownames(norm_matrix)

pval_results <- sapply(genes, function(gene) {
  x <- mat_resp_dense[gene, ]
  y <- mat_nonresp_dense[gene, ]
  wt <- wilcox.test(x, y, exact = FALSE)
  return(wt$p.value)
})

padj <- p.adjust(pval_results, method = "BH")

de_results <- data.frame(
  gene    = genes,
  log2fc  = log2fc,
  pval    = pval_results,
  padj    = padj
)

# Add significance label
de_results$significant <- de_results$padj < 0.05 & abs(de_results$log2fc) > 0.25

head(de_results[order(de_results$padj), ])

# Optional: clean and sort before exporting
de_results_export <- de_results %>%
  arrange(padj) %>%
  mutate(
    neg_log10_padj = -log10(padj),          # Pre-compute for volcano plot
    direction = case_when(
      log2fc >  0.25 & padj < 0.05 ~ "Up in Responder",
      log2fc < -0.25 & padj < 0.05 ~ "Up in Non-responder",
      TRUE                          ~ "NS"  # Not significant
    )
  )

# Export
write.csv(
  de_results_export,
  file = "stats/de_results_responder_vs_nonresponder_all.csv",
  row.names = FALSE,
  quote     = FALSE
)