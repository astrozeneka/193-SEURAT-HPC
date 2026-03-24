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
#                 ┌───────────────┴───────────────────┐
#                 ▼                                   ▼
#          [DESeq2 branch]                    [CMS branch] ← NEW 501
#       lines 51-95 (existing)         raw counts → TPM → log2 → Entrez
#                                      → getFunctionalSpectra → DeepCC label

# ══════════════════════════════════════════════════════════════════════════════
# CMS CLASSIFICATION BRANCH  (DeepCC functional spectra)
# ══════════════════════════════════════════════════════════════════════════════
library(biomaRt)
library(org.Hs.eg.db)
library(DeepCC)
library(reticulate)
use_condaenv("r-tf23", required = TRUE)
library(tensorflow)

# Load rds
pseudo_bulk <- readRDS(file = "rdata/pseudo_bulk_counts.rds")

# ── Step 1. Gene lengths via biomaRt (symbol → median transcript length) ─────
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")
tx_lengths <- getBM(
  attributes = c("hgnc_symbol", "transcript_length"),
  filters    = "hgnc_symbol",
  values     = rownames(pseudo_bulk),
  mart       = mart
)
gene_lengths <- aggregate(transcript_length ~ hgnc_symbol, data = tx_lengths, FUN = median)

# ── Step 2. Raw counts → TPM ──────────────────────────────────────────────────
counts_to_tpm <- function(counts, lengths) {
  rpk <- counts / (lengths / 1000)
  t(t(rpk) / colSums(rpk) * 1e6)
}

common_genes <- intersect(rownames(pseudo_bulk), gene_lengths$hgnc_symbol)
lengths_vec  <- gene_lengths$transcript_length[match(common_genes, gene_lengths$hgnc_symbol)]
tpm_mat      <- counts_to_tpm(as.matrix(pseudo_bulk[common_genes, ]), lengths_vec)

# ── Step 3. log2(TPM + 1)  — no ComBat: single cohort, no batch effect ────────
log2tpm <- log2(tpm_mat + 1)

# ── Step 4. Map gene symbols → Entrez IDs, remove duplicates, transpose ───────
entrez_map <- AnnotationDbi::select(org.Hs.eg.db,
  keys    = rownames(log2tpm),
  columns = "ENTREZID",
  keytype = "SYMBOL"
)
entrez_map <- entrez_map[!is.na(entrez_map$ENTREZID) & !duplicated(entrez_map$SYMBOL), ]

log2tpm_entrez <- log2tpm[entrez_map$SYMBOL, ]
rownames(log2tpm_entrez) <- entrez_map$ENTREZID
log2tpm_entrez <- log2tpm_entrez[!duplicated(rownames(log2tpm_entrez)), ]

# Transpose: samples as rows, Entrez IDs as columns (DeepCC format)
expr <- t(log2tpm_entrez) |> as.data.frame()

# ── Step 5. Functional spectra ────────────────────────────────────────────────
fs <- getFunctionalSpectra(expr, geneSets = "MSigDBv7")
saveRDS(fs, file = "rdata/FS_larc.rds")

# ── Step 6. CMS labels ────────────────────────────────────────────────────────
cms_model <- load_DeepCC_model("../174-CMS-CLASSIFIER/deepcc_model/CRC_TCGA")
labels    <- get_DeepCC_label(cms_model, fs, cutoff = 0.5)

write.csv(
  data.frame(sample = rownames(fs), CMS = labels),
  file      = "cms/CMS_labels_larc.csv",
  row.names = FALSE
)
