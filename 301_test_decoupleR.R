library(Seurat)
library(decoupleR)
# install with BiocManager::install("decoupleR")

larc_obj <- readRDS("larc_datasets/larc_merged_6k.rds")

larc_obj <- NormalizeData(
  larc_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  assay = "RNA"   # CosMX objects use "Nanostring", not "RNA"
)

# ── Subsample for quick test ────────────────────────────────────────────────
set.seed(42)
n_cells <- 500
sampled_cells <- sample(Cells(larc_obj), n_cells)
larc_sub <- subset(larc_obj, cells = sampled_cells)

# Extract log-normalized matrix (genes × cells)
log_mat <- GetAssayData(larc_sub, assay = "RNA", layer = "data")

# ── Load CollecTRI prior network ────────────────────────────────────────────
net <- get_collectri(organism = "human", split_complexes = FALSE)
head(net)

# ── Run TF activity inference (ULM) ────────────────────────────────────────
# minsize: minimum number of targets a TF must have in the matrix to be tested
tf_acts <- run_ulm(
  mat     = as.matrix(log_mat),
  net     = net,
  .source = "source",
  .target = "target",
  .mor    = "mor",
  minsize = 5
)

head(tf_acts)