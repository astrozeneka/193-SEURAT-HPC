library(Seurat)
library(decoupleR)
# install with BiocManager::install("decoupleR")

larc_obj <- readRDS("larc_datasets/larc_merged_6k.rds")

# ── Filter to annotated cells ────────────────────────────────────────────────
annotated_cells <- do.call(c, lapply(
  list.files("../191-SPLIT-COSMIX-SAMPLES/splitted/roi_annotated", pattern = "\\.csv$", full.names = TRUE),
  function(f) { df <- read.csv(f, row.names = 1); rownames(df)[df$selection_mask == 1] }
))
larc_obj <- subset(larc_obj, cells = annotated_cells)

larc_obj <- NormalizeData(
  larc_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  assay = "RNA"   # CosMX objects use "Nanostring", not "RNA"
)

# ── Subsample for quick test ────────────────────────────────────────────────
set.seed(42)
n_cells <- 1000
# sampled_cells <- sample(Cells(larc_obj), n_cells)
#larc_sub <- subset(larc_obj, cells = sampled_cells)
larc_sub <- larc_obj

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

# Export as CSV
write.csv(tf_acts, "tf_activities_ulm.csv", row.names = FALSE)