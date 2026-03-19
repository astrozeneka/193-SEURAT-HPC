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


# ── Merge ──────────────────────────────────────────────────────────────────
larc_obj <- merge(
  x = larc_a,
  y = larc_b,
  add.cell.ids = NULL,   # NULL because your IDs are already unique
  merge.data  = TRUE,    # carries over any existing normalized layers too
  project     = "LARC_6K"
)

# Merge QC
# Total cell count should equal sum of both objects
ncol(larc_a) + ncol(larc_b) == ncol(larc_obj)   # should be TRUE

# Confirm gene set is identical (should be 6000 for both)
nrow(larc_a) == nrow(larc_b)   # TRUE → safe to merge
nrow(larc_obj)                 # should still be 6000