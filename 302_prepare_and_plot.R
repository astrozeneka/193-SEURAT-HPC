library(tidyverse)
library(ggplot2)
library(ggrepel)

# ── 1. Load data ──────────────────────────────────────────────────────────────

ulm_res <- read_csv("your_ulm_results.csv")

# Your CSV columns: statistic, source, condition, score, p_value
# 'condition' = cell barcode

# ── 2. Build a metadata table ─────────────────────────────────────────────────
# You need to map each cell barcode → patient_id + response_group
# Pull this from your Seurat object metadata separately and save, OR
# reconstruct here if barcodes encode sample info

# Example: if barcodes follow a pattern like "c_1_107_cellnum"
# you may need to load Seurat metadata separately:
#
#   meta <- read_csv("seurat_metadata.csv") %>%
#     select(barcode, patient_id, response_group)
#
# For now, assuming you have it as a separate file:

meta <- read_csv("seurat_metadata.csv") %>%
  rename(condition = barcode) %>%          # match the 'condition' column name
  select(condition, patient_id, response_group)

# ── 3. Join metadata to ULM results ──────────────────────────────────────────

ulm_annotated <- ulm_res %>%
  inner_join(meta, by = "condition")

# Quick sanity check
cat("Cells matched:", n_distinct(ulm_annotated$condition), "\n")
cat("TFs found:",     n_distinct(ulm_annotated$source),    "\n")
cat("Patients found:", n_distinct(ulm_annotated$patient_id), "\n")

# ── 4. Pseudo-bulk: mean TF activity per patient per TF ───────────────────────

pseudobulk <- ulm_annotated %>%
  group_by(source, patient_id, response_group) %>%
  summarise(
    mean_score   = mean(score,   na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    n_cells      = n(),
    .groups = "drop"
  )

# ── 5. Option 1 — Per-TF activity distribution plot ──────────────────────────
# Shows patient-level mean scores per TF, split by response group
# Useful for a focused subset of TFs (pick top N by variance or known biology)

# 5a. Identify top TFs by variance across patients (data-driven selection)
top_tfs <- pseudobulk %>%
  group_by(source) %>%
  summarise(score_var = var(mean_score), .groups = "drop") %>%
  slice_max(score_var, n = 20) %>%
  pull(source)

# OR: manually specify TFs of interest based on your IF biology
# top_tfs <- c("STAT3", "IRF4", "FOXP3", "MYC", "HIF1A", "NFKB1", "TGFB1")

# 5b. Plot
pseudobulk %>%
  filter(source %in% top_tfs) %>%
  mutate(
    source         = factor(source, levels = top_tfs),  # preserve ordering
    response_group = factor(response_group, levels = c("responder", "non_responder"))
  ) %>%
  ggplot(aes(x = response_group, y = mean_score, color = response_group)) +

  # Individual patient points
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.8) +

  # Median bar
  stat_summary(fun = median, geom = "crossbar",
               width = 0.4, linewidth = 0.6, color = "black") +

  scale_color_manual(values = c("responder"     = "#E05C5C",
                                "non_responder" = "#5C7AE0")) +
  facet_wrap(~ source, scales = "free_y", ncol = 5) +
  labs(
    x     = NULL,
    y     = "Mean TF Activity Score (ULM)",
    title = "TF Activity Distribution by Response Group",
    color = "Group"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 9),
    axis.text.x      = element_text(angle = 30, hjust = 1),
    legend.position  = "bottom"
  )

ggsave("tf_activity_distribution.pdf", width = 14, height = 10)

# ── 6. Option 2 — Volcano plot (group-level comparison) ───────────────────────

tf_stats <- pseudobulk %>%
  group_by(source) %>%
  summarise(
    mean_R   = mean(mean_score[response_group == "responder"],     na.rm = TRUE),
    mean_NR  = mean(mean_score[response_group == "non_responder"], na.rm = TRUE),
    score_diff = mean_R - mean_NR,
    p_value  = tryCatch(
      wilcox.test(
        mean_score[response_group == "responder"],
        mean_score[response_group == "non_responder"],
        exact = FALSE
      )$p.value,
      error = function(e) NA_real_
    ),
    n_R  = sum(response_group == "responder"),
    n_NR = sum(response_group == "non_responder"),
    .groups = "drop"
  ) %>%
  filter(!is.na(p_value)) %>%
  mutate(
    p_adj          = p.adjust(p_value, method = "BH"),
    neg_log10_padj = -log10(p_adj),
    significance   = case_when(
      p_adj < 0.05 & score_diff >  0.5 ~ "Higher in Responders",
      p_adj < 0.05 & score_diff < -0.5 ~ "Higher in Non-responders",
      TRUE                              ~ "NS"
    ),
    label = ifelse(significance != "NS", source, NA)
  )

# Volcano plot
ggplot(tf_stats, aes(x = score_diff, y = neg_log10_padj)) +
  geom_point(aes(color = significance), size = 2.5, alpha = 0.75) +
  scale_color_manual(values = c(
    "Higher in Responders"     = "#E05C5C",
    "Higher in Non-responders" = "#5C7AE0",
    "NS"                       = "grey70"
  )) +
  geom_text_repel(
    aes(label = label), size = 3,
    max.overlaps = 20, box.padding = 0.4
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-0.5, 0.5),  linetype = "dashed", color = "grey40") +
  labs(
    x     = "Mean TF Activity Difference (Responders − Non-responders)",
    y     = "-log10(BH-adjusted p-value)",
    title = "TF Activity: Responders vs. Non-responders",
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")

ggsave("tf_activity_volcano.pdf", width = 8, height = 7)

# ── 7. Export the stats table ─────────────────────────────────────────────────

write_csv(tf_stats, "tf_group_comparison_stats.csv")