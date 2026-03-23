import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np

df = pd.read_csv("stats/de_results_pseudobulk_deseq2.csv")

PADJ_THRESH = 0.05
LFC_THRESH  = 1.0

df["neg_log10_padj"] = -np.log10(df["padj"].clip(lower=1e-300))

df["significant"] = (df["padj"] < PADJ_THRESH) & (df["log2FoldChange"].abs() > LFC_THRESH)
df["direction"]   = np.where(df["log2FoldChange"] > 0, "Up", "Down")

if __name__ == '__main__':

    sig  = df[df["significant"] == True]
    ns   = df[df["significant"] == False]
    up   = sig[sig["direction"] == "Up"]
    down = sig[sig["direction"] == "Down"]

    fig, ax = plt.subplots(figsize=(9, 7))

    ax.scatter(ns["log2FoldChange"],   -np.log10(ns["pvalue"]),   color="#aaaaaa", s=10, alpha=0.5, linewidths=0, label="NS")
    ax.scatter(up["log2FoldChange"],   -np.log10(up["pvalue"]),   color="#e63946", s=14, alpha=0.8, linewidths=0, label="Up")
    ax.scatter(down["log2FoldChange"], -np.log10(down["pvalue"]), color="#457b9d", s=14, alpha=0.8, linewidths=0, label="Down")

    # Label top significant genes by -log10(padj)
    top = sig.nlargest(20, "neg_log10_padj")
    for _, row in top.iterrows():
        color = "#e63946" if row["direction"] == "Up" else "#457b9d"
        ax.text(
            row["log2FoldChange"], -np.log10(row["pvalue"]),
            row["gene"], fontsize=6.5, color=color, ha="center", va="bottom",
            path_effects=[pe.withStroke(linewidth=2, foreground="white")],
        )

    ax.axvline(0,            color="black", lw=0.6, ls="--", alpha=0.4)
    ax.axvline( LFC_THRESH,  color="grey",  lw=0.6, ls=":",  alpha=0.5)
    ax.axvline(-LFC_THRESH,  color="grey",  lw=0.6, ls=":",  alpha=0.5)
    ax.axhline(-np.log10(PADJ_THRESH), color="grey", lw=0.6, ls=":", alpha=0.5)

    ax.set_xlabel("Log2 Fold Change", fontsize=12)
    ax.set_ylabel("-log10(p-value)", fontsize=12)
    ax.set_title("Pseudobulk DESeq2 — Volcano Plot", fontsize=13, fontweight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    plt.savefig("stats/303_volcano_pseudobulk.png", dpi=150)