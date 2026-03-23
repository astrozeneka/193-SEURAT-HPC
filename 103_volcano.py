import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np

slug = "all"
df = pd.read_csv(f"stats/de_results_responder_vs_nonresponder_{slug}.csv")

if __name__ == '__main__':

    sig = df[df["significant"] == True]
    ns  = df[df["significant"] == False]
    up_r  = sig[sig["direction"] == "Up in Responder"]
    up_nr = sig[sig["direction"] == "Up in Non-responder"]

    fig, ax = plt.subplots(figsize=(9, 7))

    ax.scatter(ns["log2fc"],  -np.log10(ns["pval"]),  color="#aaaaaa", s=10, alpha=0.5, linewidths=0, label="NS")
    ax.scatter(up_r["log2fc"],  -np.log10(up_r["pval"]),  color="#e63946", s=14, alpha=0.8, linewidths=0, label="Up in Responder")
    ax.scatter(up_nr["log2fc"], -np.log10(up_nr["pval"]), color="#457b9d", s=14, alpha=0.8, linewidths=0, label="Up in Non-responder")

    # Label top significant genes by -log10(pval)
    top = sig.nlargest(20, "neg_log10_padj")
    for _, row in top.iterrows():
        color = "#e63946" if row["direction"] == "Up in Responder" else "#457b9d"
        ax.text(
            row["log2fc"], -np.log10(row["pval"]),
            row["gene"], fontsize=6.5, color=color, ha="center", va="bottom",
            path_effects=[pe.withStroke(linewidth=2, foreground="white")],
        )

    ax.axvline(0, color="black", lw=0.6, ls="--", alpha=0.4)
    ax.set_xlabel("Log2 Fold Change", fontsize=12)
    ax.set_ylabel("-log10(p-value)", fontsize=12)
    ax.set_title("Responder vs Non-responder — Volcano Plot", fontsize=13, fontweight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    plt.savefig(f"stats/103_volcano_{slug}.png", dpi=150)
