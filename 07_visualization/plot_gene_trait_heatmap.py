#!/usr/bin/env python3
"""
07_visualization/plot_gene_trait_heatmap.py

Heatmap of signed cross-tissue significance (-log10 ACAT p-value x sign of
mean Z-score) for the top N integrated S-PrediXcan genes per trait, across
all traits.

Usage:
    python plot_gene_trait_heatmap.py \
        --integrative_results integrative_ACAT_results.tsv.gz \
        --output figures/gene_trait_heatmap.png --top_n 5
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--integrative_results", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--top_n", type=int, default=5, help="Top N genes per trait to include")
    args = ap.parse_args()

    integ = pd.read_csv(args.integrative_results, sep="\t")

    top_per_trait = integ.sort_values("p_acat").groupby("trait").head(args.top_n)
    top_genes = top_per_trait.gene_name.unique().tolist()

    sub = integ[integ.gene_name.isin(top_genes)].copy()
    sub["signed_logp"] = -np.log10(sub.p_acat) * np.sign(sub.mean_zscore)

    mat = sub.pivot_table(index="gene_name", columns="trait", values="signed_logp")
    mat = mat.loc[top_genes]

    fig, ax = plt.subplots(figsize=(10, max(4, 0.35 * len(mat))))
    vmax = np.nanmax(np.abs(mat.values))
    sns.heatmap(mat, cmap="RdBu_r", center=0, vmin=-vmax, vmax=vmax,
                annot=True, fmt=".1f", cbar_kws={"label": "signed -log10(p_ACAT)"},
                linewidths=0.5, linecolor="white", ax=ax)
    ax.set_title(f"Top integrated gene x trait associations\n"
                 f"(union of top {args.top_n} genes per trait, ACAT-combined across tissues)", fontsize=11)
    ax.set_xlabel("Trait (GWAS)")
    ax.set_ylabel("Gene")
    plt.xticks(rotation=35, ha="right")
    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
