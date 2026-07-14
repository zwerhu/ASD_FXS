#!/usr/bin/env python3
"""
07_visualization/plot_moa_enrichment.py

Bar plot of the top N mechanism-of-action (MOA) terms by nominal Fisher's
exact test p-value, from 06_drug_repurposing/03_significance_and_moa_enrichment.py.

Usage:
    python plot_moa_enrichment.py \
        --moa_enrichment cmap_moa_enrichment.tsv \
        --output figures/moa_enrichment.png --top_n 15
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--moa_enrichment", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--top_n", type=int, default=15)
    args = ap.parse_args()

    moa = pd.read_csv(args.moa_enrichment, sep="\t")
    top = moa.sort_values("pvalue").head(args.top_n).iloc[::-1]

    fig, ax = plt.subplots(figsize=(11, 8))
    ax.barh(top.MOA, -np.log10(top.pvalue), color="#e07a5f", edgecolor="black")
    for i, (_, row) in enumerate(top.iterrows()):
        ax.text(-np.log10(row.pvalue) + 0.02, i, f"{int(row.sig_count)}/{int(row.total_count)}", va="center", fontsize=10)
    ax.axvline(-np.log10(0.05), color="gray", linestyle="--", linewidth=1, label="nominal p=0.05")
    ax.set_xlabel("-log10(nominal p-value)  [Fisher's exact test]", fontsize=13)
    any_fdr_sig = (moa.fdr < 0.25).any() if "fdr" in moa.columns else False
    subtitle = "some terms reach FDR<0.25" if any_fdr_sig else "none survive FDR correction \u2014 exploratory leads only"
    ax.set_title(f"Top {args.top_n} MOA enrichment among significant reversal candidates\n"
                 f"({subtitle}, numbers = sig/total)", fontsize=13)
    ax.tick_params(axis='y', labelsize=11)
    ax.tick_params(axis='x', labelsize=11)
    ax.legend(fontsize=11)
    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
