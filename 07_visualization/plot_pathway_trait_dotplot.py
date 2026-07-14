#!/usr/bin/env python3
"""
07_visualization/plot_pathway_trait_dotplot.py

Dot plot of Hallmark GSEA normalised enrichment scores (NES) for pathways
reaching FDR < threshold in at least one trait, shown across all traits for
comparison; dot size encodes nominal significance, colour encodes NES
direction, and a coloured outline flags trait-specific FDR significance.

Usage:
    python plot_pathway_trait_dotplot.py \
        --gsea_results gsea_ALL_hallmark_results.tsv.gz \
        --output figures/pathway_trait_dotplot.png --fdr_threshold 0.25
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gsea_results", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--fdr_threshold", type=float, default=0.25)
    args = ap.parse_args()

    gsea = pd.read_csv(args.gsea_results, sep="\t")
    gsea["NOM p-val"] = gsea["NOM p-val"].astype(float)
    gsea["FDR q-val"] = gsea["FDR q-val"].astype(float)

    sig_terms = gsea.loc[gsea["FDR q-val"] < args.fdr_threshold, "Term"].unique()
    sub = gsea[gsea.Term.isin(sig_terms)].copy()
    sub["neglog10p"] = -np.log10(sub["NOM p-val"].clip(lower=1e-4))
    sub["Term_short"] = sub.Term.str.replace("HALLMARK_", "", regex=False)

    traits = sorted(sub.trait.unique())
    terms = sorted(sub.Term_short.unique())
    sub["trait_i"] = sub.trait.apply(lambda t: traits.index(t))
    sub["term_i"] = sub.Term_short.apply(lambda t: terms.index(t))

    fig, ax = plt.subplots(figsize=(1.4 * len(traits) + 3, 0.6 * len(terms) + 2))
    vmax = np.nanmax(np.abs(sub.NES))
    sc = ax.scatter(sub.trait_i, sub.term_i, s=sub.neglog10p * 80 + 30, c=sub.NES,
                     cmap="RdBu_r", vmin=-vmax, vmax=vmax, edgecolors="black", linewidths=0.7)
    fdr_sig = sub[sub["FDR q-val"] < args.fdr_threshold]
    ax.scatter(fdr_sig.trait_i, fdr_sig.term_i, s=fdr_sig.neglog10p * 80 + 30, facecolors="none",
               edgecolors="lime", linewidths=2.2)

    ax.set_xticks(range(len(traits))); ax.set_xticklabels(traits, rotation=35, ha="right")
    ax.set_yticks(range(len(terms))); ax.set_yticklabels(terms)
    ax.set_xlim(-0.5, len(traits) - 0.5); ax.set_ylim(-0.5, len(terms) - 0.5)
    plt.colorbar(sc, ax=ax, label="NES (Normalized Enrichment Score)")
    ax.set_title(f"GSEA Hallmark pathways x trait\n(dot size = -log10(nominal p); "
                 f"green outline = FDR<{args.fdr_threshold})", fontsize=11)
    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
