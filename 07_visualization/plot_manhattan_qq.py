#!/usr/bin/env python3
"""
07_visualization/plot_manhattan_qq.py

Manhattan plot and quantile-quantile (QQ) plot for a meta-analysis output
(as produced by 02_meta_analysis/meta_analysis_ivw.py). Non-suggestive SNPs
are randomly downsampled for plotting speed and output file size; all
suggestive-or-stronger SNPs (p < 1e-3 by default) are always retained.

Usage:
    python plot_manhattan_qq.py \
        --meta_results meta_groupA_ASD.tsv.gz \
        --output figures/manhattan_qq.png \
        --lambda_gc 1.136
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--meta_results", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--lambda_gc", type=float, default=None,
                     help="Genomic control lambda to show in QQ plot title (computed from data if omitted)")
    ap.add_argument("--downsample_background_frac", type=float, default=0.05)
    ap.add_argument("--background_neglogp_threshold", type=float, default=3.0)
    args = ap.parse_args()

    df = pd.read_csv(args.meta_results, sep="\t", dtype={"chromosome": "str"},
                      usecols=["chromosome", "position", "p_meta", "z_meta"])
    df = df[df.chromosome != "X"]
    df["chromosome"] = df.chromosome.astype(int)
    df = df.sort_values(["chromosome", "position"]).reset_index(drop=True)

    chrom_sizes = df.groupby("chromosome").position.max()
    chrom_offsets = chrom_sizes.cumsum().shift(1).fillna(0)
    df["cum_pos"] = df.position + df.chromosome.map(chrom_offsets)
    df["neglogp"] = -np.log10(df.p_meta.clip(lower=1e-300))

    np.random.seed(0)
    keep_sig = df.neglogp >= args.background_neglogp_threshold
    keep_bg = np.random.rand(len(df)) < args.downsample_background_frac
    plot_df = df[keep_sig | keep_bg]
    print(f"Plotting {len(plot_df):,} of {len(df):,} SNPs", flush=True)

    fig, axes = plt.subplots(2, 1, figsize=(16, 10), gridspec_kw={"height_ratios": [2.3, 1]})

    ax = axes[0]
    colors = ["#1f77b4", "#ff7f0e"]
    for i, chrom in enumerate(sorted(df.chromosome.unique())):
        sub = plot_df[plot_df.chromosome == chrom]
        ax.scatter(sub.cum_pos, sub.neglogp, s=4, color=colors[i % 2], rasterized=True)
    ax.axhline(-np.log10(5e-8), color="red", linestyle="--", linewidth=1, label="Genome-wide significance (5x10$^{-8}$)")
    ax.axhline(-np.log10(1e-5), color="gray", linestyle="--", linewidth=0.8, label="Suggestive (1x10$^{-5}$)")
    chrom_mid = df.groupby("chromosome").cum_pos.median()
    ax.set_xticks(chrom_mid); ax.set_xticklabels(chrom_mid.index, fontsize=9)
    ax.set_ylabel("-log10(p)", fontsize=13)
    ax.set_xlabel("Chromosome", fontsize=13)
    ax.set_title("Manhattan plot", fontsize=14)
    ax.legend(fontsize=10, loc="upper right")
    ax.set_ylim(0, max(plot_df.neglogp.max() * 1.05, 8))

    ax2 = axes[1]
    obs = np.sort(df.p_meta.values)
    n = len(obs)
    exp = -np.log10(np.arange(1, n + 1) / (n + 1))
    obs_log = -np.log10(np.clip(obs, 1e-300, 1))
    idx = np.unique(np.round(np.logspace(0, np.log10(max(n - 1, 1)), 20000)).astype(int))
    idx = idx[idx < n]
    ax2.scatter(exp[idx], obs_log[idx], s=3, color="#1f77b4", rasterized=True)
    maxv = max(exp.max(), obs_log.max())
    ax2.plot([0, maxv], [0, maxv], color="red", linewidth=1)

    lam = args.lambda_gc
    if lam is None and "z_meta" in df.columns:
        lam = np.median(df.z_meta.to_numpy() ** 2) / 0.4549
    elif lam is None:
        lam = float("nan")
    ax2.set_xlabel("Expected -log10(p)", fontsize=13)
    ax2.set_ylabel("Observed -log10(p)", fontsize=13)
    ax2.set_title(f"QQ plot (\u03bb={lam:.3f})", fontsize=13)

    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
