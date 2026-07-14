#!/usr/bin/env python3
"""
07_visualization/plot_flowchart.py

Generates the pipeline overview flowchart (stages: data acquisition -> QC ->
meta-analysis -> genetic correlation -> transcriptome imputation ->
cross-tissue integration -> connectivity mapping -> significance/MOA ->
drug-disease-gene-tissue integration -> candidate compounds).

This is a static, hand-laid-out diagram (not data-driven) -- edit the STAGES
list below to adapt it to a different pipeline.

Usage:
    python plot_flowchart.py --output figures/flowchart.png
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

STAGES = [
    ("GWAS summary statistics acquisition\n(OpenGWAS/GWAS Catalog: 8 traits, VCF format)", "#fff3cd", True),
    ("QC & harmonisation\nVCF\u2192TSV conversion; trait provenance verification;\n"
     "exclusion of case-case (CC-GWAS) contrasts from primary analysis", "#f8d7da", False),
    ("Fixed-effects inverse-variance-weighted meta-analysis\n"
     "(ASD case-control: iPSYCH-PGC 2017 + FinnGen; palindromic SNP exclusion,\n"
     "allele harmonisation, heterogeneity Q/I\u00b2)", "#d1e7dd", False),
    ("LD-aware genetic-correlation screen\n"
     "(distance-pruned Z-score concordance + independent lead-locus test\n"
     "vs. MZ-twin symptom-discordance GWAS)", "#d1e7dd", False),
    ("S-PrediXcan transcriptome imputation\n"
     "(GTEx v8 MASHR models, 49 tissues x 8 traits;\n"
     "varID re-keying, numpy2 compatibility fixes)", "#cfe2ff", True),
    ("Cross-tissue integration & pathway enrichment\n"
     "ACAT p-value combination across tissues; MHC exclusion;\n"
     "Hallmark GSEA per trait", "#e2d9f3", False),
    ("CMap L1000 connectivity mapping\n"
     "(Level5 GCTX \u2192 trt_cp/exemplar filtering \u2192 compound aggregation;\n"
     "5-method reversal scoring: Pearson/Spearman/KS-signed/extreme x2)", "#ffe5d9", True),
    ("Significance calling & MOA enrichment\n"
     "empirical z-score vs. compound library; Fisher's exact MOA test;\n"
     "Broad Repurposing Hub annotation", "#ffe5d9", False),
    ("Drug-Disease-Gene-Tissue integration\n"
     "(leading-gene extraction; GTEx tissue linkage; Sankey visualisation)", "#e2d9f3", False),
]
FINAL = ("Candidate reversal compounds\nfor FXS/ASD drug repurposing", "#d4edda")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    n = len(STAGES)
    box_h, gap = 1.4, 0.8
    total_h = n * (box_h + gap) + box_h + 2  # + final box + top/bottom margin

    fig, ax = plt.subplots(figsize=(13, 17))
    ax.set_xlim(0, 10)
    ax.set_ylim(1.8, total_h - 0.5)
    ax.axis("off")

    def box(x, y, w, h, text, color, fontsize=13, fontweight="normal"):
        b = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08,rounding_size=0.12",
                            facecolor=color, edgecolor="black", linewidth=1.2, zorder=3)
        ax.add_patch(b)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fontsize,
                fontweight=fontweight, zorder=4)

    def arrow(x1, y1, x2, y2):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="-|>", color="black", lw=1.6, shrinkA=2, shrinkB=2), zorder=2)

    y = total_h - 0.5 - box_h
    for text, color, is_bold in STAGES:
        box(0.5, y, 9, box_h, text, color, fontsize=14 if is_bold else 13,
            fontweight="bold" if is_bold else "normal")
        arrow(5, y, 5, y - gap)
        y -= (box_h + gap)

    box(1.5, y - box_h + gap, 7, box_h, FINAL[0], FINAL[1], fontsize=15, fontweight="bold")

    plt.savefig(args.output, dpi=150, bbox_inches="tight", pad_inches=0.15)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
