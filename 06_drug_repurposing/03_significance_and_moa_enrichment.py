#!/usr/bin/env python3
"""
06_drug_repurposing/03_significance_and_moa_enrichment.py

Two steps:

1. Significance calling: standardises the KS-signed connectivity score for
   each compound to an empirical Z-score relative to the distribution of
   KS-signed scores across all compounds tested against the same trait
   (i.e. relative standing within the tested compound library, following
   current CMap/CLUE.io convention). A per-gene hypothesis-testing framework
   was evaluated during method development and rejected: genome-wide Pearson
   correlation between a single-compound and a polygenic disease signature
   is inherently small in magnitude, and gene-level non-independence makes a
   per-gene t-test-based p-value both statistically invalid and, empirically,
   uninformative at any reasonable corrected threshold (nothing survives
   FDR<0.25 even before correction is applied properly).

   A compound-trait pair is called a significant reversal candidate if its
   empirical Z-score < --z_threshold AND all five connectivity metrics are
   concordantly negative.

2. MOA enrichment: matches compounds to mechanism-of-action (MOA)
   annotations from the Broad Institute Drug Repurposing Hub
   (https://repo-hub.broadinstitute.org, data version 2020-03-24; a
   convenient reformatted mirror is available at
   https://github.com/ncats/drug_rep, data/broad_reformatted.txt), splits
   combination MOA annotations into elementary terms, and tests each term
   for enrichment among significant candidates via Fisher's exact test with
   Benjamini-Hochberg FDR correction.

Usage:
    python 03_significance_and_moa_enrichment.py \
        --connectivity_results cmap_drug_repurposing_results.tsv.gz \
        --moa_table broad_moa.txt \
        --z_threshold -2.0 \
        --output_significant cmap_significant_drugs.tsv \
        --output_moa_enrichment cmap_moa_enrichment.tsv
"""
import argparse

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

METRIC_COLS = ["pearson_all", "spearman_all", "ks_signed_mean", "extreme_pearson_mean", "extreme_spearman_mean"]


def call_significance(df, z_threshold):
    df = df.copy()
    df["ks_z"] = df.groupby("trait")["ks_signed_mean"].transform(lambda x: (x - x.mean()) / x.std())
    df["n_methods_negative"] = (df[METRIC_COLS] < 0).sum(axis=1)
    sig = df[(df.ks_z < z_threshold) & (df.n_methods_negative == len(METRIC_COLS))].copy()
    sig.sort_values(["trait", "ks_z"], inplace=True)
    return df, sig


def moa_enrichment(sig, all_df, moa_table_path, min_term_count=5):
    moa = pd.read_csv(moa_table_path, sep="\t", low_memory=False)
    moa_map = moa.drop_duplicates(subset="CompoundName")[["CompoundName", "MOA"]].copy()
    moa_map["compound_lower"] = moa_map.CompoundName.str.lower()

    all_df = all_df.copy()
    all_df["compound_lower"] = all_df.compound.str.lower()
    all_df = all_df.merge(moa_map[["compound_lower", "MOA"]], on="compound_lower", how="left")

    sig = sig.copy()
    sig["compound_lower"] = sig.compound.str.lower()
    sig = sig.merge(moa_map[["compound_lower", "MOA"]], on="compound_lower", how="left")

    has_moa = all_df[all_df.MOA.notna()].copy()
    sig_has_moa = sig[sig.MOA.notna()].copy()
    print(f"Background (all compound-trait pairs with MOA annotation): {len(has_moa):,}")
    print(f"Foreground (significant + MOA-annotated): {len(sig_has_moa):,}")

    # Split combination MOA annotations ("X|Y") into elementary terms
    has_moa_exp = has_moa.assign(MOA=has_moa.MOA.str.split("|")).explode("MOA")
    sig_exp = sig_has_moa.assign(MOA=sig_has_moa.MOA.str.split("|")).explode("MOA")

    results = []
    counts_bg = has_moa_exp.MOA.value_counts()
    for moa_term in counts_bg.index:
        if counts_bg[moa_term] < min_term_count:
            continue
        a = (sig_exp.MOA == moa_term).sum()
        b = len(sig_exp) - a
        c = counts_bg[moa_term] - a
        d = len(has_moa_exp) - a - b - c
        odds, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
        results.append(dict(MOA=moa_term, sig_count=a, total_count=counts_bg[moa_term], odds_ratio=odds, pvalue=p))

    enrich = pd.DataFrame(results)
    enrich["fdr"] = multipletests(enrich.pvalue, method="fdr_bh")[1]
    enrich.sort_values("pvalue", inplace=True)
    return enrich


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--connectivity_results", required=True)
    ap.add_argument("--moa_table", required=True, help="Broad Repurposing Hub annotation table")
    ap.add_argument("--z_threshold", type=float, default=-2.0)
    ap.add_argument("--output_significant", required=True)
    ap.add_argument("--output_all_with_stats", default=None)
    ap.add_argument("--output_moa_enrichment", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.connectivity_results, sep="\t")
    print(f"Total compound-trait pairs: {len(df):,}", flush=True)

    df_stats, sig = call_significance(df, args.z_threshold)
    print(f"Significant reversal candidates (Z<{args.z_threshold}, all 5 metrics concordant): {len(sig):,}")
    print(sig.groupby("trait").size())

    sig.to_csv(args.output_significant, sep="\t", index=False)
    if args.output_all_with_stats:
        df_stats.to_csv(args.output_all_with_stats, sep="\t", index=False, compression="gzip")

    print("\nRunning MOA enrichment...", flush=True)
    enrich = moa_enrichment(sig, df_stats, args.moa_table)
    enrich.to_csv(args.output_moa_enrichment, sep="\t", index=False)
    print(f"Tested {len(enrich)} elementary MOA terms")
    print(enrich.head(15).to_string(index=False))


if __name__ == "__main__":
    main()
