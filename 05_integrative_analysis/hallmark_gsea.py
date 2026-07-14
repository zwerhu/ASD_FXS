#!/usr/bin/env python3
"""
05_integrative_analysis/hallmark_gsea.py

Runs pre-ranked Gene Set Enrichment Analysis (GSEA; Subramanian et al. 2005)
against the 50 MSigDB Hallmark gene sets for each trait in the cross-tissue
ACAT-integrated results table, ranking genes by signed, ACAT-combined
significance (-log10(p) x sign of the mean cross-tissue Z-score).

Gene sets are loaded from a local JSON file containing MSigDB v7.1 gene sets
released under CC-BY (Hallmark + GO collections; KEGG/BioCarta are excluded
from this mirror for licensing reasons). A convenient CC-BY mirror is
available at: https://github.com/numpde/genesets
(genesets/msigdb/parsed/v7.1/genesets.json.zip)

Usage:
    python hallmark_gsea.py \
        --integrative_results integrative_ACAT_results.tsv.gz \
        --genesets_json genesets.json \
        --output_dir gsea_out \
        --combined_output gsea_ALL_hallmark_results.tsv.gz
"""
import argparse
import json
import os

import gseapy as gp
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--integrative_results", required=True)
    ap.add_argument("--genesets_json", required=True)
    ap.add_argument("--output_dir", default="gsea_out")
    ap.add_argument("--combined_output", required=True)
    ap.add_argument("--permutation_num", type=int, default=1000)
    ap.add_argument("--fdr_report_threshold", type=float, default=0.25)
    args = ap.parse_args()

    integ = pd.read_csv(args.integrative_results, sep="\t")

    with open(args.genesets_json) as f:
        gs_all = json.load(f)
    hallmark = {k: v["symbols"] for k, v in gs_all.items() if k.startswith("HALLMARK")}
    print(f"Loaded {len(hallmark)} Hallmark gene sets", flush=True)

    os.makedirs(args.output_dir, exist_ok=True)

    all_results = []
    for trait in integ.trait.unique():
        sub = integ[integ.trait == trait].copy()
        sub["rank_metric"] = -np.log10(sub.p_acat) * np.sign(sub.mean_zscore)
        sub = sub.dropna(subset=["rank_metric"])
        sub = sub.sort_values("rank_metric", ascending=False)
        rnk = sub[["gene_name", "rank_metric"]].drop_duplicates(subset="gene_name")
        print(f"\n=== {trait}: {len(rnk)} genes ranked ===", flush=True)

        try:
            pre_res = gp.prerank(rnk=rnk, gene_sets=hallmark, min_size=5, max_size=1000,
                                  permutation_num=args.permutation_num, outdir=None, seed=42,
                                  threads=1, verbose=False)
            res_df = pre_res.res2d
            res_df["trait"] = trait
            all_results.append(res_df)
            res_df.to_csv(f"{args.output_dir}/{trait}_hallmark_gsea.tsv", sep="\t", index=False)
            sig = res_df[res_df["FDR q-val"].astype(float) < args.fdr_report_threshold].sort_values("FDR q-val")
            print(f"  Pathways FDR<{args.fdr_report_threshold}: {len(sig)}", flush=True)
            if len(sig):
                print(sig[["Term", "NES", "NOM p-val", "FDR q-val"]].head(10).to_string(index=False), flush=True)
        except Exception as e:
            print(f"  FAILED for {trait}: {e}", flush=True)

    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        combined.to_csv(args.combined_output, sep="\t", index=False, compression="gzip")
        print(f"\nCombined GSEA results: {len(combined)} rows saved to {args.combined_output}", flush=True)


if __name__ == "__main__":
    main()
