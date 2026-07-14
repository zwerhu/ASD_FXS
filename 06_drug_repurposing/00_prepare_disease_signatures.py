#!/usr/bin/env python3
"""
06_drug_repurposing/00_prepare_disease_signatures.py

Builds a gene x trait matrix of mean cross-tissue Z-scores from the
cross-tissue ACAT-integrated S-PrediXcan results, restricted to the gene
universe shared with the CMap L1000 "best-inferred" gene space. This matrix
is the "disease signature" input consumed by 02_connectivity_analysis.py
and 04_drug_disease_gene_tissue_links.py.

Usage:
    python 00_prepare_disease_signatures.py \
        --integrative_results integrative_ACAT_results.tsv.gz \
        --cmap_gene_info GSE92742_Broad_LINCS_gene_info.txt.gz \
        --output disease_signatures_for_cmap.tsv.gz \
        --gene_list_output overlap_gene_list.txt
"""
import argparse

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--integrative_results", required=True)
    ap.add_argument("--cmap_gene_info", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--gene_list_output", default=None,
                     help="Optional: also write the overlap gene list, one symbol per line "
                          "(useful for pre-filtering CMap data on a separate machine)")
    args = ap.parse_args()

    gene_info = pd.read_csv(args.cmap_gene_info, sep="\t")
    cmap_genes = set(gene_info.pr_gene_symbol)
    print(f"CMap gene universe: {len(cmap_genes):,} genes", flush=True)

    integ = pd.read_csv(args.integrative_results, sep="\t")
    print(f"Integrative results: {integ.trait.nunique()} traits, {integ.gene_name.nunique():,} genes", flush=True)

    mat = integ.pivot_table(index="gene_name", columns="trait", values="mean_zscore")
    overlap_genes = sorted(set(mat.index) & cmap_genes)
    mat_overlap = mat.loc[overlap_genes]
    print(f"\nOverlap with CMap gene space: {len(mat_overlap):,} / {len(mat):,} genes", flush=True)
    for trait in mat.columns:
        n_valid = mat_overlap[trait].notna().sum()
        print(f"  {trait}: {n_valid:,} genes with valid disease Z-score in CMap space")

    mat_overlap.to_csv(args.output, sep="\t", compression="gzip")
    print(f"\nSaved: {args.output}")

    if args.gene_list_output:
        with open(args.gene_list_output, "w") as f:
            f.write("\n".join(overlap_genes))
        print(f"Saved gene list: {args.gene_list_output} ({len(overlap_genes)} genes)")


if __name__ == "__main__":
    main()
