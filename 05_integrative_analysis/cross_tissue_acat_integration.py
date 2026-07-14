#!/usr/bin/env python3
"""
05_integrative_analysis/cross_tissue_acat_integration.py

Combines S-PrediXcan gene-tissue-trait association results across GTEx
tissues into a single per-gene-per-trait statistic using the Aggregated
Cauchy Association Test (ACAT; Liu et al. 2019, Am J Hum Genet), which is
robust to arbitrary and unknown correlation among the combined p-values --
a necessary property given the substantial cross-tissue correlation inherent
to GTEx expression prediction models trained on overlapping reference
genotypes. Naive alternatives (Fisher's method, Bonferroni-across-tissues)
assume independence and are not appropriate here.

Genes within the extended MHC region are excluded via a cytogenetic-band
proxy (chr6p21), since S-PrediXcan gene-level output does not carry genomic
coordinates.

Usage:
    python cross_tissue_acat_integration.py \
        --input spredixcan_ALL_RESULTS.tsv.gz \
        --mhc_gene_list mhc_gene_list.txt \
        --output integrative_ACAT_results.tsv.gz \
        --trait_column trait
"""
import argparse

import numpy as np
import pandas as pd


def acat(pvals):
    """Aggregated Cauchy Association Test: combine p-values, robust to correlation."""
    pvals = np.clip(np.asarray(pvals, dtype=np.float64), 1e-300, 1 - 1e-16)
    w = 1.0 / len(pvals)
    stat = np.sum(w * np.tan((0.5 - pvals) * np.pi))
    return 0.5 - np.arctan(stat) / np.pi


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="Combined S-PrediXcan results (all tissues x traits)")
    ap.add_argument("--mhc_gene_list", required=True, help="Newline-delimited gene symbol list to exclude")
    ap.add_argument("--output", required=True)
    ap.add_argument("--trait_column", default="trait")
    args = ap.parse_args()

    print("Loading combined S-PrediXcan results...", flush=True)
    df = pd.read_csv(args.input, sep="\t")
    trait_col = args.trait_column
    print(f"Loaded: {len(df):,} rows, {df[trait_col].nunique()} traits, {df.tissue.nunique()} tissues", flush=True)

    with open(args.mhc_gene_list) as f:
        mhc_genes = set(l.strip() for l in f)

    before = len(df)
    df = df[~df.gene_name.isin(mhc_genes)]
    print(f"Excluded MHC region (chr6p21 band proxy): {before - len(df):,} rows removed "
          f"({df.gene_name.nunique()} genes remain)", flush=True)

    df = df.dropna(subset=["pvalue", "zscore"])
    df = df[(df.pvalue > 0) & (df.pvalue <= 1)]
    print(f"After dropping NA/invalid p-values: {len(df):,} rows", flush=True)

    print("Combining p-values across tissues per (trait, gene) via ACAT ...", flush=True)
    records = []
    for (trait, gene), g in df.groupby([trait_col, "gene_name"]):
        p_acat = acat(g.pvalue.values)
        mean_z = g.zscore.mean()
        n_tissues_tested = g.tissue.nunique()
        best_row = g.loc[g.pvalue.idxmin()]
        records.append({
            "trait": trait, "gene_name": gene, "gene": best_row.gene,
            "p_acat": p_acat, "mean_zscore": mean_z, "n_tissues_tested": n_tissues_tested,
            "best_tissue": best_row.tissue, "best_p": best_row.pvalue, "best_zscore": best_row.zscore
        })

    integ = pd.DataFrame(records)
    integ.sort_values("p_acat", inplace=True)
    integ.to_csv(args.output, sep="\t", index=False, compression="gzip")

    print(f"\nIntegrated table: {len(integ):,} gene-trait combinations", flush=True)
    for trait in integ.trait.unique():
        sub = integ[integ.trait == trait]
        n_sig = (sub.p_acat < 0.05 / len(sub)).sum()
        print(f"  {trait}: {len(sub):,} genes tested, {n_sig} Bonferroni-significant (ACAT-combined)")

    print("\nTop 20 overall:")
    print(integ.head(20)[["trait", "gene_name", "p_acat", "mean_zscore", "n_tissues_tested", "best_tissue"]]
          .to_string(index=False))


if __name__ == "__main__":
    main()
