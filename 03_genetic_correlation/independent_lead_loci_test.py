#!/usr/bin/env python3
"""
03_genetic_correlation/independent_lead_loci_test.py

Corrected, LD-aware version of the "reference trait's suggestive loci"
concordance test in ld_pruned_correlation.py.

The naive unpruned test in ld_pruned_correlation.py can produce dramatically
inflated apparent concordance when a trait's suggestive hits are dominated by
a small number of linkage-disequilibrium blocks -- in the accompanying
manuscript, an apparently overwhelming concordance signal (p < 1e-12) using
all suggestive SNPs collapsed to a null result (p > 0.4) once restricted to
genuinely independent lead loci. This script performs that correction:

  1. Extract the reference trait's suggestive loci (p < threshold).
  2. Collapse to one LEAD SNP per 1 Mb window (the lowest-p SNP in each
     window), yielding a genuinely independent locus set.
  3. Look up each lead SNP's Z-score in the comparison trait, harmonising
     alleles.
  4. Report sign concordance (binomial test) and Pearson correlation at
     this reduced, non-redundant locus set.

Usage:
    python independent_lead_loci_test.py \
        --reference meta_groupA_ASD.tsv.gz \
        --comparison tables/ieu-b-5151.table.gz \
        --label ieu-b-5151 \
        --suggestive_p 1e-5 --window_bp 1000000
"""
import argparse

import numpy as np
import pandas as pd
from scipy import stats


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--comparison", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--suggestive_p", type=float, default=1e-5)
    ap.add_argument("--window_bp", type=int, default=1_000_000)
    ap.add_argument("--output", default=None, help="Optional path to save the lead-loci table")
    args = ap.parse_args()

    meta = pd.read_csv(args.reference, sep="\t",
                        usecols=["variant_id", "chromosome", "position", "p_meta", "z_meta",
                                  "effect_allele", "non_effect_allele"],
                        dtype={"chromosome": "str"})
    top = meta[meta.p_meta < args.suggestive_p].copy()
    top["bin"] = top.position // args.window_bp
    lead = top.sort_values("p_meta").groupby(["chromosome", "bin"], as_index=False).first()
    print(f"Lead independent loci (p<{args.suggestive_p}): {len(lead)}", flush=True)

    if args.output:
        lead.to_csv(args.output, sep="\t", index=False)

    ids = set(lead.variant_id)
    rows = []
    for chunk in pd.read_csv(args.comparison, sep="\t", chunksize=1_000_000,
                              usecols=["variant_id", "effect_allele", "non_effect_allele", "zscore"]):
        rows.append(chunk[chunk.variant_id.isin(ids)])
    comp = pd.concat(rows, ignore_index=True)

    merged = lead.merge(comp, on="variant_id", suffixes=("_ref", "_comp"))
    same = merged.effect_allele_ref == merged.effect_allele_comp
    swapped = merged.effect_allele_ref == merged.non_effect_allele_comp
    z_comp = merged.zscore.copy()
    z_comp[swapped] = -z_comp[swapped]
    merged["z_comp_h"] = z_comp
    merged = merged[same | swapped]

    n = len(merged)
    if n < 3:
        print(f"Too few independent lead loci with comparison-trait data (N={n}) for a meaningful test")
        return

    n_concord = int((np.sign(merged.z_meta) == np.sign(merged.z_comp_h)).sum())
    p_sign = stats.binomtest(n_concord, n, 0.5).pvalue
    r, r_p = stats.pearsonr(merged.z_meta, merged.z_comp_h)

    print(f"\n=== {args.label}: independent lead loci (N={n}) ===")
    print(merged[["variant_id", "chromosome", "position", "z_meta", "z_comp_h"]].to_string(index=False))
    print(f"Sign concordance: {n_concord}/{n} = {100 * n_concord / n:.1f}%  binomial p={p_sign:.3g}")
    print(f"Pearson r: {r:.3f}  p={r_p:.3g}")


if __name__ == "__main__":
    main()
