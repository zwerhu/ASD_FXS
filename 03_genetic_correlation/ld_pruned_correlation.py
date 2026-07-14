#!/usr/bin/env python3
"""
03_genetic_correlation/ld_pruned_correlation.py

LD-reference-free proxy for genetic correlation between a meta-analysed GWAS
trait and a second GWAS trait (e.g. an MZ-twin symptom-discordance GWAS),
for use when a proper LD Score Regression reference panel is not available.

IMPORTANT: this is NOT a replacement for LDSC. It is a pragmatic proxy
that should be interpreted alongside, not instead of, a properly powered
bivariate LDSC analysis whenever one can be run.

Method:
  1. Harmonise alleles between the two GWAS (same/swapped orientation
     matching; palindromic A/T and C/G variants excluded as strand-ambiguous).
  2. Distance-based pruning: keep one SNP per 1 Mb window per chromosome,
     selected independently of either trait's significance (first SNP by
     position) to avoid winner's-curse-style bias, as an approximate
     LD-independence set.
  3. Report Pearson correlation and sign concordance (binomial test) of
     Z-scores at the pruned, approximately independent SNP set.
  4. Also report an UNPRUNED analysis restricted to the reference trait's own
     suggestive hits (p < 1e-5) for comparison -- this is included
     deliberately to illustrate how dramatically unpruned analyses can be
     inflated by linkage disequilibrium (see independent_lead_loci_test.py
     for the corrected version of this specific test).

Usage:
    python ld_pruned_correlation.py \
        --reference meta_groupA_ASD.tsv.gz \
        --comparison tables/ieu-b-5151.table.gz \
        --label ieu-b-5151 \
        --output_prefix corr_pruned
"""
import argparse
import gc

import numpy as np
import pandas as pd
from scipy import stats

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def load_reference(path):
    meta = pd.read_csv(path, sep="\t",
                        usecols=["variant_id", "chromosome", "position", "effect_allele",
                                  "non_effect_allele", "z_meta", "p_meta"],
                        dtype={"variant_id": "string", "chromosome": "str", "position": "int32",
                               "effect_allele": "str", "non_effect_allele": "str",
                               "z_meta": "float32", "p_meta": "float64"})
    meta.set_index("variant_id", inplace=True)
    meta = meta[~meta.index.duplicated(keep="first")]
    return meta


def stream_comparison(path, id_set, chunksize=500_000):
    chunks = []
    reader = pd.read_csv(path, sep="\t",
                          usecols=["variant_id", "effect_allele", "non_effect_allele", "zscore"],
                          dtype={"variant_id": "string", "effect_allele": "str",
                                 "non_effect_allele": "str", "zscore": "float32"},
                          chunksize=chunksize)
    n_seen = 0
    for chunk in reader:
        n_seen += len(chunk)
        chunk = chunk[chunk.variant_id.isin(id_set)]
        chunk.dropna(inplace=True)
        if len(chunk):
            chunks.append(chunk)
    comp = pd.concat(chunks, ignore_index=True)
    del chunks
    gc.collect()
    comp.rename(columns={"effect_allele": "ea_t", "non_effect_allele": "oa_t", "zscore": "z_comp"}, inplace=True)
    comp.set_index("variant_id", inplace=True)
    comp = comp[~comp.index.duplicated(keep="first")]
    print(f"Comparison-file rows kept after overlap filter: {len(comp):,} (of {n_seen:,} total)", flush=True)
    return comp


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", required=True, help="Meta-analysis tsv.gz (output of meta_analysis_ivw.py)")
    ap.add_argument("--comparison", required=True, help="Comparison trait table.gz")
    ap.add_argument("--label", required=True, help="Short label for output filenames")
    ap.add_argument("--output_prefix", default="corr_pruned")
    ap.add_argument("--window_bp", type=int, default=1_000_000, help="Pruning window size in bp")
    ap.add_argument("--suggestive_p", type=float, default=1e-5)
    args = ap.parse_args()

    print(f"=== {args.label}: reference vs comparison trait ===", flush=True)
    print("Loading reference (meta-analysis) trait ...", flush=True)
    meta = load_reference(args.reference)
    print(f"Reference SNPs: {len(meta):,}", flush=True)
    id_set = set(meta.index)

    print(f"Streaming comparison trait: {args.comparison} ...", flush=True)
    comp = stream_comparison(args.comparison, id_set)

    m = meta.join(comp, how="inner")
    del comp
    gc.collect()
    print(f"Merged overlap: {len(m):,}", flush=True)

    ea1 = m.effect_allele.to_numpy(dtype=str); oa1 = m.non_effect_allele.to_numpy(dtype=str)
    ea2 = m.ea_t.to_numpy(dtype=str); oa2 = m.oa_t.to_numpy(dtype=str)
    same = (ea1 == ea2) & (oa1 == oa2)
    swapped = (ea1 == oa2) & (oa1 == ea2)
    comp_oa1 = np.array([COMPLEMENT.get(x, "?") for x in oa1])
    palindromic = (comp_oa1 == ea1)
    matched = same | swapped
    keep = matched & (~palindromic)
    flip = swapped & keep
    print(f"Allele match: same={same.sum():,} swapped={swapped.sum():,} "
          f"mismatched(dropped)={(len(m) - matched.sum()):,} "
          f"palindromic(dropped)={(palindromic & matched).sum():,}", flush=True)

    m2 = m.loc[keep].copy()
    z_comp = m2.z_comp.to_numpy(dtype=np.float64)
    z_comp[flip[keep]] *= -1
    m2["z_comp_h"] = z_comp
    print(f"Final harmonised overlap set: {len(m2):,}", flush=True)

    # ---- Distance-based pruning ----
    m2 = m2.reset_index()
    m2.sort_values(["chromosome", "position"], inplace=True)
    m2["bin"] = m2.position // args.window_bp
    pruned = m2.groupby(["chromosome", "bin"], as_index=False).first()
    print(f"Pruned (approx-independent) SNP count: {len(pruned):,}", flush=True)

    r, r_p = stats.pearsonr(pruned.z_meta, pruned.z_comp_h)
    n_concord = int((np.sign(pruned.z_meta) == np.sign(pruned.z_comp_h)).sum())
    n_tot = len(pruned)
    sign_p = stats.binomtest(n_concord, n_tot, 0.5).pvalue

    print("\n--- Genome-wide pruned-SNP results (approximately LD-independent) ---")
    print(f"Pruned SNPs: {n_tot:,}")
    print(f"Pearson r (Z-score concordance): {r:.4f}  (p={r_p:.3g})")
    print(f"Sign concordance: {n_concord:,}/{n_tot:,} = {100 * n_concord / n_tot:.2f}%  (binomial p={sign_p:.3g})")

    # ---- Naive unpruned top-loci check (included for comparison / caution) ----
    top = m2[m2.p_meta < args.suggestive_p]
    if len(top) >= 3:
        r_top, r_top_p = stats.pearsonr(top.z_meta, top.z_comp_h)
        n_c_top = int((np.sign(top.z_meta) == np.sign(top.z_comp_h)).sum())
        sign_p_top = stats.binomtest(n_c_top, len(top), 0.5).pvalue
        print(f"\n--- CAUTION: reference trait's suggestive loci (p<{args.suggestive_p}), UNPRUNED, N={len(top)} ---")
        print(f"Pearson r: {r_top:.4f} (p={r_top_p:.3g})")
        print(f"Sign concordance: {n_c_top}/{len(top)} = {100 * n_c_top / len(top):.2f}% (binomial p={sign_p_top:.3g})")
        print("NOTE: this unpruned result can be severely inflated by linkage disequilibrium if suggestive")
        print("loci cluster in a small number of LD blocks. Always cross-check with")
        print("independent_lead_loci_test.py before interpreting this number.")
    else:
        print(f"\nToo few suggestive loci overlapping ({len(top)}) for a meaningful targeted test")

    out_cols = ["variant_id", "chromosome", "position", "z_meta", "z_comp_h", "p_meta"]
    pruned[out_cols].to_csv(f"{args.output_prefix}_{args.label}.tsv.gz", sep="\t", index=False, compression="gzip")
    print(f"\nSaved pruned comparison set: {args.output_prefix}_{args.label}.tsv.gz")


if __name__ == "__main__":
    main()
