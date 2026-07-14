#!/usr/bin/env python3
"""
02_meta_analysis/meta_analysis_ivw.py

Memory-efficient fixed-effects inverse-variance-weighted (IVW) meta-analysis
of two GWAS summary statistics tables (as produced by
01_data_preparation/vcf_to_table.py), harmonising alleles between studies
and excluding palindromic (A/T, C/G) and orientation-ambiguous variants.

Designed for the case where one input file is substantially larger than the
other (e.g. FinnGen register-based cohorts with 15M+ variants): the smaller
study is loaded fully into memory and used to build an rsID lookup set, while
the larger study is streamed in chunks and filtered to the overlap on the
fly, keeping peak memory bounded regardless of the larger file's size.

Usage:
    python meta_analysis_ivw.py \
        --study1 tables/ieu-a-1185.table.gz \
        --study2 tables/finn-b-KRA_PSY_AUTISM_EXMORE.table.gz \
        --output meta_groupA_ASD.tsv.gz

Notes:
  - study1 should be the SMALLER of the two files for optimal memory use.
  - Input tables must have columns: variant_id, chromosome, position,
    effect_allele, non_effect_allele, effect_size, standard_error
    (this is exactly the schema produced by vcf_to_table.py).
"""
import argparse
import gc
import sys

import numpy as np
import pandas as pd
from scipy import stats

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
REQUIRED_COLS = ["variant_id", "chromosome", "position", "effect_allele",
                  "non_effect_allele", "effect_size", "standard_error"]


def load_study1(path):
    print(f"Loading study 1 (smaller, loaded fully): {path} ...", flush=True)
    s1 = pd.read_csv(path, sep="\t", usecols=REQUIRED_COLS,
                      dtype={"variant_id": "string", "chromosome": "str", "position": "int32",
                             "effect_allele": "str", "non_effect_allele": "str",
                             "effect_size": "float32", "standard_error": "float32"})
    s1.dropna(inplace=True)
    s1 = s1[(s1.effect_allele.str.len() == 1) & (s1.non_effect_allele.str.len() == 1)]
    s1.rename(columns={"effect_allele": "ea1", "non_effect_allele": "oa1",
                        "effect_size": "beta1", "standard_error": "se1"}, inplace=True)
    s1.set_index("variant_id", inplace=True)
    s1 = s1[~s1.index.duplicated(keep="first")]
    print(f"Study 1 usable SNPs: {len(s1):,}", flush=True)
    return s1


def stream_study2(path, id_set, chunksize=500_000):
    print(f"Streaming study 2 (chunked) and filtering to overlap: {path} ...", flush=True)
    chunks, n_seen = [], 0
    cols2 = ["variant_id", "effect_allele", "non_effect_allele", "effect_size", "standard_error"]
    reader = pd.read_csv(path, sep="\t", usecols=cols2,
                          dtype={"variant_id": "string", "effect_allele": "str",
                                 "non_effect_allele": "str", "effect_size": "float32",
                                 "standard_error": "float32"},
                          chunksize=chunksize)
    for i, chunk in enumerate(reader):
        n_seen += len(chunk)
        chunk = chunk[chunk.variant_id.isin(id_set)]
        chunk.dropna(inplace=True)
        if len(chunk):
            chunks.append(chunk)
        if (i + 1) % 6 == 0:
            kept = sum(len(c) for c in chunks)
            print(f"  ...processed {n_seen:,} study-2 rows, kept {kept:,} so far", flush=True)

    s2 = pd.concat(chunks, ignore_index=True)
    del chunks
    gc.collect()
    s2 = s2[(s2.effect_allele.str.len() == 1) & (s2.non_effect_allele.str.len() == 1)]
    s2.rename(columns={"effect_allele": "ea2", "non_effect_allele": "oa2",
                        "effect_size": "beta2", "standard_error": "se2"}, inplace=True)
    s2.set_index("variant_id", inplace=True)
    s2 = s2[~s2.index.duplicated(keep="first")]
    print(f"Study 2 filtered rows: {len(s2):,} (from {n_seen:,} total)", flush=True)
    return s2


def harmonise_and_meta_analyse(s1, s2):
    print("Merging on rsID index ...", flush=True)
    m = s1.join(s2, how="inner")
    del s1, s2
    gc.collect()
    print(f"Overlapping SNPs: {len(m):,}", flush=True)

    ea1 = m.ea1.to_numpy(dtype=str); oa1 = m.oa1.to_numpy(dtype=str)
    ea2 = m.ea2.to_numpy(dtype=str); oa2 = m.oa2.to_numpy(dtype=str)

    same = (ea1 == ea2) & (oa1 == oa2)
    swapped = (ea1 == oa2) & (oa1 == ea2)
    comp_oa1 = np.array([COMPLEMENT.get(x, "?") for x in oa1])
    palindromic = (comp_oa1 == ea1)  # A/T or C/G pairs -> strand-ambiguous, excluded

    n_same, n_swap = int(same.sum()), int(swapped.sum())
    n_mismatch = len(m) - n_same - n_swap
    print(f"Allele match: same={n_same:,} swapped={n_swap:,} mismatched(dropped)={n_mismatch:,}", flush=True)

    matched = same | swapped
    n_palin = int(palindromic[matched].sum())
    print(f"Palindromic SNPs among matched (excluded): {n_palin:,}", flush=True)

    keep = matched & (~palindromic)
    flip = swapped & keep

    m2 = m.loc[keep].copy()
    beta2_arr = m2.beta2.to_numpy(dtype=np.float64)
    beta2_arr[flip[keep]] *= -1
    ea1_kept, oa1_kept = ea1[keep], oa1[keep]
    print(f"Final harmonised SNP set: {len(m2):,}", flush=True)

    del m, ea1, oa1, ea2, oa2, same, swapped, comp_oa1, palindromic, matched, keep, flip
    gc.collect()

    w1 = 1.0 / (m2.se1.to_numpy(dtype=np.float32) ** 2)
    w2 = 1.0 / (m2.se2.to_numpy(dtype=np.float32) ** 2)
    beta1_arr = m2.beta1.to_numpy(dtype=np.float32)
    beta_meta = (w1 * beta1_arr + w2 * beta2_arr) / (w1 + w2)
    se_meta = np.sqrt(1.0 / (w1 + w2))
    z_meta = beta_meta / se_meta
    p_meta = 2 * stats.norm.sf(np.abs(z_meta))
    Q = w1 * (beta1_arr - beta_meta) ** 2 + w2 * (beta2_arr - beta_meta) ** 2
    I2 = np.maximum(0, (Q - 1) / Q) * 100
    het_p = stats.chi2.sf(Q, df=1)

    out = pd.DataFrame({
        "variant_id": m2.index.to_numpy(), "chromosome": m2.chromosome.to_numpy(),
        "position": m2.position.to_numpy(), "effect_allele": ea1_kept, "non_effect_allele": oa1_kept,
        "beta_meta": beta_meta, "se_meta": se_meta, "z_meta": z_meta, "p_meta": p_meta,
        "Q": Q, "I2": I2, "het_p": het_p,
        "beta_study1": beta1_arr, "se_study1": m2.se1.to_numpy(),
        "beta_study2": beta2_arr, "se_study2": m2.se2.to_numpy(),
    })
    out.sort_values("p_meta", inplace=True)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--study1", required=True, help="Smaller study table.gz (loaded fully)")
    ap.add_argument("--study2", required=True, help="Larger study table.gz (streamed in chunks)")
    ap.add_argument("--output", required=True, help="Output meta-analysis tsv.gz path")
    ap.add_argument("--chunksize", type=int, default=500_000)
    args = ap.parse_args()

    s1 = load_study1(args.study1)
    id_set = set(s1.index)
    print(f"Lookup set built: {len(id_set):,} rsIDs", flush=True)
    s2 = stream_study2(args.study2, id_set, args.chunksize)
    out = harmonise_and_meta_analyse(s1, s2)
    out.to_csv(args.output, sep="\t", index=False, compression="gzip")

    z_meta = out.z_meta.to_numpy()
    lam_meta = np.median(z_meta ** 2) / 0.4549
    n_gws = int((out.p_meta < 5e-8).sum())
    n_high_het = int((out.het_p < 0.05).sum())

    print("\n=== META-ANALYSIS SUMMARY ===", flush=True)
    print(f"Variants meta-analyzed: {len(out):,}")
    print(f"Genomic control lambda: {lam_meta:.3f}")
    print(f"Genome-wide significant (p<5e-8): {n_gws}")
    print(f"Nominally significant heterogeneity (het_p<0.05): {n_high_het:,} "
          f"({100 * n_high_het / len(out):.2f}%)")
    print("\nTop 10 hits:")
    print(out.head(10)[["variant_id", "chromosome", "position", "beta_meta", "se_meta", "p_meta", "I2"]]
          .to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
