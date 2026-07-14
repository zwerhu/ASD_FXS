#!/usr/bin/env python3
"""
04_transcriptome_imputation/recode_gwas_to_varid.py

Re-keys a GWAS summary statistics table from rsID to GTEx variant ID
(chr_pos_ref_alt_b38) using the mapping built by build_rsid_to_varid_map.py,
so that S-PrediXcan can be run with --model_db_snp_key varID (required for
the covariance lookup to work correctly; see build_rsid_to_varid_map.py
docstring for why this is necessary).

Processes the input file in chunks to bound memory usage regardless of file
size.

Usage:
    python recode_gwas_to_varid.py --study ieu-a-1185 \
        --input_dir tables --output_dir tables_recoded \
        --mapping rsid_to_varid_map.tsv.gz
"""
import argparse
import gc

import pandas as pd

OUTPUT_COLS = ["variant_id", "chromosome", "position", "effect_allele", "non_effect_allele",
               "pvalue", "zscore", "effect_size", "standard_error"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--study", required=True, help="Study short name (expects <input_dir>/<study>.table.gz)")
    ap.add_argument("--input_dir", default="tables")
    ap.add_argument("--output_dir", default="tables_recoded")
    ap.add_argument("--mapping", required=True, help="rsid_to_varid_map.tsv.gz from build_rsid_to_varid_map.py")
    ap.add_argument("--chunksize", type=int, default=1_000_000)
    args = ap.parse_args()

    mapping = pd.read_csv(args.mapping, sep="\t", usecols=["rsid", "varID"])
    mapping.set_index("rsid", inplace=True)
    mapping = mapping[~mapping.index.duplicated(keep="first")]

    chunks, n_seen = [], 0
    in_path = f"{args.input_dir}/{args.study}.table.gz"
    for chunk in pd.read_csv(in_path, sep="\t", chunksize=args.chunksize, dtype={"chromosome": "str"}):
        n_seen += len(chunk)
        chunk = chunk.join(mapping, on="variant_id", how="inner")
        if len(chunk):
            chunk = chunk.drop(columns=["variant_id"]).rename(columns={"varID": "variant_id"})[OUTPUT_COLS]
            chunks.append(chunk)

    out = pd.concat(chunks, ignore_index=True)
    del chunks
    gc.collect()

    out_path = f"{args.output_dir}/{args.study}.recoded.table.gz"
    out.to_csv(out_path, sep="\t", index=False, compression="gzip")
    print(f"{args.study}: {n_seen:,} -> {len(out):,} recoded rows -> {out_path}", flush=True)


if __name__ == "__main__":
    main()
