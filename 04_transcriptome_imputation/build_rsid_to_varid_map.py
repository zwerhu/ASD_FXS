#!/usr/bin/env python3
"""
04_transcriptome_imputation/build_rsid_to_varid_map.py

GTEx v8 MASHR prediction models and their SNP covariance matrices are keyed
by a GRCh38 "variant ID" (chr_pos_ref_alt_b38), not by rsID -- even though
the covariance file's columns are labelled RSID1/RSID2. If a GWAS is in
GRCh37 coordinates and indexed by rsID (as is typical for older summary
statistics releases), naively loading the model with rsID as the join key
will match the SNP *weights* table (which does carry rsIDs) but silently
fail to match the *covariance* table, producing 0% SNP usage without error.

This script avoids needing an explicit genome-build liftover by building a
rsID -> varID translation table directly from the union of all downloaded
tissue models' own weight tables, which conveniently carry both identifiers
for the same SNP.

Usage:
    python build_rsid_to_varid_map.py --models_dir mashr_models/mashr \
                                       --output rsid_to_varid_map.tsv.gz
"""
import argparse
import glob
import sqlite3

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models_dir", required=True, help="Directory containing mashr_<tissue>.db files")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    dbs = sorted(glob.glob(f"{args.models_dir}/*.db"))
    print(f"{len(dbs)} tissue model files found", flush=True)

    frames = []
    for db in dbs:
        con = sqlite3.connect(db)
        w = pd.read_sql("SELECT rsid, varID, ref_allele, eff_allele FROM weights WHERE rsid LIKE 'rs%'", con)
        con.close()
        frames.append(w)

    allw = pd.concat(frames, ignore_index=True)
    print(f"Total rows before dedup: {len(allw):,}", flush=True)
    allw = allw.drop_duplicates(subset=["rsid"])
    print(f"Unique rsID entries: {len(allw):,}", flush=True)

    allw.to_csv(args.output, sep="\t", index=False, compression="gzip")
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
