#!/usr/bin/env python3
"""
01_data_preparation/batch_convert_vcf.py

Batch-converts a directory of GWAS-VCF(.vcf.gz) summary statistics files
into gzipped TSV tables using vcf_to_table.py, skipping files that have
already been converted. Designed to be resumable, since individual
conversions of large files (10-20 million variants) can take several
minutes each.

Usage:
    python batch_convert_vcf.py --input_dir /path/to/vcf_files \
                                 --output_dir /path/to/tables \
                                 --converter vcf_to_table.py
"""
import argparse
import glob
import os
import subprocess
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_dir", required=True, help="Directory containing *_vcf.gz files")
    ap.add_argument("--output_dir", required=True, help="Directory to write *.table.gz outputs")
    ap.add_argument("--converter", default="vcf_to_table.py", help="Path to vcf_to_table.py")
    ap.add_argument("--pattern", default="*_vcf.gz", help="Glob pattern for input files")
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    files = sorted(glob.glob(os.path.join(args.input_dir, args.pattern)))
    if not files:
        sys.exit(f"No files matching {args.pattern} found in {args.input_dir}")

    print(f"Found {len(files)} VCF files to convert")
    for i, vcf_path in enumerate(files, 1):
        base = os.path.basename(vcf_path).replace("_vcf.gz", "").replace(".vcf.gz", "")
        out_path = os.path.join(args.output_dir, f"{base}.table.gz")
        if os.path.exists(out_path):
            print(f"[{i}/{len(files)}] SKIP (exists): {base}")
            continue
        print(f"[{i}/{len(files)}] Converting: {base}")
        subprocess.run([sys.executable, args.converter, vcf_path, out_path], check=True)
    print("All conversions complete.")


if __name__ == "__main__":
    main()
