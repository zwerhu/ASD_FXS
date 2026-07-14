#!/usr/bin/env python3
"""
06_drug_repurposing/merge_filter_gctx.py

Local pre-processing script for the CMap LINCS L1000 Level 5 dataset
(GSE92742). Intended to be run on a machine with sufficient local disk space
(the merged, decompressed Level 5 GCTX file is ~20-25 GB), since the full
dataset is far too large to reassemble in a disk-constrained analysis
sandbox.

Two-step pipeline:
  1. merge_and_decompress: streams a set of gzip-split file chunks
     (level5_part_aa, level5_part_ab, ... as produced by `split` on the
     original .gctx.gz download) directly into a decompressed .gctx file,
     using a constant-memory streaming zlib decompressor -- never holding
     the full compressed or decompressed byte stream in memory at once.
  2. filter: uses cmapPy's subsetted GCTX reader (which avoids loading the
     full 473,648 x 12,328 matrix into memory) to filter to:
       - pert_type == "trt_cp" (small-molecule compound treatments)
       - is_exemplar == 1 (QC-selected representative signature)
       - optionally, a specific set of cell lines (e.g. ["NPC", "NEU"] for
         a neural-lineage-restricted analysis)
     then aggregates to one representative signature per compound (highest
     transcriptional activity score, "tas"), and optionally restricts genes
     to a supplied gene list, before writing a small, portable TSV matrix.

Usage:
    python merge_filter_gctx.py \
        --parts_dir . --parts_pattern "level5_part_*" \
        --merged_gctx level5_merged.gctx \
        --sig_info GSE92742_Broad_LINCS_sig_info.txt.gz \
        --sig_metrics GSE92742_Broad_LINCS_sig_metrics.txt.gz \
        --gene_info GSE92742_Broad_LINCS_gene_info.txt.gz \
        --gene_list overlap_gene_list.txt \
        --output_matrix level5_filtered_matrix.tsv.gz \
        --cell_lines NPC NEU        # omit this flag to keep all cell lines
"""
import argparse
import glob
import os
import sys
import zlib


def merge_and_decompress(parts_dir, parts_pattern, merged_gctx):
    if os.path.exists(merged_gctx):
        print(f"[SKIP] {merged_gctx} already exists, skipping merge/decompress step")
        return

    parts = sorted(glob.glob(os.path.join(parts_dir, parts_pattern)))
    if not parts:
        sys.exit(f"No files matching {parts_pattern} found in {parts_dir}")
    print(f"Found {len(parts)} chunks, first/last:")
    for p in parts[:3]:
        print(" ", p)
    print("  ...")
    for p in parts[-3:]:
        print(" ", p)

    CHUNK = 8 * 1024 * 1024  # 8MB read/write blocks

    def feed_stream():
        for p in parts:
            with open(p, "rb") as fin:
                while True:
                    buf = fin.read(CHUNK)
                    if not buf:
                        break
                    yield buf

    def gzip_stream_decompress(byte_iter):
        # wbits = 16 + MAX_WBITS decodes the gzip container format (not raw deflate)
        d = zlib.decompressobj(16 + zlib.MAX_WBITS)
        for chunk in byte_iter:
            out = d.decompress(chunk)
            if out:
                yield out
        tail = d.flush()
        if tail:
            yield tail

    total_out = 0
    with open(merged_gctx, "wb") as fout:
        for out_chunk in gzip_stream_decompress(feed_stream()):
            fout.write(out_chunk)
            total_out += len(out_chunk)
            if total_out % (500 * 1024 * 1024) < len(out_chunk):
                print(f"  ...decompressed {total_out / 1e9:.2f} GB")

    print(f"Merge/decompress complete: {merged_gctx} ({total_out / 1e9:.2f} GB)")


def filter_and_aggregate(args):
    import pandas as pd
    from cmapPy.pandasGEXpress.parse_gctx import parse

    print("Reading signature metadata and determining target signature list...")
    sig_info = pd.read_csv(args.sig_info, sep="\t", low_memory=False)
    sig_metrics = pd.read_csv(args.sig_metrics, sep="\t")
    m = sig_info.merge(sig_metrics[["sig_id", "is_exemplar", "tas"]], on="sig_id", how="left")
    cand = m[(m.pert_type == "trt_cp") & (m.is_exemplar == 1)].copy()
    print(f"Candidate signatures (trt_cp + exemplar): {len(cand):,}")

    if args.cell_lines:
        before = len(cand)
        cand = cand[cand.cell_id.isin(args.cell_lines)]
        print(f"Filtered to cell lines {args.cell_lines}: {before:,} -> {len(cand):,} signatures "
              f"(covering {cand.pert_iname.nunique():,} compounds)")

    # ---- Aggregate to one representative signature per compound (highest TAS) ----
    cand = cand.sort_values("tas", ascending=False)
    rep = cand.drop_duplicates(subset="pert_iname", keep="first")
    print(f"After per-compound aggregation (highest-TAS representative): {len(rep):,} compounds")
    target_ids = rep.sig_id.tolist()
    sigid_to_name = dict(zip(rep.sig_id, rep.pert_iname))

    rid = None
    if args.gene_list and os.path.exists(args.gene_list):
        with open(args.gene_list) as f:
            gene_symbols = set(l.strip() for l in f if l.strip())
        print(f"Also filtering to {len(gene_symbols)} genes (requires gene_info for symbol->pr_gene_id mapping)")
        if args.gene_info and os.path.exists(args.gene_info):
            gi = pd.read_csv(args.gene_info, sep="\t")
            rid = gi.loc[gi.pr_gene_symbol.isin(gene_symbols), "pr_gene_id"].astype(str).tolist()
            print(f"Mapped to {len(rid)} row IDs")
        else:
            print(f"[WARNING] gene_info file not found ({args.gene_info}), skipping gene filter, keeping all genes")

    print("Reading GCTX subset via cmapPy (does not load the full matrix into memory)...")
    gctoo = parse(args.merged_gctx, cid=target_ids, rid=rid)
    print(f"Read complete, matrix shape: {gctoo.data_df.shape}")

    df = gctoo.data_df.copy()
    df.columns = [sigid_to_name.get(c, c) for c in df.columns]
    if args.gene_info and os.path.exists(args.gene_info):
        gi = pd.read_csv(args.gene_info, sep="\t", dtype={"pr_gene_id": str})
        id_to_symbol = dict(zip(gi.pr_gene_id, gi.pr_gene_symbol))
        df.index = [id_to_symbol.get(i, i) for i in df.index]

    if args.output_gctx:
        from cmapPy.pandasGEXpress.write_gctx import write
        write(gctoo, args.output_gctx)
        print(f"Wrote filtered GCTX: {args.output_gctx}")

    df.to_csv(args.output_matrix, sep="\t", compression="gzip")
    print(f"Wrote filtered matrix: {args.output_matrix} ({df.shape[0]} genes x {df.shape[1]} compounds)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--parts_dir", default=".")
    ap.add_argument("--parts_pattern", default="level5_part_*")
    ap.add_argument("--merged_gctx", default="level5_merged.gctx")
    ap.add_argument("--sig_info", required=True)
    ap.add_argument("--sig_metrics", required=True)
    ap.add_argument("--gene_info", default=None)
    ap.add_argument("--gene_list", default=None, help="Optional newline-delimited gene symbol allowlist")
    ap.add_argument("--output_matrix", required=True)
    ap.add_argument("--output_gctx", default=None, help="Optional path to also save the filtered GCTX")
    ap.add_argument("--cell_lines", nargs="*", default=None,
                     help="Restrict to these cell_id values (e.g. NPC NEU for a neural-lineage-only "
                          "analysis). Omit to keep all cell lines.")
    args = ap.parse_args()

    merge_and_decompress(args.parts_dir, args.parts_pattern, args.merged_gctx)
    filter_and_aggregate(args)


if __name__ == "__main__":
    main()
