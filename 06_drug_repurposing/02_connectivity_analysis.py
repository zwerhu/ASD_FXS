#!/usr/bin/env python3
"""
06_drug_repurposing/02_connectivity_analysis.py

Computes five complementary connectivity metrics between disease/trait
transcriptomic signatures and a library of CMap compound signatures,
following the classical Connectivity Map methodology (Lamb et al. 2006) and
its extensions:

  1. Pearson correlation, all genes.
  2. Spearman correlation, all genes.
  3. Signed Kolmogorov-Smirnov (KS) connectivity score, computed from the
     top/bottom N disease genes' positions within the rank-ordered compound
     signature (averaged across gene-set sizes N = 50, 100, 250, 500).
  4. "Extreme" Pearson correlation restricted to the top/bottom N disease
     genes (same thresholds, averaged).
  5. "Extreme" Spearman correlation over the same restricted gene set
     (same thresholds, averaged).

All metrics are combined into a per-compound "combined_rank" (average of the
per-metric ranks within each trait), where a lower combined_rank indicates
a stronger, more consistent reversal signal (compound signature opposes the
disease signature).

Processes compounds in batches to bound peak memory usage regardless of the
number of compounds or traits being analysed; only ever holds one batch's
worth of the compound matrix in memory at a time.

Usage:
    python 02_connectivity_analysis.py \
        --disease_signatures disease_signatures_for_cmap.tsv.gz \
        --cmap_matrix cmap/level5_filtered_matrix.parquet \
        --output cmap_drug_repurposing_results.tsv.gz \
        --batch_size 1500 --thresholds 50 100 250 500
"""
import argparse
import gc

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy import stats


def column_rank(mat):
    """Independent ascending rank (1 = smallest value) within each column, vectorised."""
    order = np.argsort(mat, axis=0).astype(np.int32)
    ranks = np.empty_like(order, dtype=np.float32)
    rows = np.arange(mat.shape[0], dtype=np.float32)
    for j in range(mat.shape[1]):
        ranks[order[:, j], j] = rows + 1
    return ranks


def pearson_vs_all(x, y):
    """x: (genes,), y: (genes, compounds) -> (compounds,) Pearson r of x against every column of y."""
    xc = x - np.nanmean(x)
    yc = y - y.mean(axis=0, keepdims=True)
    num = (xc[:, None] * yc).sum(axis=0)
    den = np.sqrt((xc ** 2).sum()) * np.sqrt((yc ** 2).sum(axis=0))
    return num / den


def ks_signed_one_threshold(disease_vec, rank_desc_mat, thres_n):
    """Classical CMap signed KS enrichment score at one gene-set size threshold."""
    nogenes = len(disease_vec)
    order_desc = np.argsort(-disease_vec)
    ind_up = order_desc[:thres_n]
    ind_down = order_desc[::-1][:thres_n]

    def one_side(indices):
        geneset2 = np.sort(rank_desc_mat[indices, :], axis=0)
        j = np.arange(1, thres_n + 1).reshape(-1, 1)
        a = np.max(j / thres_n - geneset2 / nogenes, axis=0)
        b = np.max(geneset2 / nogenes - (j - 1) / thres_n, axis=0)
        return np.where(a > b, a, -b)

    es_up = one_side(ind_up)
    es_down = one_side(ind_down)
    return np.where(np.sign(es_up) != np.sign(es_down), es_up - es_down, 0.0)


def extreme_corr_one_threshold(dz_v, drug_mat, disease_order_desc, thres_n, method):
    ind_up = disease_order_desc[:thres_n]
    ind_down = disease_order_desc[::-1][:thres_n]
    idx = np.concatenate([ind_up, ind_down])
    dz_sub = dz_v[idx]
    if method == "pearson":
        return pearson_vs_all(dz_sub, drug_mat[idx, :])
    sub = drug_mat[idx, :]
    sub_rank = column_rank(sub)
    dz_rank = stats.rankdata(dz_sub).astype(np.float32)
    return pearson_vs_all(dz_rank, sub_rank)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--disease_signatures", required=True,
                     help="TSV.gz with genes as rows, traits as columns (Z-scores)")
    ap.add_argument("--cmap_matrix", required=True,
                     help="Parquet file with genes as index, compounds as columns (Z-scores)")
    ap.add_argument("--output", required=True)
    ap.add_argument("--batch_size", type=int, default=1500)
    ap.add_argument("--thresholds", type=int, nargs="+", default=[50, 100, 250, 500])
    args = ap.parse_args()

    print("Loading disease signatures...", flush=True)
    disease = pd.read_csv(args.disease_signatures, sep="\t", index_col=0)

    print("Reading CMap matrix gene order and compound column names (metadata only)...", flush=True)
    pf = pq.ParquetFile(args.cmap_matrix)
    all_cols = list(pf.schema_arrow.names)
    gene_col = "gene_name"
    drug_cols = [c for c in all_cols if c != gene_col]
    n_drugs_total = len(drug_cols)
    print(f"Total compounds: {n_drugs_total}", flush=True)

    gene_order = pd.read_parquet(args.cmap_matrix, columns=[]).index.to_numpy()
    n_genes = len(gene_order)
    print(f"Genes: {n_genes}", flush=True)

    disease_aligned = disease.reindex(gene_order)
    traits = disease_aligned.columns.tolist()

    trait_prep = {}
    for trait in traits:
        dz = disease_aligned[trait].to_numpy(dtype=np.float32)
        valid = ~np.isnan(dz)
        dz_v = dz[valid]
        order_desc = np.argsort(-dz_v)
        dz_rank_full = stats.rankdata(dz_v).astype(np.float32)
        trait_prep[trait] = dict(valid=valid, dz_v=dz_v, order_desc=order_desc, dz_rank=dz_rank_full,
                                  n_valid=int(valid.sum()))
        print(f"  {trait}: {valid.sum()}/{n_genes} valid genes", flush=True)

    all_results = {trait: [] for trait in traits}
    n_batches = (n_drugs_total + args.batch_size - 1) // args.batch_size

    for b in range(n_batches):
        start, end = b * args.batch_size, min((b + 1) * args.batch_size, n_drugs_total)
        cols = drug_cols[start:end]
        print(f"\nBatch {b + 1}/{n_batches}: compounds {start}-{end}", flush=True)
        chunk_df = pd.read_parquet(args.cmap_matrix, columns=cols)
        chunk_vals = chunk_df.to_numpy(dtype=np.float32)
        del chunk_df
        chunk_rank = column_rank(chunk_vals)

        for trait in traits:
            tp = trait_prep[trait]
            valid, dz_v, n_valid = tp["valid"], tp["dz_v"], tp["n_valid"]

            cv = chunk_vals[valid, :]
            cr = column_rank(cv) if n_valid < n_genes else chunk_rank
            cr_desc = n_valid + 1 - cr

            pear_all = pearson_vs_all(dz_v, cv)
            spear_all = pearson_vs_all(tp["dz_rank"], cr)

            ks_list, ext_p_list, ext_s_list = [], [], []
            for n in args.thresholds:
                ks_list.append(ks_signed_one_threshold(dz_v, cr_desc, n))
                ext_p_list.append(extreme_corr_one_threshold(dz_v, cv, tp["order_desc"], n, "pearson"))
                ext_s_list.append(extreme_corr_one_threshold(dz_v, cv, tp["order_desc"], n, "spearman"))

            all_results[trait].append(pd.DataFrame({
                "compound": cols,
                "pearson_all": pear_all,
                "spearman_all": spear_all,
                "ks_signed_mean": np.mean(ks_list, axis=0),
                "extreme_pearson_mean": np.mean(ext_p_list, axis=0),
                "extreme_spearman_mean": np.mean(ext_s_list, axis=0),
            }))
        del chunk_vals, chunk_rank
        gc.collect()

    print("\nAggregating all batches...", flush=True)
    final_frames = []
    metric_cols = ["pearson_all", "spearman_all", "ks_signed_mean", "extreme_pearson_mean", "extreme_spearman_mean"]
    for trait in traits:
        df = pd.concat(all_results[trait], ignore_index=True)
        df["trait"] = trait
        for col in metric_cols:
            df[col + "_rank"] = df[col].rank(method="average")
        df["combined_rank"] = df[[c for c in df.columns if c.endswith("_rank")]].mean(axis=1)
        df.sort_values("combined_rank", inplace=True)
        final_frames.append(df)
        print(f"{trait} top 5 reversal candidates: {df.head(5).compound.tolist()}", flush=True)

    allres = pd.concat(final_frames, ignore_index=True)
    allres.to_csv(args.output, sep="\t", index=False, compression="gzip")
    print(f"\nDone. Results saved to {args.output} ({len(allres):,} rows)", flush=True)


if __name__ == "__main__":
    main()
