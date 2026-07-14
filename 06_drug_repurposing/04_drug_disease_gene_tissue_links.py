#!/usr/bin/env python3
"""
06_drug_repurposing/04_drug_disease_gene_tissue_links.py

For a selected set of significant, MOA-annotated compound-trait pairs
(typically the top reversal candidate per trait), extracts the genes that
contribute most to the reversal signal (largest disease_Z x compound_Z
product with opposite sign) and links each gene to the GTEx tissue in which
its imputed expression association was strongest (from the cross-tissue
integrative S-PrediXcan results), producing a table suitable for network or
Sankey visualisation (see 07_visualization/plot_sankey_diagram.py).

Usage:
    python 04_drug_disease_gene_tissue_links.py \
        --significant_drugs cmap_significant_drugs.tsv \
        --moa_table broad_moa.txt \
        --disease_signatures disease_signatures_for_cmap.tsv.gz \
        --cmap_matrix cmap/level5_filtered_matrix.parquet \
        --integrative_results integrative_ACAT_results.tsv.gz \
        --top_drugs_per_trait 1 --top_genes_per_drug 3 \
        --output cmap_drug_disease_gene_tissue.tsv
"""
import argparse

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--significant_drugs", required=True)
    ap.add_argument("--moa_table", required=True)
    ap.add_argument("--disease_signatures", required=True)
    ap.add_argument("--cmap_matrix", required=True)
    ap.add_argument("--integrative_results", required=True,
                     help="Cross-tissue ACAT results, used to look up each gene's best-supporting tissue")
    ap.add_argument("--top_drugs_per_trait", type=int, default=1)
    ap.add_argument("--top_genes_per_drug", type=int, default=3)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    sig = pd.read_csv(args.significant_drugs, sep="\t")
    moa = pd.read_csv(args.moa_table, sep="\t", low_memory=False)
    moa_map = moa.drop_duplicates(subset="CompoundName")[["CompoundName", "MOA"]]
    moa_map["compound_lower"] = moa_map.CompoundName.str.lower()
    sig["compound_lower"] = sig.compound.str.lower()
    sig_named = sig.merge(moa_map[["compound_lower", "MOA"]], on="compound_lower", how="inner")
    print(f"Significant candidates with resolvable MOA: {len(sig_named)}", flush=True)

    selected = sig_named.sort_values("ks_z").groupby("trait").head(args.top_drugs_per_trait)
    print(selected[["trait", "compound", "MOA", "ks_z"]].to_string(index=False))

    cmap = pd.read_parquet(args.cmap_matrix)
    disease = pd.read_csv(args.disease_signatures, sep="\t", index_col=0)
    disease = disease.reindex(cmap.index)

    integ = pd.read_csv(args.integrative_results, sep="\t")
    gene_tissue = integ.set_index(["trait", "gene_name"])["best_tissue"]

    rows = []
    for _, r in selected.iterrows():
        trait, compound = r.trait, r.compound
        if compound not in cmap.columns:
            continue
        dz = disease[trait]
        dg = cmap[compound]
        valid = dz.notna() & dg.notna()
        dzv, dgv = dz[valid], dg[valid]
        # Contribution: more negative dz*dg (opposite direction, large magnitude) -> larger
        # positive "contribution" score after sign flip, i.e. bigger reversal contribution.
        contrib = -(dzv * dgv)
        top_genes = contrib.sort_values(ascending=False).head(args.top_genes_per_drug)
        for gene, c in top_genes.items():
            tissue = gene_tissue.get((trait, gene), "NA")
            rows.append(dict(trait=trait, compound=compound, MOA=r.MOA, gene=gene,
                              disease_z=dzv[gene], drug_z=dgv[gene], contribution=c, best_tissue=tissue))

    net_df = pd.DataFrame(rows)
    net_df.to_csv(args.output, sep="\t", index=False)
    print(f"\n{len(net_df)} drug-disease-gene-tissue links written to {args.output}")


if __name__ == "__main__":
    main()
