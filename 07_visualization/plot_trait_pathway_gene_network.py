#!/usr/bin/env python3
"""
07_visualization/plot_trait_pathway_gene_network.py

Network diagram connecting traits to their significantly enriched Hallmark
pathways (FDR < threshold), and each pathway to its top leading-edge genes.

Usage:
    python plot_trait_pathway_gene_network.py \
        --gsea_results gsea_ALL_hallmark_results.tsv.gz \
        --output figures/trait_pathway_gene_network.png \
        --fdr_threshold 0.25 --max_genes_per_pathway 8
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from matplotlib.lines import Line2D


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gsea_results", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--fdr_threshold", type=float, default=0.25)
    ap.add_argument("--max_genes_per_pathway", type=int, default=8)
    args = ap.parse_args()

    gsea = pd.read_csv(args.gsea_results, sep="\t")
    gsea["FDR q-val"] = gsea["FDR q-val"].astype(float)
    sig = gsea[gsea["FDR q-val"] < args.fdr_threshold].copy()
    sig["Term_short"] = sig.Term.str.replace("HALLMARK_", "", regex=False)

    G = nx.Graph()
    trait_nodes, pathway_nodes, gene_nodes = set(), set(), set()

    for _, row in sig.iterrows():
        trait, pathway = row.trait, row.Term_short
        genes = row.Lead_genes.split(";")[:args.max_genes_per_pathway]
        trait_nodes.add(trait)
        pathway_nodes.add(pathway)
        G.add_edge(trait, pathway, kind="trait_pathway", nes=row.NES)
        for g in genes:
            gene_nodes.add(g)
            G.add_edge(pathway, g, kind="pathway_gene")

    pos = nx.spring_layout(G, seed=7, k=0.9, iterations=200)
    fig, ax = plt.subplots(figsize=(14, 11))

    nx.draw_networkx_edges(G, pos, edgelist=[(u, v) for u, v, d in G.edges(data=True) if d["kind"] == "pathway_gene"],
                            edge_color="#cccccc", width=1.0, ax=ax)
    trait_pathway_edges = [(u, v) for u, v, d in G.edges(data=True) if d["kind"] == "trait_pathway"]
    edge_colors = ["#d62728" if G[u][v]["nes"] < 0 else "#1f77b4" for u, v in trait_pathway_edges]
    nx.draw_networkx_edges(G, pos, edgelist=trait_pathway_edges, edge_color=edge_colors, width=3, ax=ax)

    nx.draw_networkx_nodes(G, pos, nodelist=list(trait_nodes), node_color="#ffb703",
                            node_size=2400, node_shape="s", edgecolors="black", ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=list(pathway_nodes), node_color="#8ecae6",
                            node_size=1800, node_shape="D", edgecolors="black", ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=list(gene_nodes), node_color="#adb5bd",
                            node_size=500, node_shape="o", edgecolors="black", ax=ax)

    nx.draw_networkx_labels(G, pos, {n: n for n in trait_nodes}, font_size=8, font_weight="bold", ax=ax)
    nx.draw_networkx_labels(G, pos, {n: n for n in pathway_nodes}, font_size=7.5, ax=ax)
    nx.draw_networkx_labels(G, pos, {n: n for n in gene_nodes}, font_size=6.5, ax=ax)

    legend_elems = [
        Line2D([0], [0], marker='s', color='w', markerfacecolor='#ffb703', markeredgecolor='black', markersize=14, label='Trait (GWAS)'),
        Line2D([0], [0], marker='D', color='w', markerfacecolor='#8ecae6', markeredgecolor='black', markersize=12, label=f'Pathway (FDR<{args.fdr_threshold})'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#adb5bd', markeredgecolor='black', markersize=10, label='Leading-edge gene'),
        Line2D([0], [0], color='#1f77b4', lw=3, label='Trait-pathway (NES > 0, enriched up)'),
        Line2D([0], [0], color='#d62728', lw=3, label='Trait-pathway (NES < 0, enriched down)'),
    ]
    ax.legend(handles=legend_elems, loc="upper left", fontsize=9, framealpha=0.9)
    ax.set_title(f"Trait - Pathway - Gene network\n(GSEA Hallmark pathways at FDR<{args.fdr_threshold}, "
                 f"top {args.max_genes_per_pathway} leading-edge genes each)", fontsize=12)
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    print(f"Saved {args.output}")
    print(f"Nodes: {G.number_of_nodes()} ({len(trait_nodes)} traits, {len(pathway_nodes)} pathways, {len(gene_nodes)} genes)")


if __name__ == "__main__":
    main()
