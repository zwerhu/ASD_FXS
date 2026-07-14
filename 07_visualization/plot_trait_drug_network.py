#!/usr/bin/env python3
"""
07_visualization/plot_trait_drug_network.py

Network showing, for each trait, its top N significant reversal candidate
compounds by empirical Z-score, plus any compounds significant in >= K
traits (highlighted, since these represent the broadest cross-trait
consistency).

Usage:
    python plot_trait_drug_network.py \
        --significant_drugs cmap_significant_drugs.tsv \
        --output figures/trait_drug_network.png \
        --top_per_trait 8 --shared_min_traits 3
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
    ap.add_argument("--significant_drugs", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--top_per_trait", type=int, default=8)
    ap.add_argument("--shared_min_traits", type=int, default=3)
    args = ap.parse_args()

    sig = pd.read_csv(args.significant_drugs, sep="\t")

    top_per_trait = sig.sort_values("ks_z").groupby("trait").head(args.top_per_trait)
    shared = sig.groupby("compound").filter(lambda g: g.trait.nunique() >= args.shared_min_traits)
    plot_data = pd.concat([top_per_trait, shared]).drop_duplicates(subset=["trait", "compound"])
    print(f"Network size: {plot_data.trait.nunique()} traits, {plot_data.compound.nunique()} compounds, "
          f"{len(plot_data)} edges", flush=True)

    G = nx.Graph()
    for _, row in plot_data.iterrows():
        G.add_edge(row.trait, row.compound, ks_z=row.ks_z)

    trait_nodes = set(plot_data.trait.unique())
    compound_nodes = set(plot_data.compound.unique())
    shared_compounds = set(shared.compound.unique())
    regular_compounds = compound_nodes - shared_compounds

    pos = nx.spring_layout(G, seed=11, k=1.6, iterations=800)
    fig, ax = plt.subplots(figsize=(22, 19))

    edges = list(G.edges(data=True))
    edge_colors = [d["ks_z"] for _, _, d in edges]
    vmax = max(abs(min(edge_colors)), abs(max(edge_colors)))
    nx.draw_networkx_edges(G, pos, edgelist=[(u, v) for u, v, d in edges], edge_color=edge_colors,
                            edge_cmap=plt.cm.RdBu_r, edge_vmin=-vmax, edge_vmax=vmax, width=2.2, ax=ax, alpha=0.85)

    nx.draw_networkx_nodes(G, pos, nodelist=list(trait_nodes), node_color="#ffb703",
                            node_size=4200, node_shape="s", edgecolors="black", linewidths=1.5, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=list(regular_compounds), node_color="#adb5bd",
                            node_size=750, node_shape="o", edgecolors="black", ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=list(shared_compounds), node_color="#e63946",
                            node_size=1150, node_shape="o", edgecolors="black", linewidths=2.2, ax=ax)

    for n in trait_nodes:
        x, y = pos[n]
        ax.text(x, y, n, fontsize=15, fontweight="bold", ha="center", va="center", zorder=10)
    for n in shared_compounds:
        x, y = pos[n]
        ax.text(x, y, n, fontsize=11, ha="center", va="center", zorder=10,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75))

    sm = plt.cm.ScalarMappable(cmap=plt.cm.RdBu_r, norm=plt.Normalize(vmin=-vmax, vmax=vmax))
    cbar = plt.colorbar(sm, ax=ax, shrink=0.6, label="ks_z (more negative = stronger reversal signal)")
    cbar.ax.tick_params(labelsize=13)

    legend_elems = [
        Line2D([0], [0], marker='s', color='w', markerfacecolor='#ffb703', markeredgecolor='black', markersize=20, label='Trait (GWAS)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#adb5bd', markeredgecolor='black', markersize=13, label='Significant compound (single trait)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#e63946', markeredgecolor='black', markersize=15, label=f'Shared across >={args.shared_min_traits} traits'),
    ]
    ax.legend(handles=legend_elems, loc="upper left", fontsize=15, framealpha=0.9)
    ax.set_title(f"Trait - Significant Drug Candidate Network\n"
                 f"(top {args.top_per_trait} by ks_z per trait + compounds shared across "
                 f">={args.shared_min_traits} traits)", fontsize=19)
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
