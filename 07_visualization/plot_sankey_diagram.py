#!/usr/bin/env python3
"""
07_visualization/plot_sankey_diagram.py

Sankey diagram: Trait -> Drug (with MOA) -> Driving Gene -> GTEx Tissue,
built from the output of
06_drug_repurposing/04_drug_disease_gene_tissue_links.py.

Requires plotly and kaleido for static PNG export
(pip install plotly "kaleido==0.2.1" -- newer kaleido versions require a
separately-installed Chrome browser to render, which older 0.2.x releases
bundle internally).

Usage:
    python plot_sankey_diagram.py \
        --links cmap_drug_disease_gene_tissue.tsv \
        --output figures/sankey_diagram.png \
        --title "Top drug per trait, top 3 driving genes each, linked to GTEx tissue"
"""
import argparse

import pandas as pd
import plotly.graph_objects as go

DEFAULT_TRAIT_COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd", "#8c564b", "#7f7f7f", "#bcbd22", "#17becf",
    "#e377c2", "#aec7e8", "#ffbb78", "#98df8a",
]


def hex_to_rgba(hex_color, alpha):
    hex_color = hex_color.lstrip("#")
    r, g, b = tuple(int(hex_color[i:i + 2], 16) for i in (0, 2, 4))
    return f"rgba({r},{g},{b},{alpha})"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--links", required=True, help="Output of 04_drug_disease_gene_tissue_links.py")
    ap.add_argument("--output", required=True)
    ap.add_argument("--title", default="Drug - Disease - Gene - Tissue Sankey Diagram")
    ap.add_argument("--width", type=int, default=2500)
    ap.add_argument("--height", type=int, default=1500)
    ap.add_argument("--font_size", type=int, default=27)
    ap.add_argument("--title_font_size", type=int, default=34)
    args = ap.parse_args()

    df = pd.read_csv(args.links, sep="\t").reset_index(drop=True)

    traits = df.trait.unique().tolist()
    drugs = df.drop_duplicates(subset="compound").compound.tolist()
    drug_labels = [f"{d} ({df[df.compound == d].MOA.iloc[0]})" for d in drugs]
    n_rows = len(df)
    tissues = sorted(df.best_tissue.unique().tolist())

    node_labels = traits + drug_labels + df.gene.tolist() + tissues
    n_traits, n_drugs = len(traits), len(drugs)
    trait_idx = {t: i for i, t in enumerate(traits)}
    drug_idx = {d: n_traits + i for i, d in enumerate(drugs)}
    gene_idx = {i: n_traits + n_drugs + i for i in range(n_rows)}
    tissue_idx = {t: n_traits + n_drugs + n_rows + i for i, t in enumerate(tissues)}

    trait_colors = {t: DEFAULT_TRAIT_COLORS[i % len(DEFAULT_TRAIT_COLORS)] for i, t in enumerate(traits)}

    node_colors = (
        [trait_colors[t] for t in traits] +
        [trait_colors[df[df.compound == d].trait.iloc[0]] for d in drugs] +
        [trait_colors[df.iloc[i].trait] for i in range(n_rows)] +
        ["#e63946" if t.startswith("Brain") else "#6c757d" for t in tissues]
    )

    source, target, value, link_color = [], [], [], []
    for i, row in df.iterrows():
        c = trait_colors[row.trait]
        w = max(row.contribution, 0.1)
        source.append(trait_idx[row.trait]); target.append(drug_idx[row.compound]); value.append(w)
        link_color.append(hex_to_rgba(c, 0.35))
        source.append(drug_idx[row.compound]); target.append(gene_idx[i]); value.append(w)
        link_color.append(hex_to_rgba(c, 0.35))
        source.append(gene_idx[i]); target.append(tissue_idx[row.best_tissue]); value.append(w)
        link_color.append(hex_to_rgba(c, 0.35))

    fig = go.Figure(data=[go.Sankey(
        arrangement="snap",
        node=dict(
            pad=22, thickness=26,
            line=dict(color="black", width=0.8),
            label=node_labels,
            color=node_colors,
            x=[0.001] * n_traits + [0.33] * n_drugs + [0.67] * n_rows + [0.999] * len(tissues),
        ),
        link=dict(source=source, target=target, value=value, color=link_color),
        textfont=dict(size=args.font_size, color="#000000", family="Arial Black")
    )])

    fig.update_layout(
        title=dict(text=args.title, font=dict(size=args.title_font_size, color="#000000")),
        font=dict(size=args.font_size - 3, color="#000000"),
        width=args.width, height=args.height,
        margin=dict(l=20, r=20, t=110, b=20),
        paper_bgcolor="white"
    )

    fig.write_image(args.output, scale=2)
    html_output = args.output.rsplit(".", 1)[0] + ".html"
    fig.write_html(html_output)
    print(f"Saved {args.output} and {html_output}")


if __name__ == "__main__":
    main()
