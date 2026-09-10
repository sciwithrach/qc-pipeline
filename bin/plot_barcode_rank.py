#!/usr/bin/env python3
import argparse
import math
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def load_counts_h5ad(h5ad_path, batch_id):
    adata = ad.read_h5ad(h5ad_path)
    if ("tscp_count" in adata.obs.columns) & (batch_id in adata.obs.columns):
        return pd.DataFrame({
            'counts': adata.obs["tscp_count"].values.astype(float),
            'plate': adata.obs[batch_id]
        })
    return np.asarray(adata.X.sum(axis=1)).flatten()


def compute_barcode_ranks(counts, batch_id):
    if isinstance(counts, pd.DataFrame):
        df = counts
        barcode_rank_dict = {}
        for batch in counts[batch_id].cat.categories:
            counts = df[df[batch_id] == batch]
            sorted_counts = np.sort(counts['counts'])[::-1]
            ranks = np.arange(len(sorted_counts), dtype=float)
            barcode_rank_dict[batch] = (ranks, sorted_counts)
        return barcode_rank_dict
    sorted_counts = np.sort(counts)[::-1]
    ranks = np.arange(len(sorted_counts), dtype=float)
    return (ranks, sorted_counts)


def plot_barcode_rank(barcode_rank_dict, batch_id, title=None):
    fig, ax = plt.subplots(figsize=(8, 5))

    colors = [
        '#5487ff', '#E36C2A', '#ED367D', '#57b444',
        '#002ea6', '#781c1c', '#9e1a54', '#005e5c'
    ]

    batches = list(barcode_rank_dict.keys())
    n_batches = len(batches)
    if n_batches > 8:
        factor = math.ceil(n_batches / 8)
        colors = colors * factor

    color_dict = {batch: colors[i] for i, batch in enumerate(batches)}

    for batch, (ranks, counts) in barcode_rank_dict.items():
        ax.plot(ranks, counts, color=color_dict[batch], linewidth=1.2, label=batch)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Barcode rank", fontsize=12)
    ax.set_ylabel("Transcript count", fontsize=12)
    ax.set_title(title or f"Barcode rank plot by {batch_id}", fontsize=13)
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, which="both", alpha=0.3, linestyle=":")

    plt.tight_layout()
    fig.savefig(f'plot_barcode_rank_{batch_id}.png', dpi=150, bbox_inches="tight")
    plt.close(fig)


def build_parser():
    p = argparse.ArgumentParser()
    p.add_argument('--adata', required=True)
    p.add_argument('--batch', required=True)
    return p


def main():
    args = build_parser().parse_args()

    counts = load_counts_h5ad(args.adata, args.batch)
    barcode_rank_dict = compute_barcode_ranks(counts, args.batch)
    plot_barcode_rank(barcode_rank_dict, batch_id=args.batch)


if __name__ == '__main__':
    main()
