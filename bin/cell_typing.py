#!/usr/bin/env python3
import argparse
import anndata as ad
import pandas as pd
import numpy as np


def csv_to_dict(csv):
    df = pd.read_csv(csv, header=None)
    return {row[0]: [gene for gene in row[1:] if pd.notna(gene)] for _, row in df.iterrows()}


def simple_classifier(adata, marker_dict, key_added='celltype'):
    scores = {}
    for celltype, genes in marker_dict.items():
        valid_genes = [g for g in genes if g in adata.var_names]
        if valid_genes:
            scores[celltype] = np.array(adata[:, valid_genes].X.mean(axis=1)).flatten()
    if not scores:
        raise ValueError("None of the marker genes are present in the dataset.")
    scores_matrix = np.vstack(list(scores.values())).T
    predicted_labels = np.array(list(scores.keys()))[np.argmax(scores_matrix, axis=1)]
    adata.obs[key_added] = pd.Categorical(predicted_labels)


def build_parser():
    p = argparse.ArgumentParser()
    p.add_argument('--adata', required=True)
    p.add_argument('--type-csv', required=True)
    p.add_argument('--descriptor', required=True)
    return p


def main():
    args = build_parser().parse_args()

    adata = ad.read_h5ad(args.adata)
    marker_dict = csv_to_dict(args.type_csv)
    simple_classifier(adata, marker_dict, key_added='predicted_celltype')
    adata.obs.to_csv(f'simple_cell_types_{args.descriptor}.csv')


if __name__ == '__main__':
    main()
