#!/usr/bin/env python3
import argparse
import anndata as ad


def build_parser():
    p = argparse.ArgumentParser()
    p.add_argument('--adata', required=True)
    return p


def main():
    args = build_parser().parse_args()

    adata = ad.read_h5ad(args.adata)

    # capture the HVG gene names before switching to raw -- raw is genuinely
    # unnormalised, all-gene counts, so it has no highly_variable annotation
    hvg_genes = adata.var_names[adata.var.highly_variable].tolist()

    adata = adata.raw.to_adata()
    adata.layers['counts'] = adata.X
    # run_scDblFinder.R reads rowData(sce)$highly_variable directly, so this
    # needs to be a real column on the all-gene object, not just .uns
    adata.var['highly_variable'] = adata.var_names.isin(hvg_genes)
    adata.uns['hvg_list'] = hvg_genes
    adata.write_h5ad('adata_scDblFinder_in.h5ad')

    adata.uns = {}
    adata = adata[:, adata.var_names.isin(hvg_genes)]
    adata.write_h5ad('adata_vaeda_in.h5ad')


if __name__ == '__main__':
    main()
