#!/usr/bin/env Rscript

library(clustree)
library(anndataR)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(stringr)

args         <- commandArgs(trailingOnly = TRUE)
adata_path   <- args[1]
celltype_csv <- args[2]

adata <- read_h5ad(adata_path)
sce   <- adata$as_SingleCellExperiment()

ctype <- read.csv(celltype_csv, header = TRUE, row.names = 1)
sce$predicted_celltype         <- ctype[rownames(ctype) %in% colnames(sce), 'predicted_celltype']
sce$predicted_celltype_summary <- str_remove(sce$predicted_celltype, '_.*')

if (grepl('ms', adata$uns[['basis']])) {
    basis <- 'topo_clusters_ms_res'
    title <- 'Topometry - multiscale latent space'
    p <- clustree(sce, prefix = basis) + ggtitle(title)
} else if (grepl('topo', adata$uns[['basis']])) {
    basis <- 'topo_clusters_res'
    title <- 'Topometry - spectral latent space'
    p <- clustree(sce, prefix = basis) + ggtitle(title)
} else {
    basis <- 'pca_leiden_res'
    title <- 'PCA-based latent space'
    p <- clustree(sce, prefix = basis) + ggtitle(title)
}

ggsave("clustree.png", p, width = 10, height = 15)
