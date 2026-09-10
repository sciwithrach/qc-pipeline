#!/usr/bin/env Rscript

set.seed(42)

library(scDblFinder)
library(anndataR)
library(BiocParallel)

args       <- commandArgs(trailingOnly = TRUE)
adata_path <- args[1]
batch      <- args[2]
cpus       <- as.integer(args[3])

bp <- MulticoreParam(cpus, RNGseed = 42)

adata <- read_h5ad(adata_path)
sce   <- adata$as_SingleCellExperiment()
counts(sce) <- assay(sce, 'counts')

sce <- scDblFinder(
    sce,
    samples   = batch,
    nfeatures = sum(rowData(sce)$highly_variable),
    BPPARAM   = bp
)

results <- data.frame(
    scDblFinder_calls  = sce$scDblFinder.class,
    scDblFinder_scores = sce$scDblFinder.score,
    row.names          = colnames(sce)
)

write.csv(results, 'doublets_scdblfinder.csv')
writeLines(capture.output(sessionInfo()), "sessionInfo_doublets.txt")
