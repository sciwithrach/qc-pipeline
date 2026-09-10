#!/usr/bin/env

# cellbender
apptainer pull containers/cellbender.sif docker://us.gcr.io/broad-dsde-methods/cellbender:0.3.2

# utricle qc
apptainer build containers/utricle-qc.sif containers/utricle-qc.def

# vaeda
apptainer build containers/vaeda.sif containers/vaeda.def

# r-clustering
apptainer build containers/r-clustering.sif containers/r-clustering.def