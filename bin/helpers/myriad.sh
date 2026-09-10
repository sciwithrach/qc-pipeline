# Myriad script to run QC pipeline for utricle using GPU acceleration.

# Requires Nextflow install as per UCL guidance:
# https://github.com/nf-core/configs/blob/master/docs/ucl_myriad.md

#!/bin/bash -l

# Request 1 GPU
#$ -l gpu=1

# Request wallclock time (h:m:s)
#$ -l h_rt=48:0:0

# Request RAM (integer M/G/T e.g. 1G)
#$ -l mem=4G

# Request TMPDIR space (default is 10G)
#$ -l tmpfs=15G

# Request cores
#$ -pe smp 36

# Set job name
#$ -N UtricleQC

# Set the working directory to somewhere in scratch
#$ -wd /home/sjjgrww/Scratch/output

# Change into temporary working directory to run script
cd $TMPDIR

# load the required modules - GPU
module unload compiles mpi
module load compilers/gnu/10.2.0
module load cuda/12.2.2/gnu-10.2.0

# load the required modules - QC pipeline
module load apptainer
module load java/temurin-17/17.0.2_8

### WHAT IS THE CORRECT WD? ###

# download and store the necessary containers for the pipeline
bash setup-containers.sh

# run the pipeline itself
nextflow run qc-pipeline.nf -profile ucl_myriad --outdir $TMPDIR --workDir

###############################

# save all outputs
tar zcvf $HOME/Scratch/qc-pipeline-$JOB_ID.tar.gz $TMPDIR

# allow enough time for the copy to complete