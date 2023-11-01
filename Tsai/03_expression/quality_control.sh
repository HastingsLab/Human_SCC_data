#!/bin/bash
#SBATCH --job-name=QC  # Job name
#SBATCH --mail-type=ALL           # notifications for job done & fail
#SBATCH --mail-user=eknodel@asu.edu # send-to address
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 48:00:00
#SBATCH --mem=40000
#SBATCH -q tempboost

newgrp combinedlab
source activate salmon_environment
#source activate cancergenomics

snakemake --snakefile quality_control_salmon.snakefile -j 30 --keep-target-files --rerun-incomplete --cluster "sbatch -n 1 -c 1 -q tempboost --mem=50000 -t 48:00:00"
