#!/bin/bash
#SBATCH --job-name=footprint  # Job name
#SBATCH --mail-type=ALL           # notifications for job done & fail
#SBATCH --mail-user=eknodel@asu.edu # send-to address
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 96:00:00
#SBATCH --mem=40000
#SBATCH -q tempboost

newgrp combinedlab
source activate cancergenomics

module load bedtools2/2.24.0

snakemake --snakefile footprint.snakefile -j 30 --keep-target-files --rerun-incomplete --cluster "sbatch -n 1 -c 1 -q tempboost --mem=50000 -t 96:00:00"
