#!/bin/bash
#SBATCH --job-name=QC  # Job name
#SBATCH --mail-type=ALL           # notifications for job done & fail
#SBATCH --mail-user=eknodel@asu.edu # send-to address
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 2-00:00:00
#SBATCH --mem=40000
#SBATCH -q public
#SBATCH -p general

source activate cancergenomics

snakemake --snakefile quality_control.snakefile -j 50 --keep-target-files --rerun-incomplete --cluster "sbatch -n 1 -c 1 -q public -p general --mem=50000 -t 0-10:00:00"
