#!/bin/bash
#SBATCH --job-name=Gatk-cnv  # Job name
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 168:00:00
#SBATCH --mem=40000
#SBATCH -p wildfire
#SBATCH -q wildfire

newgrp combinedlab

source activate CNV_env_2
#source activate cancergenomics
#module load picard/latest

snakemake --snakefile GATK_CNV.snakefile -j 190 --keep-target-files --rerun-incomplete --cluster "sbatch -p wildfire -q wildfire -n 1 -c 1 --mem=50000 -t 96:00:00"
