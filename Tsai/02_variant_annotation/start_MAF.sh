#!/bin/bash
#SBATCH --job-name=Gatk  # Job name
#SBATCH --mail-type=ALL           # notifications for job done & fail
#SBATCH --mail-user=eknodel@asu.edu # send-to address
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 0-00:30:00
#SBATCH --mem=40000
#SBATCH -p general
#SBATCH -q public 

#source activate var_call_env
#source activate cancergenomics_new_pvacseq
source activate vep_env

snakemake --snakefile MAF.snakefile -j 71 --keep-target-files --rerun-incomplete --cluster "sbatch -p general -q public -n 1 -c 1 --mem=50000 -t 0-00:30:00"
