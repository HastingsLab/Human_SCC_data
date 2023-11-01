#!/bin/bash
#SBATCH --job-name=NeoScore  # Job name
#SBATCH --mail-type=ALL           # notifications for job done & fail
#SBATCH --mail-user=eknodel@asu.edu # send-to address
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 12:00:00
#SBATCH --mem=40000
#SBATCH -q tempboost

newgrp combinedlab

#source activate cancergenomics
source activate var_call_env


PERL5LIB=/packages/6x/vcftools/0.1.12b/lib/per15/site_perl

snakemake --snakefile  Neoscore.snakefile -j 200 --keep-target-files --rerun-incomplete --cluster "sbatch -n 1 -c 1 --mem=50000 -t 12:00:00 -q tempboost"
