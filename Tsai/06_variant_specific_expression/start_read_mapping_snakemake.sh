#!/bin/bash
#SBATCH --job-name=HISAT2_mapping  # Job name
#SBATCH --mail-type=ALL           # notifications for job done & fail
#SBATCH --mail-user=eknodel@asu.edu # send-to address
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH --mem=1024

newgrp combinedlab

#source activate cancergenomics
source activate read_mapping_environment
PERL5LIB=/packages/6x/vcftools/0.1.12b/lib/per15/site_perl

snakemake --snakefile hisat2_read_mapping_ensembl.snakefile -j 15 --keep-target-files --rerun-incomplete --cluster "sbatch -n 1 -c 8 -t 24:00:00"
