#!/bin/bash
#SBATCH --job-name=Make_config  # Job name
#SBATCH --mail-type=ALL           # notifications for job done & fail
#SBATCH --mail-user=eknodel@asu.edu # send-to address
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 0-10:00:00
#SBATCH --mem=40000
#SBATCH -q public
#SBATCH -p general

python generate_config.py --fastq_path /data/CEM/shared/neoantigens/fastq/ \
                          --sample_info /home/eknodel/Zheng_CancerGenomics/00_misc/Zheng_metadata_for_pipeline.csv \
                          --ref_dir /data/CEM/shared/public_data/references/1000genomes_GRCh38_reference_genome \
                          --ref_basename GRCh38_full_analysis_set_plus_decoy_hla \
                          --varscan_path /home/eknodel/programs/VarScan.v2.3.9.jar \
                          --gatk_path /home/eknodel/programs/gatk-4.1.7.0/gatk \
                          --strelka /home/eknodel/programs/strelka-2.9.10-0/bin/configureStrelkaSomaticWorkflow.py \
                          --bam_readcount bam-readcount \
                          --perl_fp_filter fpfilter-2.pl
