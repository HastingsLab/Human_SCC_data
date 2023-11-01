# Setting up filesnames here:
from os.path import join
import os

configfile: "somatic_mutation_calling_config.json"

rule all:
    input:
        expand("variants/{subject}_merged_vep.maf", subject=config["all_subjects"]), # combine all mutations


rule maf:
	input:
		"variants/{subject}_merged_vep.vcf"
	output:
		"variants/{subject}_merged_vep.maf"
	params:
		sample = lambda wildcards: config[wildcards.subject]["tumor"]
	shell:
		"""
		perl ~/programs/mskcc-vcf2maf-754d68a/vcf2maf.pl --inhibit-vep --input-vcf {input} --output-maf {output} --ref-fasta /data/CEM/shared/public_data/references/VEP_annotations/homo_sapiens/101_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz --vep-path ~/Workspace/.conda/envs/vep_env/bin/vep --ncbi-build GRCh38 --tumor-id {params.sample}
		"""


