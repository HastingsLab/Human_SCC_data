# Setting up filesnames here:
from os.path import join
import os

configfile: "SCC.config.json"

rule all:
    input:
        expand("variants/{subject}_ns_merged_vep.maf", subject=config["patients_with_ns"]), # combine all mutations
	expand("/scratch/eknodel/TSAI/gatk_mutect2/{subject}.somatic.ns.filtered.pass_sorted_vep.vcf", subject=config["patients_with_ns"]),
	expand("variants/{subject}_ns_raw_vep.maf", subject=config["patients_with_ns"])

rule maf:
	input:
		"variants/{subject}_ns_merged_vep.vcf"
	output:
		"variants/{subject}_ns_merged_vep.maf"
	params:
		sample = lambda wildcards: config[wildcards.subject]["tumor_sample"]
	shell:
		"""
		perl ~/programs/mskcc-vcf2maf-754d68a/vcf2maf.pl --inhibit-vep --input-vcf {input} --output-maf {output} --ref-fasta /data/CEM/shared/public_data/references/VEP_annotations/homo_sapiens/101_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz --vep-path ~/Workspace/.conda/envs/vep_env/bin/vep --ncbi-build GRCh38 --tumor-id {params.sample}
		"""

rule vep_raw:
	input:
		"/scratch/eknodel/TSAI/gatk_mutect2/{subject}.somatic.ns.filtered.pass_sorted.vcf"
	output:
		"/scratch/eknodel/TSAI/gatk_mutect2/{subject}.somatic.ns.filtered.pass_sorted_vep.vcf"
	shell:
		"""
		vep -i {input} --format vcf --assembly GRCh38 --cache --dir_cache ~/external_scripts --offline --vcf -o {output} --force_overwrite --plugin Wildtype --plugin Frameshift --symbol --terms SO --plugin Downstream
		"""

rule maf_raw:
	input:
		"/scratch/eknodel/TSAI/gatk_mutect2/{subject}.somatic.ns.filtered.pass_sorted_vep.vcf"
	output:
		"variants/{subject}_ns_raw_vep.maf"
	params:
		sample = lambda wildcards: config[wildcards.subject]["tumor_sample"]
	shell:
		"""
		perl ~/programs/mskcc-vcf2maf-754d68a/vcf2maf.pl --inhibit-vep --input-vcf {input} --output-maf {output} --ref-fasta /data/CEM/shared/public_data/references/VEP_annotations/homo_sapiens/101_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz --vep-path /home/eknodel/Workspace/.conda/envs/vep_env/bin/ --ncbi-build GRCh38 --tumor-id {params.sample}
		"""
