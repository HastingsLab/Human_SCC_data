#! importing join
from os.path import join

configfile: "SCC.config.json"

# Directories
EXPRESSION_DIR = "/scratch/eknodel/TSAI/RNA_fastqs/ASE_readcounts/"
SORTED_BAM_AL_DIR = "/scratch/eknodel/TSAI/RNA_fastqs/sorted_bam/"

# Samples
subject = subject=config["patients_with_ns"]

rule all:
	input:
        #expand("variants/{subject}_saliva_merged_vep_gt_select.vcf.idx", subject=subject),
		expand("../02_variant_annotation/variants/{subject}_ns_merged_vep_gt_select.vcf.idx", subject=subject),
        expand(EXPRESSION_DIR + "{subject}_gatk_filtered_expression.out", EXPRESSION_DIR=EXPRESSION_DIR, subject=subject)

rule remove_duplicates:
    input:    
        var = os.path.join("../02_variant_annotation/variants/{subject}_ns_merged_vep_gt.vcf"),
        ref = "/data/CEM/shared/public_data/references/1000genomes_GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    output:
        os.path.join("../02_variant_annotation/variants/{subject}_ns_merged_vep_gt_select.vcf")
    params: 
        gatk = "/home/eknodel/gatk-4.1.7.0/gatk"
    shell:
        """
        {params.gatk} SelectVariants --V {input.var} --select-type-to-include SNP -O {output};
        """

rule index_vcf:
    input:
        var = os.path.join("../02_variant_annotation/variants/{subject}_ns_merged_vep_gt_select.vcf") 
    output:
        var = os.path.join("../02_variant_annotation/variants/{subject}_ns_merged_vep_gt_select.vcf.idx")
    params: 
        gatk = "/home/eknodel/gatk-4.1.7.0/gatk"
    shell:
        """
        {params.gatk} IndexFeatureFile -I {input.var};
        """

rule gatk_asereadcounter:
    input: 
        BAM = lambda wildcards: os.path.join(SORTED_BAM_AL_DIR, config[wildcards.subject]["RNA"] + "_RNA_HISAT2_genome_aligned_sortedbycoord_mkdup_RG.bam"),
        var = os.path.join("../02_variant_annotation/variants/{subject}_ns_merged_vep_gt_select.vcf"),
        index = os.path.join("../02_variant_annotation/variants/{subject}_ns_merged_vep_gt_select.vcf.idx")
    output:
        var = os.path.join(EXPRESSION_DIR, "{subject}_gatk_filtered_expression.tsv"),
    params:
        gatk = "/home/eknodel/gatk-4.1.7.0/gatk",
        ref = "/data/CEM/shared/public_data/references/1000genomes_GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    shell: 
        """
        {params.gatk} ASEReadCounter --output {output.var} --input {input.BAM} --R {params.ref} --variant {input.var};  
        """

rule format_output:
    input:
        var = os.path.join(EXPRESSION_DIR, "{subject}_gatk_filtered_expression.tsv"),
    output:
        var = os.path.join(EXPRESSION_DIR, "{subject}_gatk_filtered_expression.out"),
    shell:
        """
        cat {input.var} | awk '{{print $1":"$2-1":"$4":"$5, $6, $7, $8, $12}}' > {output.var};
        """

