# Setting up filesnames here:
from os.path import join
import os

samples = ["2", "3", "4", "5", "6", "7", "8", "9", "10"]

rule all:
    input:
        expand("variants/{sample}_merged.vcf", sample=samples), # combine all mutations
        expand("variants/{sample}_merged_vep.vcf", sample=samples), # run VEP
        expand("peptides/{sample}_vep.17.peptide", sample=samples) # run pvacseq

rule combine_strelka:
    input:
        strelka = os.path.join("/scratch/tjuan1/Ji_Cancergenomics/strelka/{sample}/results/variants/somatic.snvs.pass.vcf.gz"),
        strelka_indels = os.path.join("/scratch/tjuan1/Ji_Cancergenomics/strelka/{sample}/results/variants/somatic.indels.pass.vcf.gz") 
    output:
        strelka = os.path.join("../01_somatic_mutation_calling/strelka/{sample}/results/variants/somatic.pass.vcf"),
    shell:
        """
        bcftools concat -a {input.strelka_indels} {input.strelka} -o {output.strelka}
        """

rule sort_strelka:
     input:
         vcf = os.path.join("../01_somatic_mutation_calling/strelka/{sample}/results/variants/somatic.pass.vcf")
     output:
         vcf = os.path.join("../01_somatic_mutation_calling/strelka/{sample}/results/variants/somatic.pass_sorted.vcf")
     shell:
         """
         picard SortVcf I={input.vcf} O={output.vcf}
         """

rule gunzip:
    input:
        gatk = os.path.join("/scratch/tjuan1/Ji_Cancergenomics/gatk_mutect2/{sample}.somatic.filtered.pass.vcf.gz")
    output:
        gatk = os.path.join("../01_somatic_mutation_calling/gatk_mutect2/{sample}.somatic.filtered.pass.vcf")
    shell: 
        """
        gunzip -c {input.gatk} > {output.gatk}
        """

rule sort_gatk:
     input:
         vcf = os.path.join("../01_somatic_mutation_calling/gatk_mutect2/{sample}.somatic.filtered.pass.vcf")
     output:
         vcf = os.path.join("../01_somatic_mutation_calling/gatk_mutect2/{sample}.somatic.filtered.pass_sorted.vcf")
     shell:
         """
         picard SortVcf I={input.vcf} O={output.vcf}
         """

rule generate_input:
    input:
        gatk = os.path.join("../01_somatic_mutation_calling/gatk_mutect2/{sample}.somatic.filtered.pass_sorted.vcf"),
        strelka = os.path.join("../01_somatic_mutation_calling/strelka/{sample}/results/variants/somatic.pass_sorted.vcf")
    output:
        vep = os.path.join("variants/{sample}_merged.vcf")
    shell: 
        """
        python merge_Mutect2_Strelka2_variants.py \
        --gatk_file {input.gatk} \
        --strelka_file {input.strelka} \
        --vep_format_fn {output.vep}
        """

rule run_vep:
    input:
        os.path.join("variants/{sample}_merged.vcf")
    output:
        os.path.join("variants/{sample}_merged_vep.vcf")
    shell:
        """
        vep -i {input} --format vcf --assembly GRCh38 --cache --cache_version 101 --dir_cache /data/CEM/shared/public_data/references/VEP_annotations --offline --vcf -o {output} --force_overwrite --symbol --plugin Wildtype --plugin Frameshift --terms SO --plugin Downstream
	""" 

rule generate_fasta:
	input:
		os.path.join("variants/{sample}_merged_vep.vcf")
	output:
		len_17 = os.path.join("peptides/{sample}_vep.17.peptide")
	params:
		sample = "{sample}"
	shell:
		"""
		python /home/eknodel/programs/pVACtools/pvactools/tools/pvacseq/generate_protein_fasta.py {input} 8 {output.len_17} -s {params.sample};
		"""
