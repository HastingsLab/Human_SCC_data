import os

configfile: "somatic_mutation_calling_config.json"

rule all:
    input:
        expand(os.path.join("/scratch/eknodel/Ji_CancerGenomics/processed_bams/{subject}_normal_footprint.out"), subject=config["all_subjects"])

rule footprint:
    input:
        normal_bam = lambda wildcards: os.path.join(
            "/scratch/tjuan1/Ji_Cancergenomics/processed_bams/", config[wildcards.subject]["normal"] + "." + config["ref_basename"] + ".sorted.bam"),
        tumor_bam = lambda wildcards: os.path.join(
            "/scratch/tjuan1/Ji_Cancergenomics/processed_bams/", config[wildcards.subject]["tumor"] + "." + config["ref_basename"] + ".sorted.bam")
    output:
        normal_footprint = os.path.join("/scratch/eknodel/Ji_CancerGenomics/processed_bams/{subject}_normal_footprint.out"),
        tumor_footprint = os.path.join("/scratch/eknodel/Ji_CancerGenomics/processed_bams/{subject}_tumor_footprint.out")
    shell:
        """
        bedtools genomecov -ibam {input.normal_bam} -max 2 > {output.normal_footprint};
        bedtools genomecov -ibam {input.tumor_bam} -max 2 > {output.tumor_footprint}
        """
