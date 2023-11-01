import os

configfile: "SCC.config.json"

rule all:
    input:
        expand(os.path.join("/scratch/eknodel/TSAI/processed_bams/{subject}_normal_footprint.out"), subject=config["patients_with_ns"])

rule footprint:
    input:
        normal_bam = lambda wildcards: os.path.join(
            "/scratch/eknodel/TSAI/processed_bams/", config[wildcards.subject]["normal_skin_sample"] + "." + config["ref_basename"] + ".sorted.bam"),
        tumor_bam = lambda wildcards: os.path.join(
            "/scratch/eknodel/TSAI/processed_bams/", config[wildcards.subject]["tumor_sample"] + "." + config["ref_basename"] + ".sorted.bam")
    output:
        normal_footprint = os.path.join("/scratch/eknodel/TSAI/processed_bams/{subject}_normal_footprint.out"),
        tumor_footprint = os.path.join("/scratch/eknodel/TSAI/processed_bams/{subject}_tumor_footprint.out")
    shell:
        """
        bedtools genomecov -ibam {input.normal_bam} -max 2 > {output.normal_footprint};
        bedtools genomecov -ibam {input.tumor_bam} -max 2 > {output.tumor_footprint}
        """
