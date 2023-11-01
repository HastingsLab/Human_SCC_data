#! importing join
from os.path import join

# Workflow for quasi-quantification of RNAseq read counts with Salmon in non-alignment-based mode.

# Configuration file
configfile: "SCC.config.json"

rule all:
    input:
        expand("/scratch/eknodel/TSAI/HLA/{sample}/winners.hla.txt", sample=config["all_samples"]),
        expand("/scratch/eknodel/TSAI/HLA/{sample}/hla_types.out", sample=config["all_samples"])

rule polysolver:
    input:
        bam = os.path.join("/scratch/eknodel/TSAI/processed_bams/{sample}." + config["ref_basename"] + ".sorted.mkdup.bam")
    output:
        output = os.path.join("/scratch/eknodel/TSAI/HLA/{sample}/winners.hla.txt"),
    params:
        bam_dir = "/scratch/eknodel/TSAI/processed_bams/", # Change path - project specific
        in_file = "/data/{sample}." + config["ref_basename"] + ".sorted.mkdup.bam", # DO NOT CHANGE PATH - singularity specific
        temp_dir = "/scratch/eknodel/tmp/", # Change path - project specific
        out_dir = "/scratch/eknodel/TSAI/HLA/", # Change path - project specific
        out_file = "/out/{sample}" # DO NOT CHANGE PATH - singularity specific
    message: "Identifying HLA types for {wildcards.sample} with polysolver."
    run:
        shell("singularity exec -C -B {params.bam_dir}:/data -B {params.temp_dir}:/tmp -B {params.out_dir}:/out /home/eknodel/polysolver-singularity_v4.sif /home/polysolver/scripts/shell_call_hla_type {params.in_file} Unknown 1 hg38 STDFQ 0 {params.out_file}")

rule format: 
    input:
        os.path.join("/scratch/eknodel/TSAI/HLA/{sample}/winners.hla.txt")
    output:
        os.path.join("/scratch/eknodel/TSAI/HLA/{sample}/hla_types.out")
    shell:
        """ 
        cat {input} | sed 's/\\t/\\n/g' | sed '/HLA/d' | sed 's/hla/HLA/g' | sed 's/a/A/g' | sed 's/b/B/g' | sed 's/c/C/g' | sed 's/_/\\t/g' | awk '{{print $1"-"$2""$3":"$4}}' | sed -z 's/\\n/,/g' | sed 's/,$/"\\n/' | sed 's/^/"hla": "/' > {output}
        """
