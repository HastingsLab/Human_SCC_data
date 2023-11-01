import os

configfile: "RNA_SCC.config.json"

rule all:
    input:
        expand("{sample}_sortme.output", sample="SRR3877297")

rule trimmed_multiqc_analysis:
        input:
            R1 = os.path.join(config["fastq_path"], "{sample}_1.fastq"),
            R2 = os.path.join(config["fastq_path"], "{sample}_2.fastq"),
            ref = os.path.join(config["ref_dir"], config["ref_basename"] + ".fa")
        output:
            "{sample}_sortme.output"  
        params:
            workdir = "/scratch/eknodel/TSAI/RNA_fastqs/sortme/{sample}"      
        shell:
            """
            sortmerna --ref {input.ref} --reads {input.R1} --reads {input.R2} --workdir {params.workdir}
            """
