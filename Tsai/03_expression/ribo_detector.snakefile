import os

configfile: "RNA_SCC.config.json"

rule all:
    input:
        expand("/scratch/eknodel/TSAI/RNA_fastqs/trimmed_fastqs/{sample}_trimmed_read1_rrna_removed.fastq.gz", sample=config["all_samples"])

rule ribo_detector:
        input:
            R1 = "/scratch/eknodel/TSAI/RNA_fastqs/trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
            R2 = "/scratch/eknodel/TSAI/RNA_fastqs/trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
        output:
            R1 = "/scratch/eknodel/TSAI/RNA_fastqs/trimmed_fastqs/{sample}_trimmed_read1_rrna_removed.fastq.gz",
            R2 = "/scratch/eknodel/TSAI/RNA_fastqs/trimmed_fastqs/{sample}_trimmed_read2_rrna_removed.fastq.gz"
        shell:
            """
            ribodetector_cpu -t 20 -l 75 -i {input.R1} {input.R2} -e rrna --chunk_size 256 -o {output.R1} {output.R2}
            """
