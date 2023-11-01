import os

configfile: "SCC.config.json"

rule all:
    input:
        os.path.join(config["ref_dir"], config["ref_basename"] + ".fa.fai"),
        os.path.join(config["ref_dir"], config["ref_basename"] + ".dict"),
        os.path.join(config["ref_dir"], config["ref_basename"] + ".fa.amb"),
        expand(os.path.join("/scratch/eknodel/TSAI/processed_bams/{sample}." + config["ref_basename"] + ".sorted.bam.bai"), sample=config["all_samples"]),
        expand(os.path.join("/scratch/eknodel/TSAI/processed_bams/{sample}." + config["ref_basename"] + ".sorted.mkdup.bam.bai"), sample=config["all_samples"])

rule prep_refs:
    input:
        os.path.join(config["ref_dir"], config["ref_basename"] + ".fa")
    output:
        fai = os.path.join(config["ref_dir"], config["ref_basename"] + ".fa.fai"),
        dict = os.path.join(config["ref_dir"], config["ref_basename"] + ".dict"),
        amb = os.path.join(config["ref_dir"], config["ref_basename"] + ".fa.amb")
    shell:
        """
        samtools faidx {input};
        samtools dict -o {output.dict} {input};
        bwa index {input}
        """

rule map:
    input:
        fq_1 = "/scratch/eknodel/TSAI/trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
        fq_2 = "/scratch/eknodel/TSAI/trimmed_fastqs/{sample}_trimmed_read2.fastq.gz",
        fai = os.path.join(config["ref_dir"], config["ref_basename"] + ".fa.fai"),
        ref = os.path.join(config["ref_dir"], config["ref_basename"] + ".fa")
    output:
        os.path.join("/scratch/eknodel/TSAI/processed_bams/{sample}." + config["ref_basename"] + ".sorted.bam")
    params:
        id = "{sample}",
        sm = "{sample}",
        lb = "{sample}",
        pu = "{sample}",
        pl = "Illumina"
    threads:
        4
    shell:
        "bwa mem -t {threads} -R "
        "'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
        "{input.ref} {input.fq_1} {input.fq_2}"
        " | samtools fixmate -O bam - - | samtools sort "
        "-O bam -o {output}"

rule index:
    input:
        os.path.join("/scratch/eknodel/TSAI/processed_bams/{sample}." + config["ref_basename"] + ".sorted.bam")
    output:
        os.path.join("/scratch/eknodel/TSAI/processed_bams/{sample}." + config["ref_basename"] + ".sorted.bam.bai")
    shell:
        """
        samtools index {input}
        """

rule mark_duplicates:
    input:
        bam = os.path.join("/scratch/eknodel/TSAI/processed_bams/{sample}." + config["ref_basename"] + ".sorted.bam"),
        bai = os.path.join("/scratch/eknodel/TSAI/processed_bams/{sample}." + config["ref_basename"] + ".sorted.bam.bai")
    output:
        bam = os.path.join("/scratch/eknodel/TSAI/processed_bams/{sample}." + config["ref_basename"] + ".sorted.mkdup.bam"),
        metrics = os.path.join("/scratch/eknodel/TSAI/processed_bams/{sample}.picard_mkdup_metrics.txt")
    threads: 4
    shell:
        "picard -Xmx14g MarkDuplicates I={input.bam} O={output.bam} M={output.metrics}"

rule index_mkdup:
    input:
        os.path.join("/scratch/eknodel/TSAI/processed_bams/{sample}." + config["ref_basename"] + ".sorted.mkdup.bam")
    output:
        os.path.join("/scratch/eknodel/TSAI/processed_bams/{sample}." + config["ref_basename"] + ".sorted.mkdup.bam.bai")
    shell:
        "samtools index {input}"
