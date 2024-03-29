import os

configfile: "SCC.config.json"

adapter_path = "/home/eknodel/Cancer_Genomics/00_misc/adapter_sequence.fa"
perl5lib_path = "/packages/6x/vcftools/0.1.12b/lib/per15/site_perl/"

rule all:
    input:
        "raw_multiqc_results/multiqc_report.html", #raw multiqc report
        "trimmed_multiqc_results/multiqc_report.html", #trimmed multiqc report

rule untar:
    input:
        fq_1 = os.path.join(config["fastq_path"], "{sample}_1.fq.tar.gz"),
        fq_2 = os.path.join(config["fastq_path"], "{sample}_2.fq.tar.gz")
    output:
        fq_1 = os.path.join(config["fastq_path"], "{sample}_1.fq"),
        fq_2 = os.path.join(config["fastq_path"], "{sample}_2.fq")
    params:
        dir = config["fastq_path"]
    shell:
        """
        tar -zxvf {input.fq_1} -C {params.dir};
        tar -zxvf {input.fq_2} -C {params.dir}
        """

rule zip: 
    input:
        fq_1 = os.path.join(config["fastq_path"], "{sample}_1.fq"),
        fq_2 = os.path.join(config["fastq_path"], "{sample}_2.fq")
    output:
        fq_1 = os.path.join(config["fastq_path"], "{sample}_1.fq.gz"),
        fq_2 = os.path.join(config["fastq_path"], "{sample}_2.fq.gz")
    shell:
        """
        gzip -c {input.fq_1} > {output.fq_1};
        gzip -c {input.fq_2} > {output.fq_2}
        """

rule fastqc_analysis:
    input:
        fq_1 = os.path.join(config["fastq_path"], "{sample}_1.fq.gz"),
        fq_2 = os.path.join(config["fastq_path"], "{sample}_2.fq.gz")
    output:
        "raw_fastqc_results/{sample}_1_fastqc.html",
        "raw_fastqc_results/{sample}_2_fastqc.html"
    params:
        perl5lib = perl5lib_path
    shell:
        """
        PERL5LIB={params.perl5lib} fastqc -o raw_fastqc_results {input.fq_1};
        PERL5LIB={params.perl5lib} fastqc -o raw_fastqc_results {input.fq_2}
        """

rule multiqc_analysis:
	input:
		expand(
			"raw_fastqc_results/{sample}_{read}_fastqc.html",
			sample=config["all_samples"],
			read=["1", "2"])
	output:
		"raw_multiqc_results/multiqc_report.html"
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"multiqc --interactive -f "
		"-o raw_multiqc_results raw_fastqc_results"

rule trim_adapters_paired_bbduk:
    input:
        fq_1 = os.path.join(config["fastq_path"], "{sample}_1.fq.gz"),
        fq_2 = os.path.join(config["fastq_path"], "{sample}_2.fq.gz")   
    output:
        out_fq_1 = "/scratch/eknodel/TSAI/trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
        out_fq_2 = "/scratch/eknodel/TSAI/trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
    params:
        adapter = adapter_path
    threads:
        2
    shell:
        "bbduk.sh -Xmx3g in1={input.fq_1} in2={input.fq_2} out1={output.out_fq_1} out2={output.out_fq_2} ref={params.adapter} qtrim=rl trimq=30 minlen=75 maq=20"

rule trimmed_fastqc_analysis:
    input:
        fq_1 = "/scratch/eknodel/TSAI/trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
        fq_2 = "/scratch/eknodel/TSAI/trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
    output:
        fq1_fastqc = "trimmed_fastqc_results/{sample}_trimmed_read1_fastqc.html",
        fq2_fastqc = "trimmed_fastqc_results/{sample}_trimmed_read2_fastqc.html"
    params: 
        perl5lib = perl5lib_path
    shell:
        """
        PERL5LIB={params.perl5lib} fastqc -o trimmed_fastqc_results {input.fq_1};
        PERL5LIB={params.perl5lib} fastqc -o trimmed_fastqc_results {input.fq_2}
        """

rule trimmed_multiqc_analysis:
	input:
		expand(
			"trimmed_fastqc_results/{sample}_trimmed_{read}_fastqc.html",
			sample=config["all_samples"],
			read=["read1", "read2"])
	output:
		"trimmed_multiqc_results/multiqc_report.html"
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"multiqc --interactive -f "
		"-o trimmed_multiqc_results trimmed_fastqc_results"
