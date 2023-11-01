#! importing join
from os.path import join

# aligning files with HISAT2, sorting and indexing BAM files with Samtools.

# Tools
HISAT2 = "hisat2"
SAMTOOLS = "samtools"

# Directories
FQ_DIR = "/scratch/eknodel/TSAI/RNA_fastqs/trimmed_fastqs/" # path to directory with trimmed FASTQ files
SAM_AL_DIR = "/scratch/eknodel/TSAI/RNA_fastqs/sam/" # path to directory for SAM alignment files
BAM_AL_DIR = "/scratch/eknodel/TSAI/RNA_fastqs/bam/" # path to directory for BAM alignment files
SORTED_BAM_AL_DIR = "/scratch/eknodel/TSAI/RNA_fastqs/sorted_bam/" # path to directory for sorted BAM alignment files
EXPRESSION_DIR = "/scratch/eknodel/TSAI/RNA_fastqs/feature_counts/" # path to directory for feature counts

# Samples
SAMPLES = ["SRR3877297",
        "SRR3877306",
        "SRR3877315",
        "SRR3877298",
        "SRR3877307",
        "SRR3877316",
        "SRR3877299",
        "SRR3877308",
        "SRR3877317",
        "SRR3877300",
        "SRR3877309",
        "SRR3877318",
        "SRR3877310",
        "SRR3877319",
        "SRR3877302",
        "SRR3877311",
        "SRR3877320",
        "SRR3877303",
        "SRR3877312",
        "SRR3877321",
        "SRR3877304",
        "SRR3877313",
        "SRR3877322",
        "SRR3877305",
        "SRR3877314"]

rule all:
	input:
		# Defining the files that snakemake will attempt to produce as an output.
		# If there is no rule defined to produce the file, or if the file already
		# exists, snakemake will throw "Nothing to be done"
		expand(SAM_AL_DIR + "{sample}_RNA_HISAT2_genome_aligned_ensembl.sam", SAM_AL_DIR=SAM_AL_DIR, sample=SAMPLES),
		expand(BAM_AL_DIR + "{sample}_RNA_HISAT2_genome_aligned_ensembl.bam", BAM_AL_DIR=BAM_AL_DIR, sample=SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_RNA_HISAT2_genome_aligned_ensembl_sortedbycoord_mkdup_RG.bam", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_RNA_HISAT2_genome_aligned_ensembl_sortedbycoord_mkdup_RG.bam.bai", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=SAMPLES),
        expand(EXPRESSION_DIR + "{sample}_HISAT2_stranded_featurecounts_ensembl.tsv", EXPRESSION_DIR=EXPRESSION_DIR, sample=SAMPLES)


rule hisat2_align_reads:
    input:
        R1 = os.path.join(FQ_DIR, "{sample}_trimmed_read1.fastq.gz"),
        R2 = os.path.join(FQ_DIR, "{sample}_trimmed_read2.fastq.gz")
    output:
        SAM = os.path.join(SAM_AL_DIR, "{sample}_RNA_HISAT2_genome_aligned_ensembl.sam")
    params:
        index = "/data/CEM/shared/public_data/references/ensemble_GRCh38.89/wholeGenome/GRCh38_wholeGenome_reference_HISAT2",
        threads = 8
    run:
        shell("hisat2 -q --phred33 -p {params.threads} -x {params.index} -s no -1 {input.R1} -2 {input.R2} -S {output.SAM}")

rule sam_to_bam:
	input:
		SAM = os.path.join(SAM_AL_DIR, "{sample}_RNA_HISAT2_genome_aligned_ensembl.sam"),
	output:
		BAM = os.path.join(BAM_AL_DIR, "{sample}_RNA_HISAT2_genome_aligned_ensembl.bam"),
	params:
	message: "Converting SAM to BAM, only outputting mapped reads."
	run:
		shell("samtools view -b -F 4 {input.SAM} > {output.BAM}")

rule sort_bam:
    input:
        BAM = os.path.join(BAM_AL_DIR, "{sample}_RNA_HISAT2_genome_aligned_ensembl.bam"),
    output:
        SORTED_BAM = os.path.join(SORTED_BAM_AL_DIR, "{sample}_RNA_HISAT2_genome_aligned_ensembl_sortedbycoord.bam"),
    params:
    message: "Sorting BAM file."
    run:
        shell("samtools sort -O bam -o {output.SORTED_BAM} {input.BAM}")

rule index_bam:
    input:
        BAM = os.path.join(SORTED_BAM_AL_DIR, "{sample}_RNA_HISAT2_genome_aligned_ensembl_sortedbycoord.bam"),
    output:
        os.path.join(SORTED_BAM_AL_DIR, "{sample}_RNA_HISAT2_genome_aligned_ensembl_sortedbycoord.bam.bai"),
    message: "Indexing sorted BAM file."
    params:
    run:
        for x in input:
            shell("samtools index {x}")

rule mark_duplicates:
    input:
        BAM = os.path.join(SORTED_BAM_AL_DIR, "{sample}_RNA_HISAT2_genome_aligned_ensembl_sortedbycoord.bam")
    output:
        BAM = os.path.join(SORTED_BAM_AL_DIR, "{sample}_RNA_HISAT2_genome_aligned_ensembl_sortedbycoord_mkdup.bam"),
        metrics = os.path.join("{sample}.picard_mkdup_metrics.txt")
    threads: 4
    shell:
        "picard -Xmx14g MarkDuplicates I={input.BAM} O={output.BAM} M={output.metrics}"

rule readgroup:
    input:
        os.path.join(SORTED_BAM_AL_DIR, "{sample}_RNA_HISAT2_genome_aligned_ensembl_sortedbycoord_mkdup.bam")
    output:
        os.path.join(SORTED_BAM_AL_DIR, "{sample}_RNA_HISAT2_genome_aligned_ensembl_sortedbycoord_mkdup_RG.bam")
    params:
        sample="{sample}"
    shell:
        """
        picard AddOrReplaceReadGroups I={input} O={output} RGID={params.sample} RGLB={params.sample} RGPL={params.sample} RGPU={params.sample} RGSM={params.sample}
        """

rule index_mkdup_bam:
	input:
		BAM = os.path.join(SORTED_BAM_AL_DIR, "{sample}_RNA_HISAT2_genome_aligned_ensembl_sortedbycoord_mkdup_RG.bam"),
	output: 
		os.path.join(SORTED_BAM_AL_DIR, "{sample}_RNA_HISAT2_genome_aligned_ensembl_sortedbycoord_mkdup_RG.bam.bai"),
	message: "Indexing sorted BAM file."
	params:
	run:
		for x in input:
			shell("samtools index {x}")

rule feature_counts:
    input: 
        BAM = os.path.join(SORTED_BAM_AL_DIR, "{sample}_RNA_HISAT2_genome_aligned_ensembl_sortedbycoord_mkdup_RG.bam"),
        #gtf = "/data/CEM/shared/public_data/references/GENCODE/gencode.v29.annotation.gene.gtf"  
        gtf = "/data/CEM/shared/public_data/references/ensemble_GRCh38.89/wholeGenome/Homo_sapiens.GRCh38.89.gtf"
    output:
        counts = os.path.join(EXPRESSION_DIR, "{sample}_HISAT2_stranded_featurecounts_ensembl.tsv")
    params:
        threads = 5
    message: "Quantifying read counts from BAM file {input.BAM} with Subread featureCounts"
    shell: 
        """
        featureCounts -T {params.threads} -O --primary -p -s 2 -t gene -g gene_id -a {input.gtf} -o {output.counts} {input.BAM}
        """
