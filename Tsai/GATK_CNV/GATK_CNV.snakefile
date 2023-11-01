# Setting up filesnames here:

from os.path import join
import os

sample = ["Te1AK",
        "Te2AK",
        "Te2SCC",
        "Te3AK",
        "Te3SCC",
        "Te4AK",
        "Te4SCC1",
        "Te4SCC2",
        "Te5AK",
        "Te5SCC",
        "Te6SCC",
        "Te8SCC",
        "Te10AK",
        "Te12AK"]
normal = ["Te1NS",
        "Te1saliva",
        "Te2NS",
        "Te3NS",
        "Te4NS",
        "Te4saliva",
        "Te5NS",
        "Te5saliva",
        "Te6NS",
        "Te6saliva",
        "Te8NS",
        "Te8saliva",
        "Te10NS",
        "Te10saliva",
        "Te12saliva"]

# Reference genome
ref_dir = "/data/CEM/shared/public_data/references/1000genomes_GRCh38_reference_genome"
ref_basename = "GRCh38_full_analysis_set_plus_decoy_hla"

# Directory Pathways
bam_dir = "/scratch/eknodel/TSAI/processed_bams/"

rule all:
    input:
        #"targets_C.preprocessed.interval_list",
        expand("cnv_files/{subject}_tumor.counts.hdf5", subject=sample),
        expand("cnv_files/{normal}_normal.counts.hdf5", normal=normal),
        expand("cnv_files/{subject}_tumor.standardizedCR.tsv", subject=sample),
        expand("cnv_files/{normal}_normal.standardizedCR.tsv", normal=normal),
        expand("cnv_plots/{subject}_tumor.denoised.png", subject=sample),
        expand("cnv_plots/{normal}_normal.denoised.png", normal=normal),
    
rule process_intervals:
    input:
        ref = os.path.join(ref_dir, ref_basename + ".fa")
    output:
        intervals = "targets_C.preprocessed.interval_list"
    shell:
        """
        /home/eknodel/gatk-4.1.7.0/gatk PreprocessIntervals \
        -L /data/CEM/shared/public_data/references/ensemble_GRCm38.68/GRCm38_68.interval_list \
        -R {input.ref} \
        --bin-length 0 \
        --padding 250 \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O {output.intervals}
        """
        #--interval-merging-rule OVERLAPPING_ONLY \

rule collect_fragment_counts:
    input: 
        bam = os.path.join(bam_dir, "{subject}." + ref_basename + ".sorted.mkdup.bam"), 
        intervals = "targets_C.preprocessed.interval_list",
        ref = os.path.join(ref_dir, ref_basename + ".fa")
    output:
        out = "cnv_files/{subject}_tumor.counts.hdf5"
    shell:
        """
        /home/eknodel/gatk-4.1.7.0/gatk CollectReadCounts \
        -I {input.bam} \
        -L {input.intervals} \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O {output.out}
        """

rule collect_fragment_counts_normal: 
    input: 
        bam = os.path.join(bam_dir, "{normal}." + ref_basename + ".sorted.mkdup.bam"),
        intervals = "targets_C.preprocessed.interval_list",
        ref = os.path.join(ref_dir, ref_basename + ".fa")
    output: 
        out = "cnv_files/{normal}_normal.counts.hdf5"
    shell:
        """
        /home/eknodel/gatk-4.1.7.0/gatk CollectReadCounts \
        -I {input.bam} \
        -L {input.intervals} \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O {output.out}
        """

rule panel_of_normals: 
    input:
        Te1NS = "cnv_files/Te1NS_normal.counts.hdf5",
        Te1saliva = "cnv_files/Te1saliva_normal.counts.hdf5",
        Te2NS = "cnv_files/Te2NS_normal.counts.hdf5",
        Te3NS = "cnv_files/Te3NS_normal.counts.hdf5",
        Te4NS = "cnv_files/Te4NS_normal.counts.hdf5",
        Te4saliva = "cnv_files/Te4saliva_normal.counts.hdf5",
        Te5NS = "cnv_files/Te5NS_normal.counts.hdf5",
        Te5saliva = "cnv_files/Te5saliva_normal.counts.hdf5",
        Te6NS = "cnv_files/Te6NS_normal.counts.hdf5",
        Te6saliva = "cnv_files/Te6saliva_normal.counts.hdf5",
        Te8NS = "cnv_files/Te8NS_normal.counts.hdf5",
        Te8saliva = "cnv_files/Te8saliva_normal.counts.hdf5",
        Te10NS = "cnv_files/Te10NS_normal.counts.hdf5",
        Te10saliva = "cnv_files/Te10saliva_normal.counts.hdf5",
        Te12saliva = "cnv_files/Te12saliva_normal.counts.hdf5"
    output: 
        "cnv_files/pon.hdf5"
    shell: 
        """
        /home/eknodel/gatk-4.1.7.0/gatk --java-options "-Xmx12g" CreateReadCountPanelOfNormals \
        -I {input.Te1NS} \
        -I {input.Te2NS} \
        -O {output} \
        """

rule standardize_tumor:
    input: 
        tumor = "cnv_files/{subject}_tumor.counts.hdf5",
        pon = "cnv_files/pon.hdf5"
    output:
        standardized = "cnv_files/{subject}_tumor.standardizedCR.tsv", 
        denoised = "cnv_files/{subject}_tumor.denoisedCR.tsv"
    shell:
        """
        /home/eknodel/gatk-4.1.7.0/gatk --java-options "-Xmx12g" DenoiseReadCounts \
        -I {input.tumor} \
        --count-panel-of-normals {input.pon} \
        --standardized-copy-ratios {output.standardized} \
        --denoised-copy-ratios {output.denoised}
        """

rule standardize_normal:
    input:
        normal = "cnv_files/{normal}_normal.counts.hdf5",
        pon = "cnv_files/pon.hdf5"
    output:
        standardized = "cnv_files/{normal}_normal.standardizedCR.tsv",
        denoised = "cnv_files/{normal}_normal.denoisedCR.tsv"
    shell:
        """
        /home/eknodel/gatk-4.1.7.0/gatk --java-options "-Xmx12g" DenoiseReadCounts \
        -I {input.normal} \
        --count-panel-of-normals {input.pon} \
        --standardized-copy-ratios {output.standardized} \
        --denoised-copy-ratios {output.denoised}
        """

rule plot_tumor:
    input:
        standardized = "cnv_files/{subject}_tumor.standardizedCR.tsv",
        denoised = "cnv_files/{subject}_tumor.denoisedCR.tsv"
    output:
        "cnv_plots/{subject}_tumor.denoised.png"
    params:
        dict = "/data/CEM/shared/public_data/references/1000genomes_GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.dict",
        prefix = "{subject}_tumor"
    shell:
        """
        /home/eknodel/gatk-4.1.7.0/gatk PlotDenoisedCopyRatios \
        --standardized-copy-ratios {input.standardized} \
        --denoised-copy-ratios {input.denoised} \
        --sequence-dictionary {params.dict} \
        --minimum-contig-length 1000000 \
        --output cnv_plots \
        --output-prefix {params.prefix}
        """

rule plot_normal:
    input:
        standardized = "cnv_files/{normal}_normal.standardizedCR.tsv",
        denoised = "cnv_files/{normal}_normal.denoisedCR.tsv"
    output:
        "cnv_plots/{normal}_normal.denoised.png"
    params:
        dict = "/data/CEM/shared/public_data/references/1000genomes_GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.dict",
        prefix = "{normal}_normal"
    shell:
        """
        /home/eknodel/gatk-4.1.7.0/gatk PlotDenoisedCopyRatios \
        --standardized-copy-ratios {input.standardized} \
        --denoised-copy-ratios {input.denoised} \
        --sequence-dictionary {params.dict} \
        --output cnv_plots \
        --output-prefix {params.prefix}
        """

