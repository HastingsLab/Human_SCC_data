import os

configfile: "SCC.config.json"

rule all:
    input:
        expand(os.path.join("/scratch/eknodel/TSAI/gatk_mutect2/", "{subject}.somatic.ns.vcf.gz"), subject=config["patients_with_ns"]),
        expand(os.path.join("/scratch/eknodel/TSAI/gatk_mutect2", "{subject}.somatic.ns.filtered.vcf.gz"), subject=config["patients_with_ns"]),
        expand(os.path.join("/scratch/eknodel/TSAI/gatk_mutect2/", "{subject}.somatic.ns.filtered.pass.vcf.gz"), subject=config["patients_with_ns"]),
        expand(os.path.join("/scratch/eknodel/TSAI/gatk_mutect2/", "{subject_sal}.somatic.saliva.vcf.gz"), subject_sal=config["patients_with_saliva"]),
        expand(os.path.join("/scratch/eknodel/TSAI/gatk_mutect2", "{subject_sal}.somatic.saliva.filtered.vcf.gz"), subject_sal=config["patients_with_saliva"]),
        expand(os.path.join("/scratch/eknodel/TSAI/gatk_mutect2/", "{subject_sal}.somatic.saliva.filtered.pass.vcf.gz"), subject_sal=config["patients_with_saliva"])

rule tumor_with_matched_normal_skin:
    input:
        ref = os.path.join(config["ref_dir"], config["ref_basename"] + ".fa"),
        normal_bam = lambda wildcards: os.path.join(
			"/scratch/eknodel/TSAI/processed_bams/", config[wildcards.subject]["normal_skin_sample"] + "." + config["ref_basename"] + ".sorted.bam"),
        tumor_bam = lambda wildcards: os.path.join(
			"/scratch/eknodel/TSAI/processed_bams/", config[wildcards.subject]["tumor_sample"] + "." + config["ref_basename"] + ".sorted.bam")
    output:
        os.path.join("/scratch/eknodel/TSAI/gatk_mutect2/", "{subject}.somatic.ns.vcf.gz")
    params:
        gatk = config["gatk_path"],
        sm = lambda wildcards: config[wildcards.subject]["normal_skin_sample"]
    shell:
        """
        {params.gatk} Mutect2 -R {input.ref} -I {input.tumor_bam} -I {input.normal_bam} -normal {params.sm} -O {output}
        """

rule filter_tumor_with_matched_normal_skin:
    input:
        ref = os.path.join(config["ref_dir"], config["ref_basename"] + ".fa"),
        unfiltered = os.path.join("/scratch/eknodel/TSAI/gatk_mutect2/", "{subject}.somatic.ns.vcf.gz")
    output:
        filtered = os.path.join("/scratch/eknodel/TSAI/gatk_mutect2/", "{subject}.somatic.ns.filtered.vcf.gz")
    params:
        gatk = config["gatk_path"]
    shell:
        """
        {params.gatk} FilterMutectCalls -R {input.ref} -V {input.unfiltered} -O {output.filtered}
        """

rule select_pass_variants_tumor_with_matched_normal_skin:
    input:
        ref = os.path.join(config["ref_dir"], config["ref_basename"] + ".fa"),
        vcf = os.path.join("/scratch/eknodel/TSAI/gatk_mutect2/", "{subject}.somatic.ns.filtered.vcf.gz")
    output:
        os.path.join("/scratch/eknodel/TSAI/gatk_mutect2/", "{subject}.somatic.ns.filtered.pass.vcf.gz")
    params:
        gatk = config["gatk_path"]
    shell:
        """
        {params.gatk} SelectVariants -R {input.ref} -V {input.vcf} --exclude-filtered -O {output}
        """

rule tumor_with_matched_saliva:
    input:
        ref = os.path.join(config["ref_dir"], config["ref_basename"] + ".fa"),
        normal_bam = lambda wildcards: os.path.join(
            "/scratch/eknodel/TSAI/processed_bams/", config[wildcards.subject_sal]["saliva_sample"] + "." + config["ref_basename"] + ".sorted.bam"),
        tumor_bam = lambda wildcards: os.path.join(
            "/scratch/eknodel/TSAI/processed_bams/", config[wildcards.subject_sal]["tumor_sample"] + "." + config["ref_basename"] + ".sorted.bam")
    output:
        os.path.join("/scratch/eknodel/TSAI/gatk_mutect2/", "{subject_sal}.somatic.saliva.vcf.gz")
    params:
        gatk = config["gatk_path"],
        sm = lambda wildcards: config[wildcards.subject_sal]["saliva_sample"]
    shell:
        """
        {params.gatk} Mutect2 -R {input.ref} -I {input.tumor_bam} -I {input.normal_bam} -normal {params.sm} -O {output}
        """

rule filter_tumor_with_matched_saliva:
    input:
        ref = os.path.join(config["ref_dir"], config["ref_basename"] + ".fa"),
        unfiltered = os.path.join("/scratch/eknodel/TSAI/gatk_mutect2/", "{subject_sal}.somatic.saliva.vcf.gz")
    output:
        filtered = os.path.join("/scratch/eknodel/TSAI/gatk_mutect2/", "{subject_sal}.somatic.saliva.filtered.vcf.gz")
    params:
        gatk = config["gatk_path"]
    shell:
        """
        {params.gatk} FilterMutectCalls -R {input.ref} -V {input.unfiltered} -O {output.filtered}
        """

rule select_pass_variants_tumor_with_saliva:
    input:
        ref = os.path.join(config["ref_dir"], config["ref_basename"] + ".fa"),
        vcf = os.path.join("/scratch/eknodel/TSAI/gatk_mutect2/", "{subject_sal}.somatic.saliva.filtered.vcf.gz")
    output:
        os.path.join("/scratch/eknodel/TSAI/gatk_mutect2/", "{subject_sal}.somatic.saliva.filtered.pass.vcf.gz")
    params:
        gatk = config["gatk_path"]
    shell:
        """
        {params.gatk} SelectVariants -R {input.ref} -V {input.vcf} --exclude-filtered -O {output}
        """
