import os

configfile: "SCC.config.json"

rule all:
    input:
        expand("/scratch/eknodel/TSAI/strelka/{subject}_ns/runWorkflow.py", subject=config["patients_with_ns"]),
        expand("/scratch/eknodel/TSAI/strelka/{subject}_ns/results/variants/somatic.snvs.vcf.gz", subject=config["patients_with_ns"]),
        expand("/scratch/eknodel/TSAI/strelka/{subject}_ns/results/variants/somatic.snvs.pass.vcf.gz", subject=config["patients_with_ns"]),
        expand("/scratch/eknodel/TSAI/strelka/{subject_sal}_saliva/runWorkflow.py", subject_sal=config["patients_with_saliva"]),
        expand("/scratch/eknodel/TSAI/strelka/{subject_sal}_saliva/results/variants/somatic.snvs.vcf.gz", subject_sal=config["patients_with_saliva"]),
        expand("/scratch/eknodel/TSAI/strelka/{subject_sal}_saliva/results/variants/somatic.snvs.pass.vcf.gz", subject_sal=config["patients_with_saliva"]),


ruleorder: config_with_matched_normal_skin > run_with_matched_normal_skin

ruleorder: config_with_saliva > run_with_saliva

rule config_with_matched_normal_skin:
    input:
        normal_bam = lambda wildcards: os.path.join(
			"/scratch/eknodel/TSAI/processed_bams/", config[wildcards.subject]["normal_skin_sample"] + "." + config["ref_basename"] + ".sorted.bam"),
        tumor_bam = lambda wildcards: os.path.join(
			"/scratch/eknodel/TSAI/processed_bams/", config[wildcards.subject]["tumor_sample"] + "." + config["ref_basename"] + ".sorted.bam"),
        ref = os.path.join(config["ref_dir"], config["ref_basename"] + ".fa")
    output:
        "/scratch/eknodel/TSAI/strelka/{subject}_ns/runWorkflow.py"
    params:
        strelka = config["strelka"],
        run_dir = "/scratch/eknodel/TSAI/strelka/{subject}_ns"
    shell:
        """
        {params.strelka} --normalBam {input.normal_bam} --tumorBam {input.tumor_bam} --referenceFasta {input.ref} --runDir {params.run_dir}
        """

rule run_with_matched_normal_skin:
    input:
        normal_bam = lambda wildcards: os.path.join(
			"/scratch/eknodel/TSAI/processed_bams/", config[wildcards.subject]["normal_skin_sample"] + "." + config["ref_basename"] + ".sorted.bam"),
        tumor_bam = lambda wildcards: os.path.join(
			"/scratch/eknodel/TSAI/processed_bams/", config[wildcards.subject]["tumor_sample"] + "." + config["ref_basename"] + ".sorted.bam"),
        ref = os.path.join(config["ref_dir"], config["ref_basename"] + ".fa")
    output:
        snvs = "/scratch/eknodel/TSAI/strelka/{subject}_ns/results/variants/somatic.snvs.vcf.gz",
        indels = "/scratch/eknodel/TSAI/strelka/{subject}_ns/results/variants/somatic.indels.vcf.gz"
    params:
        run = "/scratch/eknodel/TSAI/strelka/{subject}_ns/runWorkflow.py"
    shell:
        """
        {params.run} -m local -j 20
        """

rule select_pass_variants_with_matched_normal_skin:
    input:
        ref = os.path.join(config["ref_dir"], config["ref_basename"] + ".fa"),
        snvs = "/scratch/eknodel/TSAI/strelka/{subject}_ns/results/variants/somatic.snvs.vcf.gz",
        indels = "/scratch/eknodel/TSAI/strelka/{subject}_ns/results/variants/somatic.indels.vcf.gz"
    output:
        snvs = "/scratch/eknodel/TSAI/strelka/{subject}_ns/results/variants/somatic.snvs.pass.vcf.gz",
        indels = "/scratch/eknodel/TSAI/strelka/{subject}_ns/results/variants/somatic.indels.pass.vcf.gz"
    params:
        gatk = config["gatk_path"]
    shell:
        """
        {params.gatk} SelectVariants -R {input.ref} -V {input.snvs} --exclude-filtered -O {output.snvs};
        {params.gatk} SelectVariants -R {input.ref} -V {input.indels} --exclude-filtered -O {output.indels}
        """

rule gunzip_with_matched_normal_skin:
    input:
        snvs = "/scratch/eknodel/TSAI/strelka/{subject}_ns/results/variants/somatic.snvs.pass.vcf.gz",
        indels = "/scratch/eknodel/TSAI/strelka/{subject}_ns/results/variants/somatic.indels.pass.vcf.gz"
    output:
        snvs = "/scratch/eknodel/TSAI/strelka/{subject}_ns/results/variants/somatic.snvs.pass.vcf",
        indels = "/scratch/eknodel/TSAI/strelka/{subject}_ns/results/variants/somatic.indels.pass.vcf"
    shell:
        """
        gunzip -c {input.snvs} > {output.snvs};
        gunzip -c {input.indels} > {output.indels}
        """

rule config_with_saliva:
    input:
        normal_bam = lambda wildcards: os.path.join(
                        "/scratch/eknodel/TSAI/processed_bams/", config[wildcards.subject_sal]["saliva_sample"] + "." + config["ref_basename"] + ".sorted.bam"),
        tumor_bam = lambda wildcards: os.path.join(
                        "/scratch/eknodel/TSAI/processed_bams/", config[wildcards.subject_sal]["tumor_sample"] + "." + config["ref_basename"] + ".sorted.bam"),
        ref = os.path.join(config["ref_dir"], config["ref_basename"] + ".fa")
    output:
        "/scratch/eknodel/TSAI/strelka/{subject_sal}_saliva/runWorkflow.py"
    params:
        strelka = config["strelka"],
        run_dir = "/scratch/eknodel/TSAI/strelka/{subject_sal}_saliva"
    shell:
        """
        {params.strelka} --normalBam {input.normal_bam} --tumorBam {input.tumor_bam} --referenceFasta {input.ref} --runDir {params.run_dir}
        """

rule run_with_saliva:
    input:
        normal_bam = lambda wildcards: os.path.join(
                        "/scratch/eknodel/TSAI/processed_bams/", config[wildcards.subject_sal]["saliva_sample"] + "." + config["ref_basename"] + ".sorted.bam"),
        tumor_bam = lambda wildcards: os.path.join(
                        "/scratch/eknodel/TSAI/processed_bams/", config[wildcards.subject_sal]["tumor_sample"] + "." + config["ref_basename"] + ".sorted.bam"),
        ref = os.path.join(config["ref_dir"], config["ref_basename"] + ".fa")
    output:
        snvs = "/scratch/eknodel/TSAI/strelka/{subject_sal}_saliva/results/variants/somatic.snvs.vcf.gz",
        indels = "/scratch/eknodel/TSAI/strelka/{subject_sal}_saliva/results/variants/somatic.indels.vcf.gz"
    params:
        run = "/scratch/eknodel/TSAI/strelka/{subject_sal}_saliva/runWorkflow.py"
    shell:
        """
        {params.run} -m local -j 20
        """

rule select_pass_variants_with_saliva:
    input:
        ref = os.path.join(config["ref_dir"], config["ref_basename"] + ".fa"),
        snvs = "/scratch/eknodel/TSAI/strelka/{subject_sal}_saliva/results/variants/somatic.snvs.vcf.gz",
        indels = "/scratch/eknodel/TSAI/strelka/{subject_sal}_saliva/results/variants/somatic.indels.vcf.gz"
    output:
        snvs = "/scratch/eknodel/TSAI/strelka/{subject_sal}_saliva/results/variants/somatic.snvs.pass.vcf.gz",
        indels = "/scratch/eknodel/TSAI/strelka/{subject_sal}_saliva/results/variants/somatic.indels.pass.vcf.gz"
    params:
        gatk = config["gatk_path"]
    shell:
        """
        {params.gatk} SelectVariants -R {input.ref} -V {input.snvs} --exclude-filtered -O {output.snvs};
        {params.gatk} SelectVariants -R {input.ref} -V {input.indels} --exclude-filtered -O {output.indels}
        """

rule gunzip_with_saliva:
    input:
        snvs = "/scratch/eknodel/TSAI/strelka/{subject_sal}_saliva/results/variants/somatic.snvs.pass.vcf.gz",
        indels = "/scratch/eknodel/TSAI/strelka/{subject_sal}_saliva/results/variants/somatic.indels.pass.vcf.gz"
    output:
        snvs = "/scratch/eknodel/TSAI/strelka/{subject_sal}_saliva/results/variants/somatic.snvs.pass.vcf",
        indels = "/scratch/eknodel/TSAI/strelka/{subject_sal}_saliva/results/variants/somatic.indels.pass.vcf"
    shell:
        """
        gunzip -c {input.snvs} > {output.snvs};
        gunzip -c {input.indels} > {output.indels}
        """
