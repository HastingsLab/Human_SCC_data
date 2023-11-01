import os

configfile: "SCC.config.json"

rule all:
    input:
        expand("/scratch/eknodel/TSAI/varscan/pileups/{sample}.pileup", sample=config["all_samples"]), #run bam_pileup
        expand("/scratch/eknodel/TSAI/varscan/{subject}.ns.varscan.snp", subject=config["patients_with_ns"]), #run VarScan
        expand("/scratch/eknodel/TSAI/varscan/{subject}.ns.varscan.indel", subject=config["patients_with_ns"]), #run VarScan
        expand("/scratch/eknodel/TSAI/varscan/{subject}.ns.varscan.snp.Somatic.hc.filter.vcf", subject=config["patients_with_ns"]), #convert to vcf
        expand("/scratch/eknodel/TSAI/varscan/{subject_sal}.saliva.varscan.snp", subject_sal=config["patients_with_saliva"]), #run VarScan
        expand("/scratch/eknodel/TSAI/varscan/{subject_sal}.saliva.varscan.indel", subject_sal=config["patients_with_saliva"]), #run VarScan
        expand("/scratch/eknodel/TSAI/varscan/{subject_sal}.saliva.varscan.snp.Somatic.hc.filter.vcf", subject_sal=config["patients_with_saliva"]) #convert to vc

rule bam_pileup: #for both normal, tumor, and saliva
    input:
        ref = os.path.join(config["ref_dir"], config["ref_basename"] + ".fa"),
        bam = os.path.join("/scratch/eknodel/TSAI/processed_bams/{sample}." + config["ref_basename"] + ".sorted.bam")
    output:
        pileup = "/scratch/eknodel/TSAI/varscan/pileups/{sample}.pileup"
    threads: 4
    shell:
        """
        samtools mpileup -f {input.ref} {input.bam} > {output.pileup}
        """

rule run_varscan_with_ns:
    input:
        normal_pileup = lambda wildcards: os.path.join(
			"/scratch/eknodel/TSAI/varscan/pileups/", config[wildcards.subject]["normal_skin_sample"] +  ".pileup"),
        tumor_pileup = lambda wildcards: os.path.join(
			"/scratch/eknodel/TSAI/varscan/pileups/", config[wildcards.subject]["tumor_sample"] + ".pileup")
    output:
        snp = "/scratch/eknodel/TSAI/varscan/{subject}.ns.varscan.snp",
        indel = "/scratch/eknodel/TSAI/varscan/{subject}.ns.varscan.indel"
    params:
        varscan = config["varscan_path"],
        basename = "/scratch/eknodel/TSAI/varscan/{subject}.ns.varscan"
    threads: 4
    shell:
        "java -jar {params.varscan} somatic {input.normal_pileup} {input.tumor_pileup} {params.basename} –min-coverage 10 –min-var-freq 0.08 –somatic-p-value 0.05"

rule isolate_calls_by_type_and_confidence_with_ns:
    input:
        varscan_snp = "/scratch/eknodel/TSAI/varscan/{subject}.ns.varscan.snp"
    output:
        varscan_snp = "/scratch/eknodel/TSAI/varscan/{subject}.ns.varscan.snp.Somatic.hc"
    params:
        varscan = config["varscan_path"]
    shell:
        """
        java -jar {params.varscan} processSomatic {input.varscan_snp}
        """

rule somatic_filter_with_ns:
    input:
        snp_somatic_hc = "/scratch/eknodel/TSAI/varscan/{subject}.ns.varscan.snp.Somatic.hc",
        indel = "/scratch/eknodel/TSAI/varscan/{subject}.ns.varscan.indel"
    output:
        snp_somatic_hc_filter = "/scratch/eknodel/TSAI/varscan/{subject}.ns.varscan.snp.Somatic.hc.filter"
    params:
        varscan = config["varscan_path"]
    shell:
        """
        java -jar {params.varscan} somaticFilter {input.snp_somatic_hc} -indel-file {input.indel} -output-file {output.snp_somatic_hc_filter}
        """

rule convert_to_vcf_with_ns:
    input:
        "/scratch/eknodel/TSAI/varscan/{subject}.ns.varscan.snp.Somatic.hc.filter"
    output:
        "/scratch/eknodel/TSAI/varscan/{subject}.ns.varscan.snp.Somatic.hc.filter.vcf"
    shell:
        """
        cat {input} | awk '{{print $1"\t" $2"\t" "." "\t" $3 "\t" $4}}' > {output}
        """

rule run_varscan_with_saliva:
    input:
        normal_pileup = lambda wildcards: os.path.join(
                        "/scratch/eknodel/TSAI/varscan/pileups/", config[wildcards.subject_sal]["saliva_sample"] +  ".pileup"),
        tumor_pileup = lambda wildcards: os.path.join(
                        "/scratch/eknodel/TSAI/varscan/pileups/", config[wildcards.subject_sal]["tumor_sample"] + ".pileup")
    output:
        snp = "/scratch/eknodel/TSAI/varscan/{subject_sal}.saliva.varscan.snp",
        indel = "/scratch/eknodel/TSAI/varscan/{subject_sal}.saliva.varscan.indel"
    params:
        varscan = config["varscan_path"],
        basename = "/scratch/eknodel/TSAI/varscan/{subject_sal}.saliva.varscan"
    threads: 4
    shell:
        "java -jar {params.varscan} somatic {input.normal_pileup} {input.tumor_pileup} {params.basename} –min-coverage 10 –min-var-freq 0.08 –somatic-p-value 0.05"

rule isolate_calls_by_type_and_confidence_with_saliva:
    input:
        varscan_snp = "/scratch/eknodel/TSAI/varscan/{subject_sal}.saliva.varscan.snp"
    output:
        varscan_snp = "/scratch/eknodel/TSAI/varscan/{subject_sal}.saliva.varscan.snp.Somatic.hc"
    params:
        varscan = config["varscan_path"]
    shell:
        """
        java -jar {params.varscan} processSomatic {input.varscan_snp}
        """

rule somatic_filter_with_saliva:
    input:
        snp_somatic_hc = "/scratch/eknodel/TSAI/varscan/{subject_sal}.saliva.varscan.snp.Somatic.hc",
        indel = "/scratch/eknodel/TSAI/varscan/{subject_sal}.saliva.varscan.indel"
    output:
        snp_somatic_hc_filter = "/scratch/eknodel/TSAI/varscan/{subject_sal}.saliva.varscan.snp.Somatic.hc.filter"
    params:
        varscan = config["varscan_path"]
    shell:
        """
        java -jar {params.varscan} somaticFilter {input.snp_somatic_hc} -indel-file {input.indel} -output-file {output.snp_somatic_hc_filter}
        """

rule convert_to_vcf_with_saliva:
    input:
        "/scratch/eknodel/TSAI/varscan/{subject_sal}.saliva.varscan.snp.Somatic.hc.filter"
    output:
        "/scratch/eknodel/TSAI/varscan/{subject_sal}.saliva.varscan.snp.Somatic.hc.filter.vcf"
    shell:
        """
        cat {input} | awk '{{print $1"\t" $2"\t" "." "\t" $3 "\t" $4}}' > {output}
        """
