# Setting up filesnames here:
from os.path import join
import os

configfile: "SCC.config.json"

rule all:
    input:
        expand("variants/{subject}_ns_merged.vcf", subject=config["patients_with_ns"]), # combine all mutations
        expand("variants/{subject}_ns_merged_vep.vcf", subject=config["patients_with_ns"]), # run VEP
        expand("peptides/{subject}_ns_vep.17.vcf", subject=config["patients_with_ns"]), # run pvacseq
        expand("peptides/{subject}_ns_vep.17.peptides", subject=config["patients_with_ns"]),
        #expand("variants/{subject_sal}_saliva_merged.vcf", subject_sal=config["patients_with_saliva"]), # combine all mutations
        #expand("variants/{subject_sal}_saliva_merged_vep.vcf", subject_sal=config["patients_with_saliva"]), # run VEP
        #expand("peptides/{subject_sal}_saliva_vep.17.vcf", subject_sal=config["patients_with_saliva"]), # run pvacseq
        #expand("peptides/{subject_sal}_saliva_vep.17.peptides", subject_sal=config["patients_with_saliva"])

rule combine_strelka_ns:
    input: 
        strelka = os.path.join("/scratch/eknodel/TSAI/strelka/{subject}_ns/results/variants/somatic.snvs.pass.vcf.gz"),      
        strelka_indels = os.path.join("/scratch/eknodel/TSAI/strelka/{subject}_ns/results/variants/somatic.indels.pass.vcf.gz")              
    output:
        strelka = os.path.join("/scratch/eknodel/TSAI/strelka/{subject}_ns/results/variants/somatic.pass.vcf")
    shell:
        """
        bcftools concat -a {input.strelka_indels} {input.strelka} -o {output.strelka}
        """

rule sort_strelka:
      input:
          vcf = os.path.join("/scratch/eknodel/TSAI/strelka/{subject}_ns/results/variants/somatic.pass.vcf")
      output:
          vcf = os.path.join("/scratch/eknodel/TSAI/strelka/{subject}_ns/results/variants/somatic.pass_sorted.vcf")
      shell:
          """
          picard SortVcf I={input.vcf} O={output.vcf}
          """

rule gunzip_ns:
    input:
        gatk = os.path.join("/scratch/eknodel/TSAI/gatk_mutect2/{subject}.somatic.ns.filtered.pass.vcf.gz")
    output:
        gatk = os.path.join("/scratch/eknodel/TSAI/gatk_mutect2/{subject}.somatic.ns.filtered.pass.vcf")
    shell: 
        """
        gunzip {input.gatk}
        """

rule sort_gatk:
      input:
          vcf = os.path.join("/scratch/eknodel/TSAI/gatk_mutect2/{subject}.somatic.ns.filtered.pass.vcf")
      output:
          vcf = os.path.join("/scratch/eknodel/TSAI/gatk_mutect2/{subject}.somatic.ns.filtered.pass_sorted.vcf")
      shell:
          """
          picard SortVcf I={input.vcf} O={output.vcf}
          """

rule generate_input_ns:
    input:
        gatk = os.path.join("/scratch/eknodel/TSAI/gatk_mutect2/{subject}.somatic.ns.filtered.pass_sorted.vcf"),
        strelka = os.path.join("/scratch/eknodel/TSAI/strelka/{subject}_ns/results/variants/somatic.pass_sorted.vcf"),
    output:
        vep = os.path.join("variants/{subject}_ns_merged.vcf")
    shell: 
        """
        python merge_Mutect2_Strelka2_variants.py \
        --gatk_file {input.gatk} \
        --strelka_file {input.strelka} \
        --vep_format_fn {output.vep}
        """

rule run_vep_ns:
    input:
        os.path.join("variants/{subject}_ns_merged.vcf")
    output:
        os.path.join("variants/{subject}_ns_merged_vep.vcf")
    shell:
        """
        vep -i {input} --format vcf --assembly GRCh38 --cache --dir_cache ~/external_scripts --offline --vcf -o {output} --force_overwrite --plugin Wildtype --plugin Frameshift --symbol --terms SO --plugin Downstream
        #vep -i {input} --cache --dir_cache ~/external_scripts --gff /data/CEM/shared/public_data/references/T2T_CHM13_v2/CHM13.v2.0.gff.gz --fasta /data/CEM/shared/public_data/references/T2T_CHM13_v2/GCA_009914755.4_CHM13_T2T_v2.0_genomic.fna --vcf -o {output} --force_overwrite --plugin Wildtype --symbol --terms SO --plugin Downstream
        """ 

rule va_annotation:
    input:
        os.path.join("variants/{subject}_ns_merged_vep.vcf")
    output:
        os.path.join("variants/{subject}_ns_merged_vep_gt.vcf")    
    params:
        sample = "{subject}"
    shell:
        """
        vcf-genotype-annotator -o {output} {input} {params.sample} 0/1;
        """

rule generate_fasta_ns:
    input:
        os.path.join("variants/{subject}_ns_merged_vep_gt.vcf")
    output:
        len_17 = os.path.join("peptides/{subject}_ns_vep.17.vcf"),
    params: 
        sample = "{subject}"
    shell:
        """
        python /home/eknodel/bin/pVACtools/pvactools/tools/pvacseq/generate_protein_fasta.py {input} 8 {output.len_17} -s {params.sample};
        """

rule format_peptides:
    input:
        os.path.join("peptides/{subject}_ns_vep.17.vcf")
    output:
        os.path.join("peptides/{subject}_ns_vep.17.peptides")
    shell:
        """
        cat {input} | sed ':a;N;$!ba;s/\\n/ /g' | sed 's/>/\\n>/g' | sed 's/ //2' | sed 's/ //2' | sed 's/ //2' | sed 's/ //2' | sed 's/ //2' | sed 's/ //2' | sed '/>WT/d' | sed 's/>MT.*ENST00/>MT\.ENST00/g' | sed 's/\./ /g' > {output};
        """

