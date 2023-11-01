# Setting up filesnames here:
from os.path import join
import os

configfile: "SCC.config.json"

# Directory Pathways
peptide_path = "/home/eknodel/Tsai_newpipeline/02_variant_annotation/peptides/"

rule all:
    input:
        expand(os.path.join(peptide_path, "{sample}.peptide_ns.wildtype"), sample=config["patients_with_ns"]),
        expand(os.path.join("NeoScore/{sample}.neoscore_ns.csv"), sample=config["patients_with_ns"])

rule prepare_input:
    input:
        peptides = os.path.join(peptide_path, "{sample}_ns_vep.17.peptide")
    output:
        peptides = os.path.join(peptide_path, "{sample}.peptide_ns.wildtype")
    shell:
        """
        paste -s -d' \n' {input.peptides} | sed '/>MT/ d' > {output.peptides}
        """

rule prepare_transcript_names:
    input: 
        peptides = os.path.join(peptide_path, "{sample}_ns_vep.17.peptide")
    output:
        peptides = os.path.join(peptide_path, "{sample}_ns_transcript_names.out")
    shell: 
        """
        paste -s -d' \n' {input.peptides} | sed 's/^.*_ENST0/ENST0/g' | sed 's/_.* / /' > {output.peptides}
        """

rule NeoScore:
    input:
        netMHCpan = os.path.join(peptide_path, "{sample}_9_netmhc_polysolver_ns.xsl"),
        netMHCstab = os.path.join(peptide_path, "{sample}_9_netmhcstab_polysolver_ns.xsl"),
        quants_file = lambda wildcards: os.path.join("/scratch/eknodel/TSAI/salmon/", config[wildcards.sample]["rna_tumor"], "quant.sf"),
        wildtype = os.path.join(peptide_path, "{sample}.peptide_ns.wildtype"),
        transcripts = os.path.join(peptide_path, "{sample}_ns_transcript_names.out")
    output:
        output1 = os.path.join("NeoScore/{sample}.neoscore_ns.csv"),
        output2 = os.path.join("NeoScore/{sample}.neoscore_restricted_ns.csv"),
        output3 = os.path.join("NeoScore/{sample}.binding_neoantigens_ns.csv")
    shell:
        """
        Rscript Neoscore.R {input.netMHCpan} {input.netMHCstab} {input.quants_file} {input.wildtype} {input.transcripts} {output.output1} {output.output2} {output.output3}
        """ 
