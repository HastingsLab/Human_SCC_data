# Setting up filenames here:
from os.path import join
# Config_file
configfile: "SCC.config.json"

# Files
peptide_path = "peptides/"

rule all:
    input:
        expand(peptide_path + "{sample}_9_netmhc_ns.xsl", sample=config["patients_with_ns"]),
        expand(peptide_path + "{sample}_9_netmhcstab_ns.xsl", sample=config["patients_with_ns"]),
        #expand(peptide_path + "{sample}_9_netmhc_saliva.xsl", sample=config["patients_with_saliva"]),
        #expand(peptide_path + "{sample}_9_netmhcstab_saliva.xsl", sample=config["patients_with_saliva"])


rule prepare_input:
    input:
        peptides = os.path.join("peptides/{sample}_ns_vep.17.peptides"),
    output:
        peptides = os.path.join("peptides/{sample}_ns_vep.17.peptides_formatted"),
    shell:
        """
        perl -lane 'print "@F[-2..-1]"' {input.peptides} | sed 's/^/>/' | sed 's/ /\\n/g' > {output.peptides};
        """

rule netMHC_ns:
    input:
        peptides = os.path.join("peptides/{sample}_ns_vep.17.peptides_formatted")
    output:
        netMHC_9  = os.path.join(peptide_path, "{sample}_9_netmhc_ns.xsl"),
    params:
        hla = lambda wildcards: config[wildcards.sample]["hla"]
    shell:
        """
        netMHCpan -a {params.hla} -f {input.peptides} -BA -s -xls -l 9  -xlsfile {output.netMHC_9};
        """

rule netMHCstab_ns:
    input:
        peptides = os.path.join("peptides/{sample}_ns_vep.17.peptides_formatted")
    output:
        netMHCstab_9  = os.path.join(peptide_path, "{sample}_9_netmhcstab_ns.xsl"),
    params:
        hla = lambda wildcards: config[wildcards.sample]["hla"]
    shell:
        """
        netMHCstabpan -a {params.hla} -f {input.peptides} -s -xls -l 9  -xlsfile {output.netMHCstab_9};
        """

rule prepare_input_saliva:
    input:
        peptides = os.path.join("peptides/{sample}_saliva_vep.17.peptides"),
    output:
        peptides = os.path.join("peptides/{sample}_saliva_vep.17.peptides_formatted"),
    shell:
        """
        cat {input.peptides} | awk '{{print $1"."$7, $8}}' | sed 's/>MT/>/' | sed 's/ /\\n/g' > {output.peptides};
        """

rule netMHC_saliva:
    input:
        peptides = os.path.join("peptides/{sample}_saliva_vep.17.peptides_formatted")
    output:
        netMHC_9  = os.path.join(peptide_path, "{sample}_9_netmhc_saliva.xsl"),
    params:
        hla = lambda wildcards: config[wildcards.sample]["hla"]
    shell:
        """
        netMHCpan -a {params.hla} -f {input.peptides} -BA -s -xls -l 9  -xlsfile {output.netMHC_9};
        """

rule netMHCstab_saliva:
    input:
        peptides = os.path.join("peptides/{sample}_saliva_vep.17.peptides_formatted")
    output:
        netMHCstab_9  = os.path.join(peptide_path, "{sample}_9_netmhcstab_saliva.xsl"),
    params:
        hla = lambda wildcards: config[wildcards.sample]["hla"]
    shell:
        """
        netMHCstabpan -a {params.hla} -f {input.peptides} -s -xls -l 9  -xlsfile {output.netMHCstab_9};
        """
