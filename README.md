# Human SCC data

Analysis of mutational burden and signature from existing scc whole exome sequencing datasets.

----

## Update information

Last Uptated: 11/1/2023

Updated by: Elizabeth Borden 

Contact: knodele@arizona.edu

## Datasets utilized

| Dataset | Paper Link | Number of samples |
| ------- | --------- | ----------------- |
| Ji et al. 2020 Cell | [Full text](https://www.sciencedirect.com/science/article/pii/S0092867420306723) | 9 | 
| ------- | --------- | ----------------- |
| Chitsazzadeh et al. 2016 Nature communications | [Full text](https://www.nature.com/articles/ncomms12601) | 7 | 
| ------- | --------- | ----------------- |
| C. Zheng et al. 2014 Cell Reports | [Full text](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4254608/) | 13 |
| ------- | --------- | ----------------- |
| Q. Zheng et al. 2021 Journal of Investigative Dermatology | [Full text](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7790860/) | 10 |
| ------- | --------- | ----------------- |

## Environments

Environments for replicating code are available in the Setup directory as .yml files for conda environments

## Pipeline

1. Mutation calling

	- Variants identified with [GATK Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2) and [Strelka2](https://github.com/Illumina/strelka)
	- To ensure high fidelity variant calls, only variants called by both software are kept

2. Variant annotation

	- Variants are annotated with the [Ensembl variant effects predictor](https://useast.ensembl.org/info/docs/tools/vep/index.html)
	- Peptides generated with [PVACseq](https://pvactools.readthedocs.io/en/latest/pvacseq.html)
	- Variants prepared for visualization with [VCF2MAF](https://github.com/mskcc/vcf2maf)
