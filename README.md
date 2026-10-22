# ADSP R4 GWAS Pipeline (Common variants)

This repository contains the code required to perform genome-wide association studies (GWAS) using Alzheimer's Disease Sequencing Project (ADSP) Release 4 36k WGS multiallelic data. The pipeline is organized into independent analysis modules covering genotype preprocessing and downstream association analyses. GCAD processed ADSP gVCF files were downloaded from the NIAGADS AWS.

### Repository structure
Each analysis module is self-contained. Read the folder-specific README.md first to follow the documented execution order and required inputs.

| Folder | Analysis | Key inputs	| Key outputs |
|--------|----------|------------|-------------|
| `phen/` | Phenotype processing	| ADSP phenotype files and documentation |	Analysis-ready phenotype and covariates files |
| `vcf_to_gds_multi/` | ADSP genotype preprocessing and GCAD QC	| Raw ADSP genotype data, phenotype files, GCAD filtering recommendations |	Filtered genotype data, QC metrics, variant annotation files |
| `qc_and_pca/` | ADSP genotype QC	| Filtered genotype data for all common variants |	Post-QC genotype data and PCs |
| `slope_prs/` |	x	|	x | x |

Each module contains:
* A detailed `README.md`
* Snakemake workflow (`Snakefile`)
* Required software environments (container and/or conda)
* Analysis-specific configuration options

### Data
ADSP genotype data and phenotype files are **not** included in this repository due to data access restrictions. Users must have authorized access to the required ADSP datasets and associated phenotype/covariate files.

### Configuration
Some input file paths may be hard-coded within the workflow and should be updated to match your local environment before running the pipeline. Analysis parameters, file names, and other user-configurable settings are defined in the CONFIG section of each script (where applicable) and can be modified as needed.
