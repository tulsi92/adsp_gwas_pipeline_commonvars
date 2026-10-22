# ADSP R4 WGS Processing Pipeline

This repository contains the code required to perform genome-wide association studies (GWAS) using Alzheimer's Disease Sequencing Project (ADSP) Release 4 WGS multiallelic data. The pipeline is organized into independent analysis modules covering genotype preprocessing and downstream association analyses. GCAD processed ADSP gVCF and phenotype files were downloaded from the NIAGADS AWS.

### Repository structure
Each analysis module is self-contained. Read the folder-specific README.md first to follow the documented execution order and required inputs.

| Folder | Analysis | Key inputs	| Key outputs |
|--------|----------|------------|-------------|
| `vcf_to_gds_multi/` | ADSP genotype preprocessing and GCAD QC	| Raw ADSP genotype data, phenotype files, GCAD filtering recommendations |	GCAD filtered genotype GDS, variant annotation files |
| `common_variant_filtering/` | GCAD filtered genotype GDS	| GCAD filtered genotype GDS |	GCAD filtered genotype GDS for all common variants |
| `qc_and_pca/` | ADSP genotype conversion, QC and PCA	| GCAD filtered genotype GDS for all common variants |	Post-QC genotype plink files and PCs |

### Prepare final dataset for downstream analysis
Once the steps are complete, we remove ADNI samples and retain only AD cases for downstream analyses as shown below:

```
module load plink

ADSP_CC_PLINK="adsp.r4.wgs_cc.{date}"
SAMPLES_KEEP="adsp_cases_excl_adni_all.10236.pheno.tsv"
PLINK_OUT="adsp.r4.wgs.cases_excl_adni"

plink \
  --bfile "${ADSP_CC_PLINK}" \
  --keep "${SAMPLES_KEEP}" \
  --make-bed \
  --out "${PLINK_OUT}"
```

Each module contains:
* A detailed `README.md`
* Snakemake workflow (`Snakefile`)
* Required software environments (container and/or conda)
* Configuration file

### Data
ADSP genotype data and phenotype files are **not** included in this repository due to data access restrictions. Users must have authorized access to the required ADSP datasets and associated phenotype/covariate files.

### Configuration
Some input file paths may be hard-coded within the workflow and should be updated to match your local environment before running the pipeline. Analysis parameters, file names, and other user-configurable settings are defined in the CONFIG section of each script (where applicable) and can be modified as needed.
