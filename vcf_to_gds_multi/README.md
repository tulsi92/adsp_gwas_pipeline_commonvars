# ADSP R4 Preprocessing

This module preprocesses ADSP Release 4 GCAD whole-genome sequencing (WGS) data for downstream genome-wide association analyses. The workflow performs sample- and variant-level filtering, prepares genotype files for PCA, annotates variants with gnomAD frequencies, and generates GDS files for downstream analyses.

## Workflow overview

<img width="2152" height="2540" alt="image" src="https://github.com/user-attachments/assets/797cba4c-770a-4811-bb77-c93d0202ebc2" />

The preprocessing workflow prepares ADSP Release 4 genotype data for downstream GWAS analyses and:

1. Decomposes and filters multiallelic variants
2. Combines the filtered multiallelic variants with previously processed biallelic variants
3. Sorts, indexes, and repairs the merged chromosome-level VCFs
4. Annotates variants with gnomAD allele frequencies and functional scores
5. Creates genome-wide variant annotation tables
6. Converts the processed VCFs to GDS format for full-data and PCA analyses

The workflow is run across autosomes 1–22. Paths and parameters are defined in `config_files/config.tsv`, while shared parameters such as the PCA MAF threshold and gnomAD version are defined in the Snakefile.

The major steps and corresponding rules are described in the table below:

| Major Step | Snakemake Rules | Purpose | Key outputs |
| --- | --- | --- | --- |
| **1. Decompose multiallelic variants** | `decompose_multiallelic_vcf` | Decompose multiallelic variants | Filtered chromosome-level multiallelic VCFs |
| **2. Multiallelic quality control and filtering** | `process_chr_vcf_multiallelic`, `extract_quality_flags`, `filter_quality_flags`, `filter_vcf_abhet`, `filter_vcf_for_maf_pca` | Extract variant quality metrics (e.g. VFLAGS and ABHet), apply initial quality filters, identify high-quality multiallelic variants and remove variants failing GCAD QC | Filtered chromosome-level decomposed VCFs |
| **3. Variant annotations** | `annotate_vcf_with_gnomad`, `extract_chr_variant_annot`, `merge_chr_variant_annot`, `extract_chr_variant_annot_gnomad`, `merge_chr_variant_annot_gnomad` | Annotate processed variants with gnomAD population allele frequencies and functional prediction scores (CADD, REVEL, SpliceAI, Pangolin, SIFT, PolyPhen, etc.). Extract functional and population frequency annotations from the annotated VCFs and merge chromosome-level files into genome-wide annotation tables | gnomAD-annotated chromosome VCFs and variant annotation tables |
| **4. Generate final processed VCFs** | `concat_biallelic_and_multiallelic_chr_vcf`, `sort_concat_chr_vcf`, `index_vcf_sorted`, `reheader_vcf_to_fix` | Merge processed biallelic and multiallelic variants, sort and index chromosome VCFs, and correct the ABHet header before downstream analyses | Final processed chromosome VCFs |
| **5. Generate GDS files** | `convert_chr_vcfs_to_gds`, `convert_chr_vcfs_to_gds_for_pca` | Convert processed VCFs into GDS format for downstream GWAS analyses and generate a separate common-variant GDS optimized for PCA | Full GDS and PCA GDS files |

## Prerequisites

### Required inputs

- ADSP Release 4 chromosome-level VCFs
- GCAD filtering recommendations
- gnomAD annotation VCFs

### Software

The workflow is implemented in Snakemake and uses:

- bcftools
- R
- SeqArray
- SeqVarTools
- tabix
- bgzip

Required software and reference files:

- SnpSift.jar
- snpEff.jar
- snpEff.config
- GRCh38.86 reference database

Container and Conda environments are provided where applicable.

---

## Configuration

Before running the pipeline:

- Update `config_files/config.tsv` with the appropriate input paths and cohort-specific settings
- Modify any remaining hard-coded paths in the workflow to match your local environment

---

## Running the pipeline

- Perform a dry run first: `snakemake -np`
- Run the entire workflow: `snakemake --profile lsf`
- Run until a specific rule: `snakemake --profile lsf --until <rule_name>`

---

## Outputs

The preprocessing workflow generates:

- Full genome GDS file
- PCA-ready GDS file
- Genome-wide variant annotation tables

These outputs are used as inputs for the `gwas/` module.

---

## Expected runtime

Processing is computationally intensive and scales with cohort size.

Typical runtimes are:

- **Chromosome preprocessing:** hours to days per chromosome
- **Variant annotation:** several hours per chromosome
- **Complete preprocessing:** several weeks

Running chromosomes in parallel is strongly recommended.
