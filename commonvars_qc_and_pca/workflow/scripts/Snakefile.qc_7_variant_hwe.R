set.seed(19690329)

## Load packages
library(tidyverse)
library(magrittr)
library(SeqVarTools)
library(SNPRelate)
library(readr)

## Inputs
# gds_fn <- "/sc/arion/projects/load/users/patelt16/projects/adsp/adsp_gwas/20240118.adsp_pipeline_v8_master_pipeline/multiallelic_gwas/vcf_to_gds_multi/output/adsp.r4.wgs/adsp.r4.wgs.biallelic_decomposed_multiallelic.chrall.gds"
# sample_pheno_tsv <- "/sc/arion/projects/load/users/patelt16/projects/adsp/adsp_gwas/20240118.adsp_pipeline_v8_master_pipeline/phen/20250519.multiallelic_run/20250519.ADSP_Umbrella_ng00067.v10.wgs.cc.all.28063.tsv"
# variant_hwe_tsv <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/adsp.r4.wgs_cc.all.variant_hwe.tsv"
# qc_summary_tsv <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/qc_metrics.tsv"
# qc_report_pdf <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/qc_report.pdf"
# snp_maf <- 0.01
# missing.rate <- 0.1
# refp.super_pops <- "all"

gds_fn <- snakemake@input$gds
sample_pheno_tsv <- snakemake@input$sample_pheno_tsv
variant_missingness_tsv <- snakemake@input$variant_missingness_tsv
variant_hwe_qq_tsv <- snakemake@output$variant_hwe_qq_tsv
variant_hwe_tsv <- snakemake@output$variant_hwe_tsv
# qc_summary_tsv <- snakemake@output$qc_summary_tsv
# qc_report_pdf <- snakemake@output$qc_report_pdf
snp_maf <- snakemake@params$snp_maf
missing.rate <- snakemake@params$missing_rate
refp.super_pops <- snakemake@params$refp_super_pops
chr <- snakemake@wildcards$chr

## Script
`%notin%` <- Negate(`%in%`)

gds <- seqOpen(gds_fn)

# seqResetFilter(gds)

## Sample pheno file
sample_pheno <- read_tsv(sample_pheno_tsv)

if (length(refp.super_pops) > 0 & 'super_pop' %notin% colnames(sample_pheno)) {
  sample_pheno <- sample_pheno |> mutate(super_pop = case_when(population == "NHW" ~ "EUR",
                                                               population == "AA" ~ "AFR",
                                                               population == "HI" ~ "AMR",
                                                               is.na(population) ~ NA,
                                                               TRUE ~ "OTHER"))
}

if (length(refp.super_pops) > 0 & refp.super_pops != "all" & 'super_pop' %in% colnames(sample_pheno)) {
  refp.super_pops_all <- unlist(strsplit(refp.super_pops, split = "_"))
  sample_pheno <- sample_pheno |>
    dplyr::filter(super_pop %in% c(refp.super_pops_all))
}

## Filter for samples in pheno file and biallelic variants
message("Filter for samples in pheno file and biallelic variants only")
seqSetFilter(gds, sample.id = sample_pheno$sample.id, variant.sel = isSNV(gds, biallelic=TRUE))

# seqSetFilterChrom(gds, include = chr, verbose = TRUE)

## Filter for variant missing rate
variant_metrics_filter <- read_tsv(variant_missingness_tsv) |>
  separate(variant.id, into = c("CHROM", "POS", "REF", "ALT"), sep = ":", remove = FALSE) |>
  dplyr::filter(CHROM == chr, !is.na(pass_missingness))

seqSetFilter(gds, variant.id = variant_metrics_filter$variant.id)

## Compute HWE p-values
message("Variant-level Hardy-Weinberg Equilibrium (HWE)/Inbreeding Coefficient Filtering")
variant.hwe <- hwe(gds)
variant.hwe.perm <- hwe(gds, permute = TRUE)

observed <- -log10(sort(variant.hwe$p))
expected <- -log10(sort(variant.hwe.perm$p))

## Create QQ plot data
qq.df <- tibble(
  expected = expected,
  observed = observed
)

## Save inbreeding coefficient metrics
variant.hwe.metrics <- tibble(
  variant.id = seqGetData(gds, "variant.id"),
  hwe_p = variant.hwe$p,
  inbreed_coeff = variant.hwe$f,
  pass_inbreed = abs(scale(variant.hwe$f)[,1]) <= 3
)

write_tsv(qq.df, variant_hwe_qq_tsv)
write_tsv(variant.hwe.metrics, variant_hwe_tsv)

seqClose(gds)
# showfile.gds(closeall=TRUE)

sessionInfo()