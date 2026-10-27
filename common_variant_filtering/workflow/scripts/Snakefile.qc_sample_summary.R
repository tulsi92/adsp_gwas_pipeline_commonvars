set.seed(19690329)

## Load packages
library(tidyverse)
library(magrittr)
library(SeqVarTools)
library(SNPRelate)
library(readr)

## Inputs
# sample_missingness_tsv <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/test/adsp.r4.wgs_cc.all.chrall.sample_missingness_filter.tsv"
# sample_hethom_tsv <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/test/adsp.r4.wgs_cc.all.chrall.sample_hethom_filter.tsv"
# sample_titv_tsv <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/test/adsp.r4.wgs_cc.all.chrall.sample_titv_filter.tsv"
# sample_inbreed_tsv <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/test/adsp.r4.wgs_cc.all.chrall.sample_inbreed_filter.tsv"
# sample_qc_summary_tsv <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/test/adsp.r4.wgs_cc.all.chrall.sample_qc_summary.tsv"
# refp.super_pops <- "all"
# snp_maf <- 0.1

sample_pheno_tsv <- snakemake@input$sample_pheno_tsv
sample_missingness_tsv <- snakemake@input$sample_missingness_tsv
sample_hethom_tsv <- snakemake@input$sample_hethom_tsv
sample_titv_tsv <- snakemake@input$sample_titv_tsv
sample_inbreed_tsv <- snakemake@input$sample_inbreed_tsv
sample_qc_summary_tsv <- snakemake@output$sample_qc_summary_tsv
refp.super_pops <- snakemake@params$refp_super_pops
snp_maf <- snakemake@params$snp_maf

## Script
`%notin%` <- Negate(`%in%`)

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

## Merge sample metrics files
sample_missingness <- read_tsv(sample_missingness_tsv)
sample_hethom <- read_tsv(sample_hethom_tsv)
sample_titv <- read_tsv(sample_titv_tsv)
sample_inbreed <- read_tsv(sample_inbreed_tsv)

sample_pheno_metrics <- sample_pheno |>
  dplyr::select(sample.id, sex, diagnosis, population) |>
  left_join(sample_missingness, by = "sample.id", suffix = c("", ".miss")) |>
  left_join(sample_hethom, by = "sample.id", suffix = c(".miss", ".hethom")) |>
  left_join(sample_titv, by = "sample.id", suffix = c(".hethom", ".titv")) |>
  left_join(sample_inbreed, by = "sample.id", suffix = c(".titv", ".inbreed"))

head(sample_pheno_metrics)

## Check pass samples
sample_pheno_metrics |> 
  dplyr::filter(missing_rate < snp_maf) |>
  dplyr::count(pass_hethom, pass_titv, pass_inbreed_coeff)

## Write output
sample_pheno_metrics |> 
  write_tsv(sample_qc_summary_tsv)

