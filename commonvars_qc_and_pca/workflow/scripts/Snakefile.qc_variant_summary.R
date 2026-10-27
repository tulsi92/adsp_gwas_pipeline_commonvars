set.seed(19690329)

## Load packages
library(tidyverse)
library(magrittr)
library(readr)

## Inputs
# variant_heterozygosity_tsv <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/test/adsp.r4.wgs_cc.all.chrall.variant_heterozygosity.tsv"
# variant_hwe_tsv <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/test/adsp.r4.wgs_cc.all.chrall.variant_hwe.tsv"
# variant_diff_missingness_tsv <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/test/adsp.r4.wgs_cc.all.chrall.variant_differential_missingness.tsv"
# variant_qc_summary_tsv <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/test/adsp.r4.wgs_cc.all.chrall.variant_qc_summary.tsv"
# snp_maf <- 0.1

variant_heterozygosity_tsv <- snakemake@input$variant_heterozygosity_tsv
variant_hwe_qq_tsv <- snakemake@input$variant_hwe_qq_tsv
variant_hwe_tsv <- snakemake@input$variant_hwe_tsv
variant_diff_missingness_tsv <- snakemake@input$variant_diff_missingness_tsv
variant_qc_report_heterozygosity_pdf <- snakemake@output$qc_report_heterozygosity_pdf
variant_qc_report_hwe_pdf <- snakemake@output$qc_report_hwe_pdf
variant_qc_summary_tsv <- snakemake@output$variant_qc_summary_tsv
snp_maf <- snakemake@params$snp_maf

## Script
`%notin%` <- Negate(`%in%`)

## Variant metrics files
variant_heterozygosity <- read_tsv(variant_heterozygosity_tsv)
variant_hwe_qq <- read_tsv(variant_hwe_qq_tsv)
variant_hwe <- read_tsv(variant_hwe_tsv)
variant_diff_missingness <- read_tsv(variant_diff_missingness_tsv)

## Plot variant heterozygosity
pdf(variant_qc_report_heterozygosity_pdf, width = 10, height = 7.5)
ggplot(variant_heterozygosity, aes(x = heterozygosity)) +
  geom_histogram(bins = 100, fill = "#31a354", color = "white") +
  geom_vline(xintercept = 0.6, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Variant-level Heterozygosity",
       x = "Heterozygosity",
       y = "Number of Variants") +
  theme_bw(base_size = 14)

dev.off()

## Plot variant hwe QQ plot
pdf(variant_qc_report_hwe_pdf, width = 7, height = 7)
ggplot(variant_hwe_qq, aes(x = expected, y = observed)) +
  geom_point(alpha = 0.4, color = "#2c7fb8") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "HWE P-value QQ Plot",
       x = expression(-log[10]("Expected P")),
       y = expression(-log[10]("Observed P"))) +
  theme_bw(base_size = 14)

## Plot variant hwe
ggplot(variant_hwe, aes(x = inbreed_coeff)) +
  geom_histogram(bins = 100, fill = "#31a354", color = "white") +
  geom_vline(xintercept = c(-3,3)*sd(variant_hwe$inbreed_coeff)+mean(variant_hwe$inbreed_coeff),
             color = "red", linetype = "dashed") +
  labs(title = "Variant-level Inbreeding Coefficient",
       x = "Inbreeding Coefficient (F)",
       y = "Number of Variants") +
  theme_bw(base_size = 14)

dev.off()

## Combine variant metrics
variant_metrics <- variant_heterozygosity |>
  left_join(variant_hwe, by = "variant.id", suffix = c(".het", ".hwe")) |>
  left_join(variant_diff_missingness, by = "variant.id", suffix = c(".hwe", ".diffmiss"))

head(variant_metrics)

## Check pass variants
variant_metrics |>
  # dplyr::filter(missing_rate < snp_maf) |>
  dplyr::count(pass_heterozygosity, pass_inbreed, DIFF_MISS_EUR < 0.05 | DIFF_MISS_AFR < 0.05 | DIFF_MISS_AMR < 0.05)

## Write output
variant_metrics |> 
  write_tsv(variant_qc_summary_tsv)

