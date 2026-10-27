## Load packages
library(tidyverse)
library(magrittr)
library(SeqArray)

# gds <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/common_all/adsp.r4.wgs_cc.all.chrall.filtered.gds"
# sample_pheno_tsv <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/common_all/adsp.r4.wgs_cc.all.sample_annot.tsv"
# sample_qc_summary_tsv <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/common_all/report/adsp.r4.wgs_cc.all.chrall.sample_qc_summary.tsv"
# wes_capture_kit_info_tsv <- "/sc/arion/projects/load/data-ext/ADSP/ADSP_Umbrella_ng00067.v7/misc/gcad.r2.wes.20504.VCPA1.1.2019.11.01.v3.qcmetrics.seq.ALL.txt"
# sample_pheno_filtered_merged_tsv <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/common_all/adsp.r4.wgs_cc.all.chrall.sample_pheno_merged_covars.tsv"

gds_fn <- snakemake@input$gds
sample_pheno_tsv <- snakemake@input$sample_pheno_tsv
sample_qc_summary_tsv <- snakemake@input$sample_qc_summary_tsv
wes_capture_kit_info_tsv <- snakemake@input$wes_capture_kit_info_tsv
sample_pheno_filtered_merged_tsv <- snakemake@output$sample_pheno_filtered_merged_tsv

sample_qc_summary <- read_tsv(sample_qc_summary_tsv)

sample_filterlist <- sample_qc_summary |>
  dplyr::filter(pass_hethom == TRUE &
                pass_titv == TRUE &
                pass_inbreed_coeff == TRUE) %$% sample.id

wes_capture_kit_info <- read_tsv(wes_capture_kit_info_tsv) |>
  dplyr::select(sample.id = Sample_name, wes_target_regions = Target_regions, wes_pipeline_type = Pipeline_type)

sample_pheno <- read_tsv(sample_pheno_tsv) |>
  dplyr::filter(sample.id %in% sample_filterlist) |>
  left_join(wes_capture_kit_info, by = "sample.id") |>
  mutate(sequencing_center_platform = paste(sequencing_center, sequencing_platform, sep = "_"))

## Filter post-PCA sample list for samples passing QC as discrepancies
gds <- seqOpen(gds_fn)

seqSetFilter(gds, sample.id = sample_pheno$sample.id)

gds_samples <- seqGetData(gds, "sample.id")

sample_pheno_filtered_merged <- sample_pheno |>
  dplyr::filter(sample.id %in% gds_samples)

sample_pheno_filtered_merged |>
  write_tsv(sample_pheno_filtered_merged_tsv)

seqClose(gds)
