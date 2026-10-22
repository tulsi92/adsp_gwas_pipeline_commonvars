
#!/usr/bin/Rscript

##load packages and set variables

library(readr)

library(dplyr)

library(stringr)

library(SeqArray)

args <- commandArgs(trailingOnly = TRUE)

adsp_meta_name <- args[1]

adsp_variant_qc_name <- args[2]

adsp_gds_name <- args[3]

adsp_variant_qc_het_hwe_inbreed_name <- args[4]

adsp_vcf_name <- args[5]

hwe_p_thresh <- as.numeric(args[6])



##import adsp metadata

adsp_meta <- read_tsv(adsp_meta_name, guess_max = 1e7)



##import and process adsp variant qc data

adsp_variant_qc <- read_tsv(adsp_variant_qc_name, guess_max = 1e9) %>%

    mutate(pass_hwe = ifelse(hwe_p < hwe_p_thresh, FALSE, TRUE))

adsp_variant_qc_het <- adsp_variant_qc %>%

    filter(pass_heterozygosity == FALSE) %>%

    mutate(why = "het") %>%

    select(variant.id, why)

adsp_variant_qc_hwe <- adsp_variant_qc %>%

    filter(pass_hwe == FALSE) %>%

    mutate(why = "hwe") %>%

    select(variant.id, why)

adsp_variant_qc_inbreed <- adsp_variant_qc %>%

    filter(pass_inbreed == FALSE) %>%

    mutate(why = "inbreed") %>%

    select(variant.id, why)

adsp_variant_qc_het_hwe_inbreed <- bind_rows(adsp_variant_qc_het, adsp_variant_qc_hwe, adsp_variant_qc_inbreed) %>%

    arrange(variant.id)

adsp_variant_qc_pass <- anti_join(adsp_variant_qc, adsp_variant_qc_het_hwe_inbreed)



##output adsp variant qc data

write_tsv(adsp_variant_qc_het_hwe_inbreed, adsp_variant_qc_het_hwe_inbreed_name)



##import adsp gds data then set sample and variant filters

adsp_gds <- seqOpen(adsp_gds_name)

seqSetFilter(adsp_gds, sample.id=adsp_meta$IID)

seqSetFilter(adsp_gds, variant.id=adsp_variant_qc_pass$variant.id)



##output adsp filtered vcf data

seqGDS2VCF(adsp_gds, adsp_vcf_name)

