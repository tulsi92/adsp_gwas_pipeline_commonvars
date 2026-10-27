library(tidyverse)
library(SeqArray)

gds_fn <- snakemake@input$gds
sample_pheno_tsv <- snakemake@input$sample_pheno_tsv
chr_vcf <- snakemake@output$chr_vcf
# snp_maf <- snakemake@params$snp_maf
# missing.rate <- snakemake@params$missing_rate

chr <- snakemake@wildcards$chr

sample_pheno <- read_tsv(sample_pheno_tsv)

gds <- seqOpen(gds_fn)

seqSetFilter(gds, sample.id = sample_pheno$sample.id)

# seqSetFilterCond(gds, snp_maf = c(0, maf), missing.rate = missing.rate)

seqSetFilterChrom(gds, chr)

seqGDS2VCF(gds, chr_vcf)

gds

seqClose(gds)

sessionInfo()
