### 2022-01-19
### Rscript to run Ruth for ancestry-aware HWE on ADSP WGS/WES
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(SeqVarTools))
suppressPackageStartupMessages(library(SeqArray))
suppressPackageStartupMessages(library(SNPRelate))

# seqParallelSetup(cluster=8)

## PCA TIBBLE FUN
PCA_tibble <- function(pc, iid) {
  pc = as.data.frame(pc)
  colnames(pc) <- paste0("PC", 1:ncol(pc))
  as_tibble(pc) %>%
    mutate(IID = as.character(iid)) %>%
    select(IID, everything())
}

## FILES
pruned_merged.vcf <- "/hpc/users/patelt16/projects/adsp/ruth/050922/output/adsp-r2-wes-test.chrall.merged.vcf.gz"
outfn.gds         <- "/hpc/users/patelt16/projects/adsp/ruth/050922/output/chrall.pruned.gds"
pca_plot_evec     <- "/hpc/users/patelt16/projects/adsp/ruth/050922/output/adsp_pca_plots.pdf"
num_pcs           <- 3
pca_algorithm     <- "randomized"
outfn.eval        <- "/hpc/users/patelt16/projects/adsp/ruth/050922/output/adsp-r2-wes-test.eval.txt"
outfn.evec        <- "/hpc/users/patelt16/projects/adsp/ruth/050922/output/adsp-r2-wes-test.evec.txt"

pruned_merged.vcf <- snakemake@input[["pruned_merged_vcf"]]
outfn.gds         <- snakemake@output[["pruned_merged_gds"]]
# pca_plot_evec     <- snakemake@output[["pca_plots_pdf"]]
num_pcs           <- snakemake@params[["num_pcs"]]
pca_algorithm     <- snakemake@params[["pca_algorithm"]]
outfn.eval        <- snakemake@output[["eigenval"]]
outfn.evec        <- snakemake@output[["eigenvec"]]

print(num_pcs)

## CONVERT VCF TO GDS
pruned.gds        <- snpgdsVCF2GDS(pruned_merged.vcf, outfn.gds, verbose = TRUE)

## RUN PCA
pruned.gds.fn     <- snpgdsOpen(pruned.gds)
pruned.gds.pca    <- snpgdsPCA(pruned.gds.fn, eigen.cnt = as.numeric(num_pcs), algorithm = pca_algorithm, verbose = TRUE, num.thread=24) #eigen.cnt=0 to get all EVs

## add config option for algorithm and output rds file
df.eval           <- data.frame(eigenval = pruned.gds.pca$eigenval)
write_tsv(df.eval, outfn.eval)

# head(round(pruned.gds.pca$varprop))
# head(round(pruned.gds.pca$varprop*100, 2)) # variance proportion %
# write_tsv(df.varprop, outfn.eval)

df.evec <- PCA_tibble(pruned.gds.pca$eigenvect, pruned.gds.pca$sample.id)
write_tsv(df.evec, outfn.evec)

pcplot <- df.evec %>%
  ggplot(aes(x = PC1, y = PC2)) +
    geom_point() + theme_bw() + labs(x = "PC1", y = "PC2")

# pdf(pca_plot_evec)
# plot(pruned.gds.pca)
# plot(pruned.gds.pca, 1:3)
# dev.off()

# print(pcplot)
# ggsave(pca_plot_evec, plot = pcplot, device="pdf")

## CLOSE GDS
# snpgdsClose(pruned.gds.pca)
showfile.gds(closeall = TRUE)
