#logging
log <- file(snakemake@log[[1]], open="wt")
sink(log, append=TRUE)
sink(log, append=TRUE, type="message")

# Restore output to console
#sink() 
#sink(type="message")

#backup R-based installation if conda didn't work or wasn't used
#we check if packages are installed first

# list of bioconductor packages
bioc_packages = c("DESeq2")

# load or install&load all
package.check <- lapply(
  bioc_packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
    BiocManager::install(x)
    library(x, character.only = TRUE)
    }
  }
)


#load libraries
library("DESeq2")

dds_snakemake.wildcards <- readRDS(snakemake@input[[1]])


# if normal DESeq function fails due to lack of dispersion use alternative way:
# "all gene-wise dispersion estimates are within 2 orders of magnitude"

# "<<-" global assignment

tryCatch(
    {
        dea_snakemake.wildcards <<- DESeq(dds_snakemake.wildcards)
    },
    error=function(b) {
                         dea_snakemake.wildcards <- estimateSizeFactors(dds_snakemake.wildcards)
                         dea_snakemake.wildcards <- estimateDispersionsGeneEst(dea_snakemake.wildcards)
                         dispersions(dea_snakemake.wildcards) <- mcols(dea_snakemake.wildcards)$dispGeneEst
                         dea_snakemake.wildcards <<- nbinomWaldTest(dea_snakemake.wildcards)
                }
)

saveRDS(dea_snakemake.wildcards, file=snakemake@output[[1]])

