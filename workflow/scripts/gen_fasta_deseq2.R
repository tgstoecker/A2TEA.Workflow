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
#bioc_packages = c("DESeq2")

# load or install&load all
#package.check <- lapply(
#  bioc_packages,
#  FUN = function(x) {
#    if (!require(x, character.only = TRUE)) {
#    BiocManager::install(x)
#    library(x, character.only = TRUE)
#    }
#  }
#)


# load libraries
library(DESeq2)

gen_fasta_samples <- read.delim("R/gen_fasta_samples.csv")

gfs_snakemake.wildcards <- gen_fasta_samples[gen_fasta_samples$species == snakemake@params[["species"]], ]
gfs_snakemake.wildcards_names <- paste0(gfs_snakemake.wildcards$sample, "_", gfs_snakemake.wildcards$unit)


counts=read.csv(snakemake@input[[1]], sep = "\t", header = TRUE, skip = 0, row.names = "Geneid")
counts <- counts[, c(-1:-5)]
colnames(counts) <- gfs_snakemake.wildcards_names


#object for storing names of samples and condition of the samples 
gfs_snakemake.wildcards$samples <- paste0(gfs_snakemake.wildcards$sample, sep="_", gfs_snakemake.wildcards$unit)
gfs_snakemake.wildcards <- gfs_snakemake.wildcards[, c("samples", "condition")]
rownames(gfs_snakemake.wildcards) = gfs_snakemake.wildcards[,1]


## Check if sample names created match with the column names of featurecounts. both ## ## should match 
all(rownames(gfs_snakemake.wildcards) %in% colnames(counts))

dds = DESeqDataSetFromMatrix(countData = counts, colData = gfs_snakemake.wildcards, design = ~ condition)
#dea = DESeq(dds)

# if normal DESeq function fails due to lack of dispersion use alternative way:
# "all gene-wise dispersion estimates are within 2 orders of magnitude"
# "<<-" global assignment

tryCatch(
    {
        dea <<- DESeq(dds)
    },
    error=function(b) {
                         dea <- estimateSizeFactors(dds)
                         dea <- estimateDispersionsGeneEst(dea)
                         dispersions(dea) <- mcols(dea)$dispGeneEst
                         dea <<- nbinomWaldTest(dea)
                }
)


saveRDS(dea, file=snakemake@output[[1]])
