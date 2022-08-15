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
#bioc_packages = c("GenomicFeatures", "tximport", "rhdf5", "DESeq2")

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

# list of cran packages
#cran_packages = c("readr")
# load or install&load all
#package.check <- lapply(
#  cran_packages,
#  FUN = function(x) {
#    if (!require(x, character.only = TRUE)) {
#      install.packages(x, dependencies = TRUE, repos = "http://cran.us.r-project.org")
#      library(x, character.only = TRUE)
#    }
#  }
#)


#load libraries
library(GenomicFeatures)
library(tximport)
library(rhdf5)
library(readr)
library(DESeq2)


txdb_snakemake.wildcards <- makeTxDbFromGFF(snakemake@input[["annotation"]], format="auto")

k_snakemake.wildcards <- keys(txdb_snakemake.wildcards, keytype = "TXNAME")
tx2gene_snakemake.wildcards <- select(txdb_snakemake.wildcards, k_snakemake.wildcards, "GENEID", "TXNAME")

#user choice for whether gene or transcript level analysis should be performed
quantification_level <- as.character(snakemake@params[["quantification_level"]])
print(quantification_level)
quantification_level == "NO"

samples <- read.delim("config/samples.tsv")
samples_snakemake.wildcards <- samples[samples$species == snakemake@params[["species"]], ]
sample_snakemake.wildcards_names <- paste0(samples_snakemake.wildcards$sample, "_", samples_snakemake.wildcards$unit)

files_snakemake.wildcards <- file.path("kallisto_quant", snakemake@wildcards, sample_snakemake.wildcards_names, "abundance.h5")

names(files_snakemake.wildcards) <- paste0(sample_snakemake.wildcards_names)


if (quantification_level == "NO") {
  #if user chose gene level analysis
  txi.kallisto_snakemake.wildcards <- tximport(files_snakemake.wildcards, type = "kallisto", tx2gene = tx2gene_snakemake.wildcards, txOut = FALSE, importer = readr::read_tsv)
} else {
  #else just output/use transcript-level - txOut = TRUE
  txi.kallisto_snakemake.wildcards <- tximport(files_snakemake.wildcards, type = "kallisto", txOut = TRUE, importer = readr::read_tsv)
}

#sampleTable <- samples[1:8,"condition", drop=FALSE]
rownames(samples_snakemake.wildcards) <- colnames(txi.kallisto_snakemake.wildcards$counts)


DESeqDataSet_snakemake.wildcards <- DESeqDataSetFromTximport(txi.kallisto_snakemake.wildcards, samples_snakemake.wildcards, ~condition)

saveRDS(DESeqDataSet_snakemake.wildcards, file=snakemake@output[[1]])
