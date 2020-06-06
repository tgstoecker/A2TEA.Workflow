library(GenomicFeatures)

txdb_snakemake.wildcards <- makeTxDbFromGFF(snakemake@params[["annotation"]], format="auto")

k_snakemake.wildcards <- keys(txdb_snakemake.wildcards, keytype = "TXNAME")
tx2gene_snakemake.wildcards <- select(txdb_snakemake.wildcards, k_snakemake.wildcards, "GENEID", "TXNAME")


library(tximport)
library(rhdf5)
library(readr)


samples <- read.delim("samples.tsv")
samples_snakemake.wildcards <- samples[samples$species == snakemake@params[["species"]], ]
sample_snakemake.wildcards_names <- paste0(samples_snakemake.wildcards$sample, "_", samples_snakemake.wildcards$unit)

files_snakemake.wildcards <- file.path("kallisto_quant", snakemake@wildcards, sample_snakemake.wildcards_names, "abundance.h5")

names(files_snakemake.wildcards) <- paste0(sample_snakemake.wildcards_names)

txi.kallisto_snakemake.wildcards <- tximport(files_snakemake.wildcards, type = "kallisto", tx2gene = tx2gene_snakemake.wildcards, txOut = FALSE, importer = readr::read_tsv)



library(DESeq2)

#sampleTable <- samples[1:8,"treatment", drop=FALSE]
rownames(samples_snakemake.wildcards) <- colnames(txi.kallisto_snakemake.wildcards$counts)


DESeqDataSet_snakemake.wildcards <- DESeqDataSetFromTximport(txi.kallisto_snakemake.wildcards, samples_snakemake.wildcards, ~treatment)

saveRDS(DESeqDataSet_snakemake.wildcards, file=snakemake@output[[1]])
