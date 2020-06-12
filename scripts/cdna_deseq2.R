library("DESeq2")


dds_snakemake.wildcards <- readRDS(snakemake@input[[1]])

dea_snakemake.wildcards <- DESeq(dds_snakemake.wildcards)

saveRDS(dea_snakemake.wildcards, file=snakemake@output[[1]])

