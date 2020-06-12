library(DESeq2)

gen_fasta_samples <- read.delim("gen_fasta_samples.csv")

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
dea = DESeq(dds)
saveRDS(dea, file=snakemake@output[[1]])
