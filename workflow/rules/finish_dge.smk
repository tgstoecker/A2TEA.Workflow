# depending on if either or both alignment & pseudoalignment was performed,
# copy final DESeq2 results inside RDS objects to final location in R/ directory 

if len(CDNA_FASTA_SPECIES) != 0:
    if len(GEN_FASTA_SPECIES) != 0:
        rule deseq2_complete:
            input:
                cdna = expand("R/deseq2/{species}/dea_cdna/dea_{species}", species=CDNA_FASTA_SPECIES),
                gen = expand("R/deseq2/{species}/dea_gen/dea_{species}", species=GEN_FASTA_SPECIES),
            output:
                expand("R/deseq2/dea_final/dea_{species}", species=RNA_SPECIES),
            shell:
                "cp R/deseq2/*/dea_cdna/* R/deseq2/dea_final/ && "
                "cp R/deseq2/*/dea_gen/* R/deseq2/dea_final/"


if len(CDNA_FASTA_SPECIES) != 0:
    if len(GEN_FASTA_SPECIES) == 0:
        rule deseq2_complete:
            input:
                cdna = expand("R/deseq2/{species}/dea_cdna/dea_{species}", species=CDNA_FASTA_SPECIES),
            output:
                expand("R/deseq2/dea_final/dea_{species}", species=RNA_SPECIES),
            shell:
                "cp R/deseq2/*/dea_cdna/* R/deseq2/dea_final/"


if len(CDNA_FASTA_SPECIES) == 0:
    if len(GEN_FASTA_SPECIES) != 0:
        rule deseq2_complete:
            input:
                gen = expand("R/deseq2/{species}/dea_gen/dea_{species}", species=GEN_FASTA_SPECIES),
            output:
                expand("R/deseq2/dea_final/dea_{species}", species=RNA_SPECIES),
            shell:
                "cp R/deseq2/*/dea_gen/* R/deseq2/dea_final/"
