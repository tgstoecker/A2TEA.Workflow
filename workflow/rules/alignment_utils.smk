# if there is at least 1 species for which a genomic fasta was supplied -> classic alignment with STAR
if len(GEN_FASTA_SPECIES) != 0:

    def get_GEN_FASTA_trimmed_fastqs(wildcards):
        """Get raw FASTQ files from unit sheet."""
        if is_single_end(wildcards.sample, wildcards.unit):
            s = GEN_FASTA_SAMPLES.loc[ (wildcards.sample, wildcards.unit), ["fq1"] ].dropna()
            if not is_gzipped(s.fq1):
                return [ f"trimmed/{s.fq1}.gz" ]
            else:
                return [ f"trimmed/{s.fq1}" ]
        else:
            u = GEN_FASTA_SAMPLES.loc[ (wildcards.sample, wildcards.unit), ["fq1", "fq2"] ].dropna()
            if not is_gzipped(u.fq1):
                return [ f"trimmed/{u.fq1}.gz", f"trimmed/{u.fq2}.gz" ]
            else:
                return [ f"trimmed/{u.fq1}", f"trimmed/{u.fq2}" ]
