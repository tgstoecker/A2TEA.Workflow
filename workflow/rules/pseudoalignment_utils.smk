# if there is at least 1 species for which a cDNA fasta was supplied -> pseudoalignment with kallisto
if len(CDNA_FASTA_SPECIES) != 0:

    def get_paired_info(wildcards):
        """Get single/paired sample information from sample sheet."""
        opt = ""
        if not is_single_end(wildcards.sample, wildcards.unit):
            return [ f"" ]
        else:
            opt += "--single "
            opt += ("--fragment-length {sample.fragment_length_mean} "
                    "--sd {sample.fragment_length_sd}").format(
                           sample=CDNA_FASTA_SAMPLES.loc[(wildcards.sample, wildcards.unit)])
            return opt


    def get_CDNA_FASTA_trimmed_fastqs(wildcards):
        """Get raw FASTQ files from unit sheet."""
        if is_single_end(wildcards.sample, wildcards.unit):
            s = CDNA_FASTA_SAMPLES.loc[ (wildcards.sample, wildcards.unit), ["fq1"] ].dropna()
            return [ f"trimmed/{s.fq1}" ]
        else:
                u = CDNA_FASTA_SAMPLES.loc[ (wildcards.sample, wildcards.unit), ["fq1", "fq2"] ].dropna()
                return [ f"trimmed/{u.fq1}", f"trimmed/{u.fq2}" ]
