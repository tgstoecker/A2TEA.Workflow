def is_single_end(sample, unit):
    """Determine whether unit is single-end."""
    fq2_present = pd.isnull(samples.loc[(sample, unit), "fq2"])
    if isinstance(fq2_present, pd.core.series.Series):
        # if this is the case, get_fastqs cannot work properly
        raise ValueError(
            f"Multiple fq2 entries found for sample-unit combination {sample}-{unit}.\n"
            "This is most likely due to a faulty units.tsv file, e.g. "
            "a unit name is used twice for the same sample.\n"
            "Try checking your units.tsv for duplicates."
        )
    return fq2_present


def get_fastqs(wildcards):
    """Get raw FASTQ files from unit sheet."""
    if is_single_end(wildcards.sample, wildcards.unit):
        s = samples.loc[ (wildcards.sample, wildcards.unit), ["fq1"] ].dropna()
        return [ f"rawreads/{s.fq1}" ]
    else:
        u = samples.loc[ (wildcards.sample, wildcards.unit), ["fq1", "fq2"] ].dropna()
        return [ f"rawreads/{u.fq1}", f"rawreads/{u.fq2}" ]

def get_trimmed_fastqs(wildcards):
    """Get raw FASTQ files from unit sheet."""
    if is_single_end(wildcards.sample, wildcards.unit):
        s = samples.loc[ (wildcards.sample, wildcards.unit), ["fq1"] ].dropna()
        return [ f"trimmed/{s.fq1}" ]
    else:
        u = samples.loc[ (wildcards.sample, wildcards.unit), ["fq1", "fq2"] ].dropna()
        return [ f"trimmed/{u.fq1}", f"trimmed/{u.fq2}" ]



##SE_samples
def get_SE_fastqs(wildcards):
    """Get raw FASTQ files from unit sheet."""
    if is_single_end(wildcards.sample, wildcards.unit):
        s = SE_samples.loc[ (wildcards.sample, wildcards.unit), ["fq1"] ].dropna()
        return [ f"rawreads/{s.fq1}" ]
    else:
        return

##PE_samples
def get_PE_fastqs(wildcards):
    """Get raw FASTQ files from unit sheet."""
    if not is_single_end(wildcards.sample, wildcards.unit):
        s = PE_samples.loc[ (wildcards.sample, wildcards.unit), ["fq1", "fq2"] ].dropna()
        return [ f"rawreads/{s.fq1}", f"rawreads/{s.fq2}"  ]
    else:
        return


def is_genomic_fasta_present(w):
    """Determine whether species is to be analysed with genomic or transcriptomic read workflow"""
    """If genomic fasta is given in species table, genomic workflow is chosen"""
    genomic_fasta_present = pd.isnull(species_table.loc[(w), "cDNA_fasta"])
    return genomic_fasta_present


def all_species_with_genomic_fasta(sp):
    gen_fasta_species = []
    for i in sp:
        if is_genomic_fasta_present(i) == True:
            gen_fasta_species.append(i)
    return gen_fasta_species


def all_species_with_cDNA_fasta(sp):
    cDNA_fasta_species = []
    for i in sp:
        if is_genomic_fasta_present(i) == False:
            cDNA_fasta_species.append(i)
    return cDNA_fasta_species


GEN_FASTA_SPECIES = all_species_with_genomic_fasta(SPECIES)
CDNA_FASTA_SPECIES = all_species_with_cDNA_fasta(SPECIES)


#create subsets of samples depending on which input is given
CDNA_FASTA_SAMPLES = samples[samples['species'].isin(CDNA_FASTA_SPECIES)]
GEN_FASTA_SAMPLES = samples[samples['species'].isin(GEN_FASTA_SPECIES)]

#write GEN_FASTA_SAMPLES to file for easy use in R later on
GEN_FASTA_SAMPLES.to_csv(r'gen_fasta_samples.csv', sep='\t', index = False)
