import pandas as pd
import os
import yaml
from snakemake.utils import validate, min_version

##### load config and define samples ####

configfile: "config.yaml"
species_table = pd.read_table("species.tsv").set_index("species")
#species = list(species_table.species.unique())
SPECIES = species_table.index.tolist()

ORTHOFINDER = "orthofinder/"
N_SPECIES = len(SPECIES)
CHUNKS_ORTHO = 2


samples = pd.read_csv("samples.tsv", dtype=str, sep="\t").set_index(["sample", "unit"], drop=False)
samples.index = samples.index.set_levels( [i.astype(str) for i in samples.index.levels])  # enforce str in index

#create SE and PE subset (makes trimmomatic much easier)
SE_samples = samples[samples["fq2"].isna()].set_index(["sample", "unit"], drop=False)
SE_samples.index = SE_samples.index.set_levels( [i.astype(str) for i in SE_samples.index.levels])  # enforce str in index

PE_samples = samples[samples["fq2"].notna()]


wildcard_constraints:
    unit="|".join(samples["unit"]),


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


#for the expansion analysis
hypotheses = pd.read_csv("hypotheses.tsv", dtype=str, sep="\t").set_index(["hypothesis"], drop=False)

HYPOTHESES = hypotheses.index.tolist()
N_HYPOTHESES = len(HYPOTHESES)
N_HYPOTHESES

#pandas actually truncates too long names.. if you don't stop it ;D
pd.set_option("display.max_colwidth", 10000)


rule all:
    input:
        "multiqc/multiqc_report.html",
        expand("R/deseq2/dea_final/dea_{species}", species=SPECIES),
        ORTHOFINDER + "complete.check",
###        expand("{hypo.expanded_in}.ok", hypo=hypotheses.itertuples()),
##        expand("checks/tea/{hypothesis}/exp_OGs_proteinnames.check", hypothesis=HYPOTHESES),
##        expand("tea/{hypothesis}/exp_OGs_proteinnames/", hypothesis=HYPOTHESES),
##        dynamic("tea/{hypothesis}/exp_OGs_proteinnames/proteinnames_{OG}.txt"),
#        expand("checks/expansion/{hypothesis}_finished.txt", hypothesis=HYPOTHESES),
        "tea/A2TEA_finished.RData",
##        expand("tea/{hypothesis}/expansion_tibble/expansion_tibble.rds", hypothesis=HYPOTHESES),
        expand("tea/{hypothesis}/ref_fasta/hypothesis_{hypothesis}_species.fa", hypothesis=HYPOTHESES),
##        expand("tea/{hypothesis}/fa_records/{OG}.fa", hypothesis=HYPOTHESES),




#this will create (if needed) the directories for the trimmomatic rules
FQC_RAW_DIR = "fastqc/raw/"
if not os.path.exists(FQC_RAW_DIR):
    os.makedirs(FQC_RAW_DIR)


FQC_TRIM_DIR = "fastqc/trimmed/"
if not os.path.exists(FQC_TRIM_DIR):
    os.makedirs(FQC_TRIM_DIR)


#rawreads fastqc
rule FastQC_raw:
    input:
        "rawreads/{file}",
    output:
        touch("checks/fastqc/raw/{file}.check"),
    threads: 4
    log:
        "logs/fastqc/raw/{file}.log"
    shell:
        "fastqc -t {threads} -o fastqc/raw/ {input}"


#trimmed reads fastqc
rule FastQC_trimmed:
    input:
        trim_check = "checks/trimmed/trim_cleanup.check",
        raw = "checks/fastqc/raw/{file}.check",
    output:
        touch("checks/fastqc/trimmed/{file}.check"),
    params:
        file = "trimmed/{file}"
    threads: 4
    log:
        "logs/fastqc/trimmed/{file}.log"
    shell:
        "fastqc -t {threads} -o fastqc/trimmed/ {params.file}"



####there has to be a better way than to have three different multiqc rules here...

if len(GEN_FASTA_SPECIES) != 0:
    if len(CDNA_FASTA_SPECIES) != 0:
        rule multiqc:
            input:
                fqc_raw_1 = expand("checks/fastqc/raw/{sample.fq1}.check", sample=SE_samples.itertuples()),
                fqc_raw_2 = expand("checks/fastqc/raw/{sample.fq1}.check", sample=PE_samples.itertuples()),
                fqc_raw_3 = expand("checks/fastqc/raw/{sample.fq2}.check", sample=PE_samples.itertuples()),
                fqc_trimmed_1 = expand("checks/fastqc/trimmed/{sample.fq1}.check", sample=SE_samples.itertuples()),
                fqc_trimmed_2 = expand("checks/fastqc/trimmed/{sample.fq1}.check", sample=PE_samples.itertuples()),
                fqc_trimmed_3 = expand("checks/fastqc/trimmed/{sample.fq2}.check", sample=PE_samples.itertuples()),
                trimmomatic_SE = expand("checks/trimmed/{sample.sample}_{sample.unit}_SE.check", sample=SE_samples.itertuples()),
                trimmomatic_PE = expand("checks/trimmed/{sample.sample}_{sample.unit}_PE.check", sample=PE_samples.itertuples()),
                kallisto = expand("kallisto_quant/{sample.species}/{sample.sample}_{sample.unit}", sample=CDNA_FASTA_SAMPLES.itertuples()),
                star = expand("star/{sample.species}/{sample.sample}_{sample.unit}_Log.final.out", sample=GEN_FASTA_SAMPLES.itertuples()),
                featureCounts = expand("featureCounts/{sample.species}/gene_level/{sample.sample}_{sample.unit}_counts.txt.summary", sample=GEN_FASTA_SAMPLES.itertuples()),
            output:
                "multiqc/multiqc_report.html"
            params:
                plot = "-ip",
                trimmomatic = "logs/trimmomatic/",
                kallisto = "logs/kallisto/quant/",
                star = "star/",
                featureCounts = "featureCounts/",
            log:
                "logs/multiqc/multiqc.log"
            shell:
                "multiqc {params.plot} -d fastqc/raw/ fastqc/trimmed/ {params.trimmomatic} {params.kallisto} {params.star} {params.featureCounts} -o multiqc/"


if len(GEN_FASTA_SPECIES) != 0:
    if len(CDNA_FASTA_SPECIES) == 0:
        rule multiqc:
            input:
                fqc_raw_1 = expand("checks/fastqc/raw/{sample.fq1}.check", sample=SE_samples.itertuples()),
                fqc_raw_2 = expand("checks/fastqc/raw/{sample.fq1}.check", sample=PE_samples.itertuples()),
                fqc_raw_3 = expand("checks/fastqc/raw/{sample.fq2}.check", sample=PE_samples.itertuples()),
                fqc_trimmed_1 = expand("checks/fastqc/trimmed/{sample.fq1}.check", sample=SE_samples.itertuples()),
                fqc_trimmed_2 = expand("checks/fastqc/trimmed/{sample.fq1}.check", sample=PE_samples.itertuples()),
                fqc_trimmed_3 = expand("checks/fastqc/trimmed/{sample.fq2}.check", sample=PE_samples.itertuples()),
                trimmomatic_SE = expand("checks/trimmed/{sample.sample}_{sample.unit}_SE.check", sample=SE_samples.itertuples()),
                trimmomatic_PE = expand("checks/trimmed/{sample.sample}_{sample.unit}_PE.check", sample=PE_samples.itertuples()),
                star = expand("star/{sample.species}/{sample.sample}_{sample.unit}_Log.final.out", sample=GEN_FASTA_SAMPLES.itertuples()),
                featureCounts = expand("featureCounts/{sample.species}/gene_level/{sample.sample}_{sample.unit}_counts.txt.summary", sample=GEN_FASTA_SAMPLES.itertuples()),
            output:
                "multiqc/multiqc_report.html"
            params:
                plot = "-ip",
                trimmomatic = "logs/trimmomatic/",
                star = "star/",
                featureCounts = "featureCounts/",
            log:
                "logs/multiqc/multiqc.log"
            shell:
                "multiqc {params.plot} -d fastqc/raw/ fastqc/trimmed/ {params.trimmomatic} {params.star} {params.featureCounts} -o multiqc/"


if len(GEN_FASTA_SPECIES) == 0:
    if len(CDNA_FASTA_SPECIES) != 0:
        rule multiqc:
            input:
                fqc_raw_1 = expand("checks/fastqc/raw/{sample.fq1}.check", sample=SE_samples.itertuples()),
                fqc_raw_2 = expand("checks/fastqc/raw/{sample.fq1}.check", sample=PE_samples.itertuples()),
                fqc_raw_3 = expand("checks/fastqc/raw/{sample.fq2}.check", sample=PE_samples.itertuples()),
                fqc_trimmed_1 = expand("checks/fastqc/trimmed/{sample.fq1}.check", sample=SE_samples.itertuples()),
                fqc_trimmed_2 = expand("checks/fastqc/trimmed/{sample.fq1}.check", sample=PE_samples.itertuples()),
                fqc_trimmed_3 = expand("checks/fastqc/trimmed/{sample.fq2}.check", sample=PE_samples.itertuples()),
                trimmomatic_SE = expand("checks/trimmed/{sample.sample}_{sample.unit}_SE.check", sample=SE_samples.itertuples()),
                trimmomatic_PE = expand("checks/trimmed/{sample.sample}_{sample.unit}_PE.check", sample=PE_samples.itertuples()),
                kallisto = expand("kallisto_quant/{sample.species}/{sample.sample}_{sample.unit}", sample=CDNA_FASTA_SAMPLES.itertuples()),
            output:
                "multiqc/multiqc_report.html"
            params:
                plot = "-ip",
                trimmomatic = "logs/trimmomatic/",
                kallisto = "logs/kallisto/quant/",
            log:
                "logs/multiqc/multiqc.log"
            shell:
                "multiqc {params.plot} -d fastqc/raw/ fastqc/trimmed/ {params.trimmomatic} {params.kallisto} -o multiqc/"


########
#Trimmomatic
########

ruleorder: Trimmomatic_PE > Trimmomatic_SE


def get_trimmed_PE_output(wildcards):
    if not is_single_end(wildcards.sample, wildcards.unit):
        s = samples.loc[ (wildcards.sample, wildcards.unit), ["fq1", "fq2"] ].dropna()
        return [ f"trimmed/{s.fq1}", f"trimmed/unpaired/unpaired_{s.fq1}", f"trimmed/{s.fq2}", f"trimmed/unpaired/unpaired_{s.fq2}"  ]


def get_trimmed_SE_output(wildcards):
    if is_single_end(wildcards.sample, wildcards.unit):
        s = samples.loc[ (wildcards.sample, wildcards.unit), ["fq1"] ].dropna()
        return [ f"trimmed/{s.fq1}"  ]



#this will create the directories needed for the trimmomatic rules
TRIM_DIR = "trimmed/unpaired/"
if not os.path.exists(TRIM_DIR):
    os.makedirs(TRIM_DIR)


rule Trimmomatic_PE:
    input:
        get_PE_fastqs,
    output:
        touch("checks/trimmed/{sample}_{unit}_PE.check"),
    log:
        "logs/trimmomatic/PE/{sample}_{unit}.log"
    params:
        illuminaclip=config["trim_illuminaclip"],
        trimmer=config["trim_options"],
        compression_level="-9",
        out = get_trimmed_PE_output,
    threads:
        config["threads_trimmomatic"]
    shell:
        "trimmomatic PE -threads {threads} {input} {params.out} {params.illuminaclip} {params.trimmer} 2> {log}"


rule Trimmomatic_SE:
    input:
        get_SE_fastqs,
    output:
        touch("checks/trimmed/{sample}_{unit}_SE.check"),
    log:
        "logs/trimmomatic/SE/{sample}_{unit}.log"
    params:
        illuminaclip=config["trim_illuminaclip"],
        trimmer=config["trim_options"],
        compression_level="-9",
        out = get_trimmed_SE_output,
    threads:
        config["threads_trimmomatic"]
    shell:
        "trimmomatic SE -threads {threads} {input} {params.out} {params.illuminaclip} {params.trimmer} 2> {log}"


rule Trimmomatic_cleanup:
    input:
        expand("checks/trimmed/{sample.sample}_{sample.unit}_SE.check", sample=SE_samples.itertuples()),
        expand("checks/trimmed/{sample.sample}_{sample.unit}_PE.check", sample=PE_samples.itertuples()),
    output:
        touch("checks/trimmed/trim_cleanup.check"),
    shell:
        "rm -r trimmed/unpaired/"



if len(GEN_FASTA_SPECIES) != 0:
#######
#STAR
#######
    rule STAR_index:
        input:
            fasta = lambda wildcards: species_table.gen_fasta[species_table.index == wildcards.species],
            annotation = lambda wildcards: species_table.annotation[species_table.index == wildcards.species],
        output:
            directory("STAR_indexes/{species}")
        message:
            "Creating STAR index"
        params:
            extra = "",
            threads= config["threads_star_index"],
            length = config["read_length_star_index"],
            size = config["limitGenomeGenerateRAM"],
            temp_dir = "STAR_tmp_{species}",
            #possible to add params here, by referring to the table like this:
            #SampleSM = lambda wildcards: list(units_table.SampleSM[units_table.Unit == wildcards.unit]),
#        log:
#            "logs/star_index/{species}/Log.out"
        shell:
            'mkdir {output} && '
            'STAR --runThreadN {params.threads} '
            '--runMode genomeGenerate '
            '--genomeDir {output} '
            '--limitGenomeGenerateRAM {params.size} '
            '--genomeFastaFiles {input.fasta} '
            '--sjdbGTFfile {input.annotation} '
            '--sjdbOverhang {params.length} '
            '--outTmpDir {params.temp_dir}'


    def get_GEN_FASTA_trimmed_fastqs(wildcards):
        """Get raw FASTQ files from unit sheet."""
        if is_single_end(wildcards.sample, wildcards.unit):
            s = GEN_FASTA_SAMPLES.loc[ (wildcards.sample, wildcards.unit), ["fq1"] ].dropna()
            return [ f"trimmed/{s.fq1}" ]
        else:
            u = GEN_FASTA_SAMPLES.loc[ (wildcards.sample, wildcards.unit), ["fq1", "fq2"] ].dropna()
            return [ f"trimmed/{u.fq1}", f"trimmed/{u.fq2}" ]


    rule STAR_align:
            input:
                trim_check = "checks/trimmed/trim_cleanup.check",
                dir = expand("STAR_indexes/{species}", species = GEN_FASTA_SPECIES),
            output:
               # see STAR manual for additional output files -
                align = "star/{species}/{sample}_{unit}_Aligned.sortedByCoord.out.bam",
#                log = "star/{species}/{sample}_{unit}_Log.final.out"
            log:
                "star/{species}/{sample}_{unit}_Log.final.out"
            threads: config["threads_star"]
            params:
                correct_genome = lambda wildcards: GEN_FASTA_SAMPLES.loc[(wildcards.sample, wildcards.unit) , 'species'],
                rest = config["STAR"],
                sample = get_GEN_FASTA_trimmed_fastqs,
            shell:
                'STAR --runThreadN {threads} '
                '--genomeDir STAR_indexes/{params.correct_genome} '
                '--readFilesIn {params.sample} '
                '--readFilesCommand zcat '
                '--outFileNamePrefix star/{params.correct_genome}/{wildcards.sample}_{wildcards.unit}_ '
                '{params.rest}'


    rule index_BAMs:
        input:
            "star/{species}/{sample}_{unit}_Aligned.sortedByCoord.out.bam"
        output:
            "star/{species}/{sample}_{unit}_Aligned.sortedByCoord.out.bam.csi"
        params:
#        threads = config["threads_index_sorted_bams_with_dups"],
            index_type = "-c"
        threads:
            config["threads_index_sorted_bams_with_dups"],
        shell:
            "samtools index {params.index_type} -@ {threads} {input}"


    def get_paired_info(wildcards):
        """Get raw FASTQ files from unit sheet."""
        if not is_single_end(wildcards.sample, wildcards.unit):
            return [ f"-p" ]
        else:
            return [ f"" ]


    rule featureCounts:
        input:
            bams="star/{species}/{sample}_{unit}_Aligned.sortedByCoord.out.bam",
        output:
            gene_level="featureCounts/{species}/gene_level/{sample}_{unit}_counts.txt",
#            transcript_level="featureCounts/{species}/transcript_level/{sample}_{unit}_counts.txt",
            gene_level_summary="featureCounts/{species}/gene_level/{sample}_{unit}_counts.txt.summary",
        params:
            #calling gtf this way is quite janky - change of the index name will lead to an error, because here specfially "species\n" is sliced away
            gtf = lambda wildcards:species_table.annotation[species_table.index == wildcards.species].to_string(index=False)[9:],
            paired= get_paired_info,
        log:
            gene_level="logs/featureCounts/{species}/{sample}_{unit}_featurecount_gene.log",
            transcript_level="logs/featureCounts/{species}/{sample}_{unit}_featurecount_transcript.log",
        threads:
            config["threads_featureCounts"]
        shell:
            "featureCounts -T {threads} {params.paired} -O -M -t exon -g gene_id -a {params.gtf} -o {output.gene_level} {input.bams} 2> {log.gene_level} "
#        "featureCounts -T {threads} {params.paired} -O -M -t exon -g transcript_id -a {params.gtf} -o {output.transcript_level} {input.bams} 2> {log.transcript_level} "


    def getCountsForSpecies(wildcards):
        counts = list()
        # only the actual counts (ending with .txt) are considered in the diectory
        # also the output is sorted - when sticking to control and treatment everything is as it should
        # control 1-4, treatment 1-4
        for c in sorted(os.listdir("featureCounts/"+wildcards.species+"/gene_level/")):
            if c.endswith(".txt"):
                counts.append(os.path.join("featureCounts/",wildcards.species,"gene_level/",c))
        return counts


    rule merge_counts:
        input:
            expand("featureCounts/{sample.species}/gene_level/{sample.sample}_{sample.unit}_counts.txt", sample=GEN_FASTA_SAMPLES.itertuples()),
        output:
            "featureCounts/{species}/gene_level/counts_merged.txt",
        params:
            counts=getCountsForSpecies,
        run:
            # Merge count files.
            frames = (pd.read_csv(fp, sep="\t", skiprows=1,
                            index_col=list(range(6)))
                for fp in params.counts)
            merged = pd.concat(frames, axis=1)
            merged.to_csv(output[0], sep="\t", index=True)


    rule gen_DESeq2:
        input:
            "featureCounts/{species}/gene_level/counts_merged.txt",
        output:
            "R/deseq2/{species}/dea_gen/dea_{species}"
        params:
            species = lambda wildcards: samples.species[samples.species == wildcards.species],
        script:
            "scripts/gen_fasta_deseq2.R"



if len(CDNA_FASTA_SPECIES) != 0:
##########
#kallisto
##########
    rule kallisto_index:
        input:
            fasta = lambda wildcards: species_table.cDNA_fasta[species_table.index == wildcards.species],
        output:
            "kallisto_indexes/{species}.idx"
        log:
            "logs/kallisto/indexes/{species}.log"
        threads: 1
        shell:
            "kallisto index -i {output} {input.fasta}"


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



    rule kallisto_quant:
        input:
            trim_check = "checks/trimmed/trim_cleanup.check",
            index = "kallisto_indexes/{species}.idx",
        output:
            directory("kallisto_quant/{species}/{sample}_{unit}")
        log:
            "logs/kallisto/quant/{species}/{sample}_{unit}.quant.log"
        params:
            paired = get_paired_info,
            input = get_CDNA_FASTA_trimmed_fastqs,
            bootstrap = "0", # 100; if we wanted to work on transcript level and make use of the bootstraps
        threads:
            config["threads_kallisto_quant"]
        shell:
            "kallisto quant -i {input.index} -o {output} -b {params.bootstrap} -t {threads} "
            "{params.paired} {params.input} 2> {log}"


########
#diff. exp
########
    rule tximport_and_setup:
        input:
            expand("kallisto_quant/{sample.species}/{sample.sample}_{sample.unit}", sample=CDNA_FASTA_SAMPLES.itertuples()),
        output:
            "R/tximport/{species}/DESeqDataSet_{species}"
        params:
            annotation = lambda wildcards: species_table.annotation[species_table.index == wildcards.species],
            species = lambda wildcards: samples.species[samples.species == wildcards.species],
        script:
            "scripts/tximport.R"


    rule cdna_DESeq2:
        input:
            "R/tximport/{species}/DESeqDataSet_{species}"
        output:
            "R/deseq2/{species}/dea_cdna/dea_{species}"
        script:
            "scripts/cdna_deseq2.R"



if len(CDNA_FASTA_SPECIES) != 0:
    if len(GEN_FASTA_SPECIES) != 0:
        rule deseq2_complete:
            input:
                cdna = expand("R/deseq2/{species}/dea_cdna/dea_{species}", species=CDNA_FASTA_SPECIES),
                gen = expand("R/deseq2/{species}/dea_gen/dea_{species}", species=GEN_FASTA_SPECIES),
            output:
                expand("R/deseq2/dea_final/dea_{species}", species=SPECIES),
            shell:
                "cp R/deseq2/*/dea_cdna/* R/deseq2/dea_final/ && "
                "cp R/deseq2/*/dea_gen/* R/deseq2/dea_final/"


if len(CDNA_FASTA_SPECIES) != 0:
    if len(GEN_FASTA_SPECIES) == 0:
        rule deseq2_complete:
            input:
                cdna = expand("R/deseq2/{species}/dea_cdna/dea_{species}", species=CDNA_FASTA_SPECIES),
            output:
                expand("R/deseq2/dea_final/dea_{species}", species=SPECIES),
            shell:
                "cp R/deseq2/*/dea_cdna/* R/deseq2/dea_final/"


if len(CDNA_FASTA_SPECIES) == 0:
    if len(GEN_FASTA_SPECIES) != 0:
        rule deseq2_complete:
            input:
                gen = expand("R/deseq2/{species}/dea_gen/dea_{species}", species=GEN_FASTA_SPECIES),
            output:
                expand("R/deseq2/dea_final/dea_{species}", species=SPECIES),
            shell:
                "cp R/deseq2/*/dea_gen/* R/deseq2/dea_final/"



#####
#Orthofinder
#####

def get_species_fasta(wildcards):
    return species_table.loc[wildcards.species]["pep_fasta"]


def get_longest_isoforms(wildcards):
    return os.path.join("FS/longest_isoforms/", os.path.split(get_species_fasta(wildcards))[1])


if config["auto_isoform_filtering"] == "YES":
    rule filter_isoforms:
        output:
            "FS/longest_isoforms/{species}.fa"
        params:
            fa = get_species_fasta,
            iso = get_longest_isoforms,
        shell:
            "python scripts/longest_isoforms.py {params.fa} && "
            "mv {params.iso} {output}"


    rule Orthofinder_link_all:
        input:
            "FS/longest_isoforms/{species}.fa",
        output:
            ORTHOFINDER + "{species}.fa"
        params: get_longest_isoforms,
        shell: 
            "ln --symbolic $(readlink --canonicalize {input}) {output}"


else:
    rule Orthofinder_link_all:
        output:
            ORTHOFINDER + "{species}.fa"
        params:
            fa = get_species_fasta,
        shell:
            "ln --symbolic $(readlink --canonicalize {params.fa}) {output}"



rule Orthofinder_prepare:
    "Split fasta files, rename species and sequences and prepare blast databases"
    input:
        fastas = expand(ORTHOFINDER + "{species}.fa", species=SPECIES)
    output:
        check = touch(ORTHOFINDER + "prepare.check"),
        fastas = expand(ORTHOFINDER + "Species{species_number}.fa", species_number=[x for x in range(0, N_SPECIES)]),
        db = expand(ORTHOFINDER + "diamondDBSpecies{database_number}.dmnd", database_number=[x for x in range(0, N_SPECIES)])
    params:
        fasta_dir = ORTHOFINDER,
        temp_dir1 = ORTHOFINDER + "OrthoFinder/Results_*",
        temp_dir2 = ORTHOFINDER + "OrthoFinder/Results_*/WorkingDirectory/"
    log:
        "logs/orthofinder/prepare.log"
    benchmark:
        "benchmarks/orthofinder/prepare.json"
#    conda:
#        "orthofinder.yml"
    shell:
        "orthofinder -t 16 "
         "--fasta {params.fasta_dir} "
         "--search diamond "
         "--only-prepare "
         "2> {log} 1>&2 && "
         "mv {params.temp_dir2}* {params.fasta_dir} && "
         "rm --recursive --force {params.temp_dir1}"


rule index_pep_FASTAs:
    input:
        ORTHOFINDER + "Species{species_number}.fa",
    output:
        ORTHOFINDER + "Species{species_number}.fa.fai",
    shell:
        "samtools faidx {input}"



# added choice of chunk usage or not ;D
with open('config.yaml', 'r') as f:
    config = yaml.load(f)

if config["chunk_usage"] == "ON":
    rule Orthofinder_split:
        "Split the headers of the input pep fastas into multiple files"
        input:
            fai = ancient(ORTHOFINDER + "Species{species_number}.fa.fai")
        output:
            expand(ORTHOFINDER + "{{species_number}}/chunks/ids_{chunk_id}.tsv", chunk_id=['{0:05d}'.format(x) for x in range(0, CHUNKS_ORTHO)])
        params:
            folder = ORTHOFINDER,
            number_of_chunks = CHUNKS_ORTHO,
            species = "{species_number}"
        log:
        #send the log here as an exception; also controls the actual output 
            "orthofinder/{species_number}/split.log"
        benchmark:
            "benchmarks/orthofinder/{species_number}/split.json"
    #    conda:
    #        "orthofinder.yml"
        shell:
            "split "
            "--number l/{params.number_of_chunks} "
            "--numeric-suffixes "
            "--suffix-length 5 "
            "--additional-suffix .tsv "
            "{input.fai} "
            "{params.folder}/{params.species}/chunks/ids_ "
            "2> {log}"


    rule Orthofinder_BLASTP:
        "Run blastp of each chunk"
        input:
            fasta = ORTHOFINDER + "Species{species_number}.fa",
            fai   = ancient(ORTHOFINDER + "Species{species_number}.fa.fai"),
            chunk = ORTHOFINDER + "{species_number}/chunks/ids_{chunk_id}.tsv",
            db    = ORTHOFINDER + "diamondDBSpecies{database_number}.dmnd",
        output:
            tsv = ORTHOFINDER + "{species_number}/{database_number}/blastp_{chunk_id}.tsv"
        log:
            "logs/orthofinder/{species_number}/{database_number}/blastp_{chunk_id}.log"
        benchmark:
            "benchmarks/orthofinder/{species_number}/{database_number}/blastp_{chunk_id}.json"
        threads: 6
    #    conda:
    #        "orthofinder.yml"
        shell:
            "cut -f 1 {input.chunk} "
            "| xargs samtools faidx {input.fasta}"
            "| diamond blastp "
            "--ultra-sensitive "
            "--db {input.db} "
            "--outfmt 6 "
            "--evalue 0.001 "
            "--out {output.tsv} "
            "--threads {threads} "
            "2> {log} 1>&2"


    rule Orthofinder_BLASTP_merge:
        "Merge results from the different blastps"
        input:
            expand(ORTHOFINDER + "{{species_number}}/{{database_number}}/blastp_{chunk_id}.tsv", chunk_id = ['{0:05d}'.format(x) for x in range(0, CHUNKS_ORTHO)])
        output:
            tsv = ORTHOFINDER + "Blast{species_number}_{database_number}.txt"
        log:
            "logs/orthofinder/{species_number}/{database_number}/blastp_merge.log"
        benchmark:
            "benchmarks/orthofinder/{species_number}/{database_number}/blastp_merge.json"
    #    conda:
    #        "orthofinder.yml"
        shell:
            "cat {input} > {output} 2> {log}"


## if chunks are not desired or feasible (anything other than "ON" under chunk_usage in the config.yaml)
else:
    rule Orthofinder_BLASTP:
        "Run blastp of each chunk"
        input:
            fasta = ORTHOFINDER + "Species{species_number}.fa",
            fai   = ancient(ORTHOFINDER + "Species{species_number}.fa.fai"),
            db    = ORTHOFINDER + "diamondDBSpecies{database_number}.dmnd",
        output:
            tsv = ORTHOFINDER + "{species_number}/{database_number}/blastp.tsv"
        log:
            "logs/orthofinder/{species_number}/{database_number}/blastp.log"
        benchmark:
            "benchmarks/orthofinder/{species_number}/{database_number}/blastp.json"
        threads: 6
    #    conda:
    #        "orthofinder.yml"
        shell:
            "diamond blastp --query {input.fasta} "
            "--ultra-sensitive "
            "--db {input.db} "
            "--outfmt 6 "
            "--evalue 0.001 "
            "--out {output.tsv} "
            "--threads {threads} "
            "2> {log} 1>&2"  


    rule Orthofinder_BLASTP_merge:
        "Merge results from the different blastps"
        input:
            expand(ORTHOFINDER + "{{species_number}}/{{database_number}}/blastp.tsv")
        output:
            tsv = ORTHOFINDER + "Blast{species_number}_{database_number}.txt"
        log:
            "logs/orthofinder/{species_number}/{database_number}/blastp_merge.log"
        benchmark:
            "benchmarks/orthofinder/{species_number}/{database_number}/blastp_merge.json"
    #    conda:
    #        "orthofinder.yml"
        shell:
            "cat {input} > {output} 2> {log}"


rule Orthofinder_orthogroups:
    "Join blastp results, normalize bitscores and run mcl"
    input:
        tsv = expand(ORTHOFINDER + "Blast{database_number}_{species_number}.txt",
                  species_number = [x for x in range(0,N_SPECIES)],
                      database_number = [x for x in range(0,N_SPECIES)])
    output:
        touch(ORTHOFINDER + "groups.check"),
    params:
        fasta_dir = ORTHOFINDER,
        temp_dir = ORTHOFINDER + "OrthoFinder/Results_*",
#        inflation = params["orthofinder"]["mcl_inflation"],
    threads: 24 # There is no reason to go beyond this value - really?
    log:
        "logs/orthofinder/groups.log"
    benchmark:
        "benchmarks/orthofinder/groups.json"
#    conda:
#        "orthofinder.yml"
    shell:
        "orthofinder "
        "--algthreads {threads} "
        "--inflation 1.5 "
        "--blast {params.fasta_dir} "
        "--only-groups "
        "2> {log} 1>&2 && "
        "mv {params.temp_dir}/ {params.fasta_dir} "


rule Orthofinder_Trees:
    input:
        rules.Orthofinder_orthogroups.output
    output:
        touch(ORTHOFINDER + "trees.check")
    params:
        orthofinder_dir = ORTHOFINDER + "Results_*",
        tree_program = "fasttree",
        msa_program = "muscle",
#        mafft as alternative
#    conda:
#        "orthofinder.yml"
    log:
        "logs/orthofinder/trees.log"
    benchmark:
        "benchmarks/orthofinder/trees.bmk"
    threads: 24
    shell:
        "orthofinder "
        "--from-groups {params.orthofinder_dir} "
        "--only-trees "
        "--method msa "
        "--msa_program {params.msa_program} "
        "--tree_program fasttree "
        "--algthreads {threads} "
        "--threads {threads} "
        "2> {log} 1>&2"


rule Orthofinder_orthologues:
    input:
        ORTHOFINDER + "trees.check"
    output:
        touch(ORTHOFINDER + "orthologues.check")
#    conda: "orthofinder.yml"
    params:
        orthofinder_dir = ORTHOFINDER + "Results_*_1",
    log:
        "logs/orthofinder/orthologues.log"
    threads: 24
    shell:
        "orthofinder "
        "--from-trees {params.orthofinder_dir} "
        "-a {threads} "
        "--threads {threads} "
        "2> {log} 1>&2"


rule Orthofinder_cleanup:
    input:
        ORTHOFINDER + "orthologues.check"
    output:
        touch(ORTHOFINDER + "cleanup.check")
    params:
        orthofinder_dir = ORTHOFINDER,
        n_species = N_SPECIES - 1,
        results_dir_1 = "Results_*_1",
        results_dir_2 = "Results_*_2",
    log:
        "logs/orthofinder/clean.log"
    shell:
        "pushd {params.orthofinder_dir} 2> {log} 1>&2 && "
        "mkdir search/ && "
        "for i in {{0..{params.n_species}}} ; do "
        "mv ${{i}} search/ ; "
        "done && "
        "mv Blast* search/ &&"
        "mv diamond* search/ && "
        ""
        "mv {params.results_dir_1}/Gene_Trees/ {params.results_dir_2}/ && "
        "mv {params.results_dir_1}/MultipleSequenceAlignments/ {params.results_dir_2}/ && "
        "mv {params.results_dir_1}/Orthogroup_Sequences/ {params.results_dir_2}/ && "
        "mv {params.results_dir_1}/WorkingDirectory/SpeciesTree_unrooted* {params.results_dir_2}/WorkingDirectory/ && "
        "mv {params.results_dir_1}/WorkingDirectory/*_ids/ {params.results_dir_2}/WorkingDirectory/ && "
        "mv {params.results_dir_2} final-results/ && "
        "rm -rf Results_*_1 && "
        "mv Results*/Comparative_Genomics_Statistics/* final-results/Comparative_Genomics_Statistics/ && "
        "mv Results*/Orthogroups/ final-results/ && "
        "rm -rf Results_* "


rule Orthofinder_complete:
    input:
        ORTHOFINDER + "cleanup.check"
    output:
        touch(ORTHOFINDER + "complete.check")


##################################################################
##################################################################
#Expansion Analysis
# can handle both multi expansion and multi comparison species, now ;D
# the strplit command works really well in this context
## this way, multiple species are propagated into R as a simple vector  

def get_hypo_num(wildcards):
    """Get compared_to entries from hypotheses(.tsv) for each hypothesis. """
    num = hypotheses.loc[ (wildcards.hypothesis), 'hypothesis']
    return num

def get_hypo_name(wildcards):
    """Get compared_to entries from hypotheses(.tsv) for each hypothesis. """
    name = hypotheses.loc[ (wildcards.hypothesis), 'name']
    return name


# for both the species we want to check for expansion as well as those being compared to:
# if more than one species is compared to than we have to split the string based on ";" 

def get_exp_species(wildcards):
    """Get compared_to entries from hypotheses(.tsv) for each hypothesis. """
    exps_1 = hypotheses.loc[ (wildcards.hypothesis), 'expanded_in']
    if exps_1.count(";") > 0:
        exps_2 = str.split(exps_1, ";")
        return exps_2
    else:
        return exps_1
    return


#if more than one species is compared to than we have to split the string based on ";"
def get_com_species(wildcards):
    """Get compared_to entries from hypotheses(.tsv) for each hypothesis. """
    ct_1 = hypotheses.loc[ (wildcards.hypothesis), 'compared_to']
    if ct_1.count(";") > 0:
        ct_2 = str.split(ct_1, ";")
        return ct_2
    else:
        return ct_1
    return


#if more than one species is compared to than we have to split the string based on ";"
#the goal in any scenarios to return a list (!) with just the species used in the hypothesis
#we also need to add the path to the longest isoform peptide fasta files and add the .fa suffix
#another trick is that we can simply use the species base name in any case,
#since the longest isoform output is always named after it!
#+the filtering for longest isoform also only retains the base name of the gene/protein

def get_all_hypothesis_species(wildcards):
    """Get compared_to entries from hypotheses(.tsv) for each hypothesis. """
    path_prefix = 'FS/longest_isoforms/'
    suffix = '.fa'
    exp = hypotheses.loc[ (wildcards.hypothesis), 'expanded_in']
    ct = hypotheses.loc[ (wildcards.hypothesis), 'compared_to']
    if ct.count(";") > 0:
        ct = str.split(ct, ";")
        ct.append(exp)
        ct = [path_prefix + x + suffix for x in ct]
        return ct
    else:
        output = []
        output.append(exp)
        output.append(ct)
        output = [path_prefix + x + suffix for x in output]
        return output
    return


rule create_hypothesis_fasta:
    input:
        orthology = ORTHOFINDER + "complete.check",
    output:
        "tea/{hypothesis}/ref_fasta/hypothesis_{hypothesis}_species.fa",
    params:
        all_hyp_species = get_all_hypothesis_species,
    shell:
        "cat {params.all_hyp_species} > {output}"


checkpoint expansion:
    input:
#have to expand because everything has to be finished at this point (at least for now)
        orthology = ORTHOFINDER + "complete.check",
# added the hypothesis fasta to be sure that it is finished here...
        hypo_fasta = "tea/{hypothesis}/ref_fasta/hypothesis_{hypothesis}_species.fa"
    params:
        num = get_hypo_num,
        name = get_hypo_name,
        expansion = get_exp_species,
        comparison = get_com_species,
        add_blast_hits = config["add_blast_hits"]
    output:
        directory("tea/{hypothesis}/exp_OGs_proteinnames/"),
        directory("checks/tea/{hypothesis}/"),
        directory("tea/{hypothesis}/expansion_tibble/"),
        "tea/{hypothesis}/extended_BLAST_hits/extended_BLAST_hits.RDS",
    threads: 1
    script:
        "scripts/expansion.R"


rule fasta_extraction:
    input:
        protein_lists = "tea/{hypothesis}/exp_OGs_proteinnames/{OG}.txt",
        hypothesis_fasta = "tea/{hypothesis}/ref_fasta/hypothesis_{hypothesis}_species.fa",
    output:
        "tea/{hypothesis}/fa_records/{OG}.fa"
    threads: 1
    shell:
        "faSomeRecords {input.hypothesis_fasta} {input.protein_lists} {output}"


rule muscle_MSA:
    input:
        "tea/{hypothesis}/fa_records/{OG}.fa"
    output:
        "tea/{hypothesis}/muscle/{OG}.afa"
    threads: 1
    shell:
        "muscle -in {input} -out {output}"


rule trimAl:
    input:
        "tea/{hypothesis}/muscle/{OG}.afa"
    output:
        "tea/{hypothesis}/trimAl/{OG}.afa"
    threads: 1
    shell:
        "trimal -automated1 -in {input} -out {output}"


rule FastTree:
    input:
        "tea/{hypothesis}/trimAl/{OG}.afa"
    output:
        "tea/{hypothesis}/trees/{OG}.tree"
    threads: 1
    shell:
        "FastTree {input} > {output}"


#CHECKPOINTS ARE SO AWESOME!
def solve_expansion(wildcards):
    checkpoint_output = checkpoints.expansion.get(**wildcards).output[0]
    file_names = expand("tea/{hypothesis}/trees/{OG}.tree", hypothesis=wildcards.hypothesis, OG=glob_wildcards(os.path.join(checkpoint_output, "{OG}.txt")).OG)
    return file_names


rule expansion_checkpoint_finish:
    input:
        solve_expansion
    output:
        "checks/expansion/{hypothesis}_finished.txt",
    shell:
        "touch {output}"

# necessary to be done here with the expression analysis;
#        expression = expand("R/deseq2/dea_final/dea_{species}", species=SPECIES),


rule final_tea_output:
    input:
        expand("checks/expansion/{hypothesis}_finished.txt", hypothesis=HYPOTHESES),
    output:
        "tea/A2TEA_finished.RData"
    script:
        "scripts/final_tea_computation.R"
