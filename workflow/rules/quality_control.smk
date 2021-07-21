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



# multiqc
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



