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
            "../scripts/gen_fasta_deseq2.R"


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


    rule tximport_and_setup:
        input:
            expand("kallisto_quant/{sample.species}/{sample.sample}_{sample.unit}", sample=CDNA_FASTA_SAMPLES.itertuples()),
        output:
            "R/tximport/{species}/DESeqDataSet_{species}"
        params:
            annotation = lambda wildcards: species_table.annotation[species_table.index == wildcards.species],
            species = lambda wildcards: samples.species[samples.species == wildcards.species],
        script:
            "../scripts/tximport.R"


    rule cdna_DESeq2:
        input:
            "R/tximport/{species}/DESeqDataSet_{species}"
        output:
            "R/deseq2/{species}/dea_cdna/dea_{species}"
        script:
            "../scripts/cdna_deseq2.R"
