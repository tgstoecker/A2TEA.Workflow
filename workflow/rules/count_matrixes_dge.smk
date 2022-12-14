# if there is at least 1 species for which a genomic fasta was supplied -> classic alignment with STAR
if len(GEN_FASTA_SPECIES) != 0:

    #check whether gene or transcript level quantification was chosen
    if config["transcript_level_quantification"] == "NO":

        rule featureCounts_gene_level:
            input:
                bams="star/{species}/{sample}_{unit}_Aligned.sortedByCoord.out.bam",
                gtf="resources/{species}.gtf"
            output:
                gene_level="featureCounts/{species}/{sample}_{unit}_counts.txt",
                gene_level_summary="featureCounts/{species}/{sample}_{unit}_counts.txt.summary",
            params:
                #calling gtf this way is quite janky - change of the index name will lead to an error, because here specfially "species\n" is sliced away
               #gtf = lambda wildcards:species_table.annotation[species_table.index == wildcards.species].to_string(index=False)[8:],
                paired= get_paired_info,
            log:
                gene_level="logs/featureCounts/{species}/{sample}_{unit}_featurecount_gene.log",
                transcript_level="logs/featureCounts/{species}/{sample}_{unit}_featurecount_transcript.log",
            threads:
                config["threads_featureCounts"]
            conda:
                "../envs/subread.yaml"
            shell:
                "featureCounts -T {threads} {params.paired} -O -M -t exon -g gene_id -a {input.gtf} -o {output.gene_level} {input.bams} 2> {log.gene_level} "

    if config["transcript_level_quantification"] == "YES":
        rule featureCounts_transcript_level:
            input:
                bams="star/{species}/{sample}_{unit}_Aligned.sortedByCoord.out.bam",
                gtf="resources/{species}.gtf"
            output:
                transcript_level="featureCounts/{species}/{sample}_{unit}_counts.txt",
                transcript_level_summary="featureCounts/{species}/{sample}_{unit}_counts.txt.summary",
            params:
                #calling gtf this way is quite janky - change of the index name will lead to an error, because here specfially "species\n" is sliced away
                #gtf = lambda wildcards:species_table.annotation[species_table.index == wildcards.species].to_string(index=False)[8:],
                paired= get_paired_info,
            log:
                transcript_level="logs/featureCounts/{species}/{sample}_{unit}_featurecount_transcript.log",
            threads:
                config["threads_featureCounts"]
            conda:
                "../envs/subread.yaml"
            shell:
                "featureCounts -T {threads} {params.paired} -O -M -t exon -g transcript_id -a {input.gtf} -o {output.transcript_level} {input.bams} 2> {log.transcript_level}"



    rule merge_counts:
        input:
            expand("featureCounts/{sample.species}/{sample.sample}_{sample.unit}_counts.txt", sample=GEN_FASTA_SAMPLES.itertuples()),
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
        log:
            "logs/deseq2/{species}.log"
        conda:
            "../envs/deseq2_tximport.yaml"
        script:
            "../scripts/gen_fasta_deseq2.R"



# if there is at least 1 species for which a cDNA fasta was supplied -> pseudoalignment with kallisto
if len(CDNA_FASTA_SPECIES) != 0:


    rule tximport_and_setup:
        input:
            quants=expand("kallisto_quant/{sample.species}/{sample.sample}_{sample.unit}", sample=CDNA_FASTA_SAMPLES.itertuples()),
            annotation="resources/{species}.gtf",
        output:
            "R/tximport/{species}/DESeqDataSet_{species}"
        params:
#            annotation = lambda wildcards: species_table.annotation[species_table.index == wildcards.species],
            species = lambda wildcards: samples.species[samples.species == wildcards.species],
            quantification_level = config["transcript_level_quantification"],
        log:
            "logs/tximport/{species}.log"
        conda:
            "../envs/deseq2_tximport.yaml"
        script:
            "../scripts/tximport.R"


    rule cdna_DESeq2:
        input:
            "R/tximport/{species}/DESeqDataSet_{species}"
        output:
            "R/deseq2/{species}/dea_cdna/dea_{species}"
        log:
            "logs/deseq2/{species}.log"
        conda:
            "../envs/deseq2_tximport.yaml"
        script:
            "../scripts/cdna_deseq2.R"
