# if there is at least 1 species for which a genomic fasta was supplied -> classic alignment with STAR
if len(GEN_FASTA_SPECIES) != 0:
#######
#STAR
#######
    rule STAR_index:
        input:
#            fasta = lambda wildcards: species_table.gen_fasta[species_table.index == wildcards.species],
            fasta = "resources/{species}.gen.fa",
#            annotation = lambda wildcards: species_table.annotation[species_table.index == wildcards.species],
            annotation = "resources/{species}.gtf"
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
        conda:
            "../envs/star.yaml"
        log:
            "logs/star_index/{species}/Log.out"
        shell:
            'mkdir {output} && '
            'STAR --runThreadN {params.threads} '
            '--runMode genomeGenerate '
            '--genomeDir {output} '
            '--limitGenomeGenerateRAM {params.size} '
            '--genomeFastaFiles {input.fasta} '
            '--sjdbGTFfile {input.annotation} '
            '--sjdbOverhang {params.length} '
            '--outTmpDir {params.temp_dir} '
            '{log}'


    rule STAR_align:
        input:
            trim_check = "checks/trimmed/trim_cleanup.check",
            dir = expand("STAR_indexes/{species}", species = GEN_FASTA_SPECIES),
        output:
            # see STAR manual for additional output files -
            align = "star/{species}/{sample}_{unit}_Aligned.sortedByCoord.out.bam",
#            log = "star/{species}/{sample}_{unit}_Log.final.out"
        log:
            "star/{species}/{sample}_{unit}_Log.final.out"
        threads: config["threads_star"]
        params:
            correct_genome = lambda wildcards: GEN_FASTA_SAMPLES.loc[(wildcards.sample, wildcards.unit) , 'species'],
            rest = config["STAR"],
            sample = get_GEN_FASTA_trimmed_fastqs,
        conda:
            "../envs/star.yaml"
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
            config["threads_index_sorted_bams_with_dups"]
        conda:
            "../envs/samtools.yaml"
        shell:
            "samtools index {params.index_type} -@ {threads} {input}"
