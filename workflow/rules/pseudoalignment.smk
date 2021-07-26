# if there is at least 1 species for which a cDNA fasta was supplied -> pseudoalignment with kallisto
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
        conda:
            "../envs/kallisto.yaml"
        shell:
            "kallisto index -i {output} {input.fasta}"


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
        conda:
            "../envs/kallisto.yaml"
        shell:
            "kallisto quant -i {input.index} -o {output} -b {params.bootstrap} -t {threads} "
            "{params.paired} {params.input} 2> {log}"
