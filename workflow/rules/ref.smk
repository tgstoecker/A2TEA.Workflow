rule annotation_handling:
    input:
        "config/species.tsv",
    output:
        temp("resources/{species}.annotation")
    params:
        species = lambda wildcards: wildcards.species,
        annotation = lambda wildcards: species_table[species_table.index == wildcards.species]['annotation'].item(),
    log: 
        "logs/annotation_handling/handling_{species}_annotation.log"
    run:
        handle_annotation(params.annotation, params.species)


rule annotation_standardization:
    input:
        "resources/{species}.annotation"
    output:
        "resources/{species}.gtf"
    log:
        "logs/annotation_standardization/{species}_annotation_standardization.log"
    conda:
        "../envs/gffread.yaml"
    shell:
        "gffread {input} -T -o {output}"


rule pep_fasta_handling:
    input:
        "config/species.tsv",
    output:
        "resources/{species}.pep.fa"
    params:
        species = lambda wildcards: wildcards.species,
        pep_fasta = lambda wildcards: species_table[species_table.index == wildcards.species]['pep_fasta'].item(),
    log:
        "logs/handling_pep_fasta/handling_{species}_pep_fasta.log"
#this rule needs a wildcard_constraint, to prevent it from competing with the "filter_isoforms" rule 
#this is done easily by preventing the wildcard "species" here, from being the string 'longest_isoform' 
#some downstream changes made this superfluous, keeping it here for if the probem should arise again later
#    wildcard_constraints:
#         species='[^(longest_isoform)][0-9a-zA-Z]*'
    run:
        handle_pep_fasta(params.pep_fasta, params.species)


rule cdna_fasta_handling:
    input:
        "config/species.tsv",
    output:
        "resources/{species}.cdna.fa"
    params:
        species = lambda wildcards: wildcards.species,
        cdna_fasta = lambda wildcards: species_table[species_table.index == wildcards.species]['cDNA_fasta'].item(),
    log:
        "logs/handling_pep_fasta/handling_{species}_cdna_fasta.log"
    run:
        handle_cdna_fasta(params.cdna_fasta, params.species)


rule gen_fasta_handling:
    input:
        "config/species.tsv",
    output:
        "resources/{species}.gen.fa"
    params:
        species = lambda wildcards: wildcards.species,
        gen_fasta = lambda wildcards: species_table[species_table.index == wildcards.species]['gen_fasta'].item(),
    log:
        "logs/handling_pep_fasta/handling_{species}_gen_fasta.log"
    run:
        handle_gen_fasta(params.gen_fasta, params.species)
