#not strictly necessary, since "with open" will overwrite the species.tsv; but I want this to be clean
rule delete_example_species_tsv:
#    input:
#       ahrd_dir =  rules.download_ahrd_wrapper.output,
    output:
        temp(touch("checks/ahrd/ahrd_clean.check"))
    shell:
        "rm AHRD_Snakemake/resources/species.tsv"


rule create_ahrd_species_tsv:
    input:
        old_rem = "checks/ahrd/ahrd_clean.check",
        #expand to subset of species for which AHRD should be run on
        ahrd_fastas = expand("resources/longest_isoforms/{species}.fa", species=AHRD_SPECIES),
    output:
        "AHRD_Snakemake/resources/species.tsv"
    run:
        create_ahrd_species_tsv(input.ahrd_fastas)


#rule run_AHRD_wrapper:

#rule ahrd_done:
#    input:
#        expand("../results/{species}.ahrd_output.tsv", species=SPECIES)


rule handle_user_func_annotations:
    input:
        "config/species.tsv"
    output:
        "resources/functional_annotation/{species}.func_annotation.tsv"
    params:
        species = lambda wildcards: wildcards.species,
        func_annotation= lambda wildcards: species_table[species_table.index == wildcards.species]['function'].item(),
    log:
        "logs/functional_annotation/non_ahrd/{species}.handle.log"
    run:
        handle_user_func_annotation_table(params.func_annotation, params.species)


rule check_user_func_annotations:
    input:
        "resources/functional_annotation/{species}.func_annotation.tsv"
    output:
        touch("checks/functional_annotation/non_ahrd/{species}.check")
    params:
        species = lambda wildcards: wildcards.species,
    log:
        "logs/functional_annotation/non_ahrd/{species}.check.log"
    run:
        check_user_func_annotation_table_validity(input[0], params.species)


#rule combine into R list object
# user supplied as well ahrd based...
rule test_user_func:
    input:
        expand("checks/functional_annotation/non_ahrd/{species}.check", species=NON_AHRD_SPECIES)
    output:
        touch("test_user_func.txt")


