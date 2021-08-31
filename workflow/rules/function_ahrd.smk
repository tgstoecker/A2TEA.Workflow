##PART 1 - dealing with running AHRD for species the user did not supply their own functional annotation file

#not strictly necessary, since "with open" will overwrite the species.tsv; but I want this to be clean
rule delete_example_species_tsv:
    output:
        temp(touch("checks/ahrd/ahrd_clean.check"))
    shell:
        "rm workflow/rules/AHRD_Snakemake/resources/species.tsv"


checkpoint create_ahrd_species_tsv:
    input:
        old_rem = "checks/ahrd/ahrd_clean.check",
        #expand to subset of species for which AHRD should be run on
        ahrd_fastas = expand("resources/longest_isoforms/{species}.fa", species=AHRD_SPECIES),
    output:
        touch("checks/ahrd/new_ahrd_species_tsv_incl_{species}")
    run:
        create_ahrd_species_tsv(input.ahrd_fastas)


#rule combine into R list object
# user supplied as well ahrd based...

#depending rule of the checkpoint needs to get function as input
#AHRD_snakemake rules are performed between this rule and the checkpoint
rule ahrd_checkpoint_end:
    input:
        species_tsv_checkpoint_end,
    output:
        touch("checks/ahrd/{species}.species_tsv_checkpoint_end.check")
#    shell:
#        touch({output})
#   script:
#        "Rscript that creates RDATA object with list object inside for each species table"


rule ahrd_done:
    input:
        expand("checks/ahrd/{species}.species_tsv_checkpoint_end.check", species=AHRD_SPECIES)
    output:
        touch("checks/ahrd/ahrd_done.txt")


############################################################################################

##PART 2 - dealing with user supplied functional annotation input

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
