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
        touch("checks/ahrd/new_ahrd_species_tsv_incl_{species}.check")
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

#can be removed once I have a rule that combines AHRD/nonAHRD handling
#rule ahrd_done:
#    input:
#        expand("checks/ahrd/{species}.species_tsv_checkpoint_end.check", species=AHRD_SPECIES)
#    output:
#        touch("checks/ahrd/ahrd_done.check")


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


###Part 3 - merging all functional annotation tables inside Rscript - list data structure with a per species func. annotation table
rule create_func_annotation_RDS:
    input:
        #ahrd_done
        ahrd_checks = expand("checks/ahrd/{species}.species_tsv_checkpoint_end.check", species=AHRD_SPECIES),
#        ahrd_done = "checks/ahrd/ahrd_done.check",
        #ahrd files
#        ahrd_files = expand("results/{species}.ahrd_output.tsv", species=AHRD_SPECIES),
        #checked user supplied func. annotation tables
        user_file_checks = expand("checks/functional_annotation/non_ahrd/{species}.check", species=NON_AHRD_SPECIES),
        #get renamed user supplied func. annotation tables
        renamed_user_files = expand("resources/functional_annotation/{species}.func_annotation.tsv", species=NON_AHRD_SPECIES),
    output:
        "results/functional_annotation/species_functional_annotation.rds"
    params:
        #ahrd files - need to be called via params since otherwise our checkpoint logic won't work (since we use easy expand() here)
        ahrd_files = expand("results/{species}.ahrd_output.tsv", species=AHRD_SPECIES),
    log:
        "logs/functional_annotation/func_annotation_combine.log"
    conda:
        "../envs/func_annotation.yaml"
    script:
        "../scripts/func_annotation_combine.R"
