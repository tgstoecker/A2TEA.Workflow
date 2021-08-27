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


#rule ahrd_done:
#    input:
#        expand("../results/{species}.ahrd_output.tsv", species=SPECIES)


#rule check user supplied tables

#rule combine into R list object
# user supplied as well ahrd based...
