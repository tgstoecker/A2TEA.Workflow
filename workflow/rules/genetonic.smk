rule create_GeneTonic_list_files:
    input:
        dea="R/deseq2/dea_final/dea_{species}",
        SFA="results/functional_annotation/species_functional_annotation.rds"
    output:
        "GeneTonic/{species}/GeneTonic_{species}_{ontology}.rds"
    params:
        species="{species}",
        ontology="{ontology}",
        DEG_FDR = config["DEG_FDR"]
    log:
        "logs/GeneTonic/{species}/gtl_{species}_{ontology}.log"
    conda:
        "../envs/genetonic.yaml"
    script:
        "../scripts/genetonic_lists.R"
