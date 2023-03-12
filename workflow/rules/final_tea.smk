rule final_tea_output:
    input:
        expand("R/deseq2/dea_final/dea_{species}", species=RNA_SPECIES),
        expand("checks/expansion/{hypothesis}_finished_add_OG_analysis.txt", hypothesis=HYPOTHESES),
        expand(rules.cafe5_complete_set.output, hypothesis=HYPOTHESES),
        SFA = "results/functional_annotation/species_functional_annotation.rds",
    output:
        "tea/A2TEA_finished.RData"
    params:
        DEG_FDR = config["DEG_FDR"]
    log: 
        "logs/final_tea/final_tea.log"
    conda:
        "../envs/final_tea.yaml"
    script:
        "../scripts/final_tea_computation.R"
