rule final_tea_output:
    input:
        expand("checks/expansion/{hypothesis}_finished.txt", hypothesis=HYPOTHESES),
        expand(rules.cafe5_complete_set.output, hypothesis=HYPOTHESES),
#        expand("{hypothesis}_cafe_check", hypothesis=HYPOTHESES),
    output:
        "tea/A2TEA_finished.RData"
    script:
        "../scripts/final_tea_computation.R"
