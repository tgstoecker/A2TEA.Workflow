#Expansion Analysis
# can handle both multi expansion and multi comparison species, now ;D

rule create_hypothesis_fasta:
    input:
        ORTHOFINDER + "complete.check",
    output:
        "tea/{hypothesis}/ref_fasta/hypothesis_{hypothesis}_species.fa",
    params:
        all_hyp_species = get_all_hypothesis_species,
    conda:
        "../envs/expansion.yaml"
    shell:
        "cat {params.all_hyp_species} > {output}"


checkpoint expansion:
    input:
#have to expand because everything has to be finished at this point (at least for now)
        orthology = ORTHOFINDER + "complete.check",
# added the hypothesis fasta to be sure that it is finished here...
        hypo_fasta = "tea/{hypothesis}/ref_fasta/hypothesis_{hypothesis}_species.fa"
    params:
        num = get_hypo_num,
        name = get_hypo_name,
        expansion = get_exp_species,
        comparison = get_com_species,
        add_blast_hits = config["add_blast_hits"],
        expansion_factor = get_expansion_factor,
        expansion_difference = get_expansion_difference,
	ploidy_normalization = get_ploidy_normalization,
        Nmin_expanded_in = get_Nmin_expanded_in,
        Nmin_compared_to = get_Nmin_compared_to,
        expanded_in_all_found = get_expanded_in_all_found,
        compared_to_all_found = get_compared_to_all_found,
    output:
        directory("tea/{hypothesis}/exp_OGs_proteinnames/"),
        directory("checks/tea/{hypothesis}/"),
        directory("tea/{hypothesis}/expansion_tibble/"),
        "tea/{hypothesis}/extended_BLAST_hits/extended_BLAST_hits.RDS",
    threads: 1
    conda:
        "../envs/expansion.yaml"
    script:
        "../scripts/expansion.R"


rule fasta_extraction:
    input:
        protein_lists = "tea/{hypothesis}/exp_OGs_proteinnames/{OG}.txt",
        hypothesis_fasta = "tea/{hypothesis}/ref_fasta/hypothesis_{hypothesis}_species.fa",
    output:
        "tea/{hypothesis}/fa_records/{OG}.fa"
    threads: 1
    conda:
        "../envs/expansion.yaml"
    shell:
        "faSomeRecords {input.hypothesis_fasta} {input.protein_lists} {output}"


rule muscle_MSA:
    input:
        "tea/{hypothesis}/fa_records/{OG}.fa"
    output:
        "tea/{hypothesis}/muscle/{OG}.afa"
    threads: 1
    conda:
        "../envs/expansion.yaml"
    shell:
        "muscle -in {input} -out {output}"


rule trimAl:
    input:
        "tea/{hypothesis}/muscle/{OG}.afa"
    output:
        "tea/{hypothesis}/trimAl/{OG}.afa"
    threads: 1
    conda:
        "../envs/expansion.yaml"
    shell:
        "trimal -automated1 -in {input} -out {output}"


rule FastTree:
    input:
        "tea/{hypothesis}/trimAl/{OG}.afa"
    output:
        "tea/{hypothesis}/trees/{OG}.tree"
    threads: 1
    conda:
        "../envs/expansion.yaml"
    shell:
        "FastTree {input} > {output}"


rule expansion_checkpoint_finish:
    input:
        solve_expansion
    output:
        "checks/expansion/{hypothesis}_finished.txt",
    shell:
        "touch {output}"

