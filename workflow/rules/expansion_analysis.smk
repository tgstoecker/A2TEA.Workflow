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
#        add_blast_hits = config["add_blast_hits"],
        add_OGs = config["add_OGs"],
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
        #for additional OGs instead of BLAST hits
        directory("tea/{hypothesis}/add_OGs_sets/id_lists/"),
        "tea/{hypothesis}/add_OGs_object/add_OG_analysis_object.RDS",
    threads: 1
    conda:
        "../envs/expansion.yaml"
    script:
        "../scripts/expansion.R"



#rule fasta_extraction:
#    input:
#        protein_lists = "tea/{hypothesis}/exp_OGs_proteinnames/{OG}.txt",
#        hypothesis_fasta = "tea/{hypothesis}/ref_fasta/hypothesis_{hypothesis}_species.fa",
#    output:
#        "tea/{hypothesis}/fa_records/{OG}.fa"
#    threads: 1
#    conda:
#        "../envs/expansion.yaml"
#    shell:
#        "faSomeRecords {input.hypothesis_fasta} {input.protein_lists} {output}"


#rule muscle_MSA:
#    input:
#        "tea/{hypothesis}/fa_records/{OG}.fa"
#    output:
#        "tea/{hypothesis}/muscle/{OG}.afa"
#    threads: 1
#    conda:
#        "../envs/expansion.yaml"
#    shell:
#        "muscle -in {input} -out {output}"


#rule trimAl:
#    input:
#        "tea/{hypothesis}/muscle/{OG}.afa"
#    output:
#        "tea/{hypothesis}/trimAl/{OG}.afa"
#    threads: 1
#    conda:
#        "../envs/expansion.yaml"
#    shell:
#        "trimal -automated1 -in {input} -out {output}"


#rule FastTree:
#    input:
#        "tea/{hypothesis}/trimAl/{OG}.afa"
#    output:
#        "tea/{hypothesis}/trees/{OG}.tree"
#    threads: 1
#    conda:
#        "../envs/expansion.yaml"
#    shell:
#        "FastTree {input} > {output}"


#rule expansion_checkpoint_finish:
#    input:
#        solve_expansion
#    output:
#        "checks/expansion/{hypothesis}_finished.txt",
#    shell:
#        "touch {output}"



#########

#since wildcard constrainsts sometimes clash with glob wildcards - see e.g. https://github.com/snakemake/snakemake/issues/482,
#my hacky workaround is to add ".add" to all filenames related to the additionalOGs analysis

## nested snakemake checkpoints are annoying at the moment 
#quick fix - inside the expansion.R script addtional empty set_num files which we can ignore but nonetheless exist are created
#this way, I don't have to worry about evaluating whether or not the user chosen number of addtional OGs were able to be added (for more details see expansion.R script)
#I stick to this trick inside the following rules checkpoint internal rules - if empty file then the biotools aren't used but rather an empty file is created
#in the final_tea script I simply ignore such cases
rule fasta_extraction_add_OG_analysis:
    input:
        protein_lists = "tea/{hypothesis}/add_OGs_sets/id_lists/{OG}/add_OGs_set_num-{set_num}.txt",
        hypothesis_fasta = "tea/{hypothesis}/ref_fasta/hypothesis_{hypothesis}_species.fa",
    output:
        "tea/{hypothesis}/add_OGs_sets/fa_records/{OG}/add_OGs_set_num-{set_num}.fa.add"
    threads: 1
    conda:
        "../envs/expansion.yaml"
    shell:
        """
        if ! [ -s {input.protein_lists} ]; then
          touch {output}
        else
          faSomeRecords {input.hypothesis_fasta} {input.protein_lists} {output}
        fi
        """


rule muscle_MSA_add_OG_analysis:
    input:
#        protein_lists = "tea/{hypothesis}/add_OGs_sets/{OG}/add_OGs_set_num-{set_num}.txt",
        fasta_files = "tea/{hypothesis}/add_OGs_sets/fa_records/{OG}/add_OGs_set_num-{set_num}.fa.add"
    output:
        "tea/{hypothesis}/add_OGs_sets/muscle/{OG}/add_OGs_set_num-{set_num}.afa.add"
    threads: 1
    conda:
        "../envs/expansion.yaml"
    shell:
        """
        if ! [ -s {input.fasta_files} ]; then
          touch {output}
        else
          muscle -in {input.fasta_files} -out {output}
        fi
        """


rule trimAl_add_OG_analysis:
    input:
        "tea/{hypothesis}/add_OGs_sets/muscle/{OG}/add_OGs_set_num-{set_num}.afa.add"
    output:
        "tea/{hypothesis}/add_OGs_sets/trimAl/{OG}/add_OGs_set_num-{set_num}.afa.add"
    threads: 1
    conda:
        "../envs/expansion.yaml"
    shell:
        """
        if ! [ -s {input} ]; then
          touch {output}
        else
          trimal -automated1 -in {input} -out {output}
        fi
        """

rule FastTree_add_OG_analysis:
    input:
        "tea/{hypothesis}/add_OGs_sets/trimAl/{OG}/add_OGs_set_num-{set_num}.afa.add"
    output:
        "tea/{hypothesis}/add_OGs_sets/trees/{OG}/add_OGs_set_num-{set_num}.tree.add"
    threads: 1
    conda:
        "../envs/expansion.yaml"
    shell:
        """
        if ! [ -s {input} ]; then
          touch {output}
        else
          FastTree {input} > {output}
        fi
        """


rule expansion_checkpoint_finish_add_OG_analysis:
    input:
        solve_expansion_add_OG_analysis,
    output:
        "checks/expansion/{hypothesis}_finished_add_OG_analysis.txt",
    shell:
        "touch {output}"


#rule merge_ultimate:
#    input:
#        expand("checks/expansion/{hypothesis}_finished.txt", hypothesis=HYPOTHESES),
#        expand("checks/expansion/{hypothesis}_finished_add_OG_analysis.txt", hypothesis=HYPOTHESES),
#    output:
#        "works.txt"
#    shell:
#        "touch {output}"
