#Expansion Analysis
# can handle both multi expansion and multi comparison species, now ;D


#initially we reduced the trees, etc. to only those species specifically part of hypothesis
# this was superseded by the idea that a A2TEA run rep. an experiment and in the final analyses all species are included although the cal.
# of expanded OGs is based only on those species that are part of the respective hypothesis
# we keep the following rule commented out as we might return to this in the future
#rule create_hypothesis_fasta:
#    input:
#        ORTHOFINDER + "complete.check",
#    output:
#        "tea/{hypothesis}/ref_fasta/hypothesis_{hypothesis}_species.fa",
#    params:
#        all_hyp_species = get_all_hypothesis_species,
#    conda:
#        "../envs/expansion.yaml"
#    shell:
#        "cat {params.all_hyp_species} > {output}"


rule create_superset_fasta:
    input:
        ORTHOFINDER + "complete.check",
    output:
        "tea/superset_fasta/superset_species.fa",
    params:
        superset_species_fa = create_all_superset_species_fa(),
    conda:
        "../envs/expansion.yaml"
    shell:
        "cat {params.superset_species_fa} > {output}"

rule expansion_computation:
    input:
#have to expand because everything has to be finished at this point (at least for now)
        orthology = ORTHOFINDER + "complete.check",
## added the hypothesis fasta to be sure that it is finished here...
##        hypo_fasta = "tea/{hypothesis}/ref_fasta/hypothesis_{hypothesis}_species.fa"
#        superset_fa = "tea/superset_fasta/superset_species.fa"
    params:
        num = get_hypo_num,
        name = get_hypo_name,
        expansion = get_exp_species,
        comparison = get_com_species,
        add_OGs = config["add_OGs"],
        expansion_factor = get_expansion_factor,
        expansion_difference = get_expansion_difference,
	ploidy_normalization = get_ploidy_normalization,
        Nmin_expanded_in = get_Nmin_expanded_in,
        Nmin_compared_to = get_Nmin_compared_to,
        expanded_in_all_found = get_expanded_in_all_found,
        compared_to_all_found = get_compared_to_all_found,
    output:
        #files describing the checkpoint-end needed wildcard "OG"
#        directory("tea/{hypothesis}/expansion_cp_target_OGs/"),
        #for additional OGs instead of BLAST hits
        directory("tea/{hypothesis}/add_OGs_sets/id_lists/"),
        "tea/{hypothesis}/add_OGs_object/add_OG_analysis_object.RDS",
        directory("tea/{hypothesis}/expansion_tibble/"),
        "tea/{hypothesis}/extended_BLAST_hits/extended_BLAST_hits.RDS",
    threads: config["threads_expansion_calc"]
    conda:
        "../envs/expansion.yaml"
    script:
        "../scripts/expansion.R"


#perhaps I can run the Rscript as a rule and only cp the id_lists dirs as .txt files into checkpoint?
#this way even if sth. fails the R part does never need to recompute?
checkpoint expansion:
    input:
#        directory("tea/{hypothesis}/add_OGs_sets/id_lists/"),
#        rules.expansion_computation.output[0],
        id_lists = "tea/{hypothesis}/add_OGs_sets/id_lists/",
        superset_fa = "tea/superset_fasta/superset_species.fa"
    output:
        directory("tea/{hypothesis}/expansion_cp_target_OGs/"),
    conda:
        "../envs/expansion.yaml"
    shell:
        "mkdir -p {output}; "
        "for dir in {input.id_lists}/*; do dir_name=${{dir##*/}}; touch {output}/$dir_name.txt; done"

#########

#since wildcard constrainsts sometimes clash with glob wildcards - see e.g. https://github.com/snakemake/snakemake/issues/482,
#my hacky workaround is to add ".add" to all filenames related to the additionalOGs analysis

## nested snakemake checkpoints are annoying at the moment 
#quick fix - inside the expansion.R script addtional empty set_num files which we can ignore but nonetheless exist are created
#this way, I don't have to worry about evaluating whether or not the user chosen number of addtional OGs were able to be added (for more details see expansion.R script)
#I stick to this trick inside the following rules checkpoint internal rules - if empty file then the biotools aren't used but rather an empty file is created
#in the final_tea script I simply ignore such cases
rule fasta_extraction:
    input:
        protein_lists = "tea/{hypothesis}/add_OGs_sets/id_lists/{OG}/add_OGs_set_num-{set_num}.txt",
#        hypothesis_fasta = "tea/{hypothesis}/ref_fasta/hypothesis_{hypothesis}_species.fa",
        superset_fa = "tea/superset_fasta/superset_species.fa",
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
          faSomeRecords {input.superset_fa} {input.protein_lists} {output}
        fi
        """


rule muscle_MSA:
    input:
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


rule trimAl:
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
          trimal -automated1 -in {input} -out {output};
          #test for edge case in which trimal would produce empty file because everything is trimmed
          if ! [ -s {output} ]; then
            echo "Would trim everything (empty file) - keep raw alignment" &&
            cp {input} {output}
          else
            #do nothing; : is null command in bash
            :
          fi
        fi
        """

rule FastTree:
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


rule expansion_checkpoint_finish:
    input:
        solve_expansion_add_OG_analysis,
    output:
        "checks/expansion/{hypothesis}_finished_add_OG_analysis.txt",
    shell:
        "touch {output}"
