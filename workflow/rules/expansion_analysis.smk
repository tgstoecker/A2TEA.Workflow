#Expansion Analysis
# can handle both multi expansion and multi comparison species, now ;D
# the strplit command works really well in this context
## this way, multiple species are propagated into R as a simple vector  

def get_hypo_num(wildcards):
    """Get compared_to entries from hypotheses(.tsv) for each hypothesis. """
    num = hypotheses.loc[ (wildcards.hypothesis), 'hypothesis']
    return num

def get_hypo_name(wildcards):
    """Get compared_to entries from hypotheses(.tsv) for each hypothesis. """
    name = hypotheses.loc[ (wildcards.hypothesis), 'name']
    return name


# for both the species we want to check for expansion as well as those being compared to:
# if more than one species is compared to than we have to split the string based on ";" 

def get_exp_species(wildcards):
    """Get compared_to entries from hypotheses(.tsv) for each hypothesis. """
    exps_1 = hypotheses.loc[ (wildcards.hypothesis), 'expanded_in']
    if exps_1.count(";") > 0:
        exps_2 = str.split(exps_1, ";")
        return exps_2
    else:
        return exps_1
    return



#if more than one species is compared to than we have to split the string based on ";"
def get_com_species(wildcards):
    """Get compared_to entries from hypotheses(.tsv) for each hypothesis. """
    ct_1 = hypotheses.loc[ (wildcards.hypothesis), 'compared_to']
    if ct_1.count(";") > 0:
        ct_2 = str.split(ct_1, ";")
        return ct_2
    else:
        return ct_1
    return


#if more than one species is compared to than we have to split the string based on ";"
#the goal in any scenarios to return a list (!) with just the species used in the hypothesis
#we also need to add the path to the longest isoform peptide fasta files and add the .fa suffix
#another trick is that we can simply use the species base name in any case,
#since the longest isoform output is always named after it!
# this also applies if the user did the isoform filtering on their own since the renamed pep fastas are linked to isoform directory in that case ;D
#+the filtering for longest isoform also only retains the base name of the gene/protein
def get_all_hypothesis_species(wildcards):
    """Get compared_to entries from hypotheses(.tsv) for each hypothesis. """
    path_prefix = 'FS/longest_isoforms/'
    suffix = '.fa'
    exp = hypotheses.loc[ (wildcards.hypothesis), 'expanded_in']
    ct = hypotheses.loc[ (wildcards.hypothesis), 'compared_to']
    if ct.count(";") > 0:
        ct = str.split(ct, ";")
        ct.append(exp)
        ct = [path_prefix + x + suffix for x in ct]
        return ct
    else:
        output = []
        output.append(exp)
        output.append(ct)
        output = [path_prefix + x + suffix for x in output]
        return output
    return


rule create_hypothesis_fasta:
    input:
        orthology = ORTHOFINDER + "complete.check",
    output:
        "tea/{hypothesis}/ref_fasta/hypothesis_{hypothesis}_species.fa",
    params:
        all_hyp_species = get_all_hypothesis_species,
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
        add_blast_hits = config["add_blast_hits"]
    output:
        directory("tea/{hypothesis}/exp_OGs_proteinnames/"),
        directory("checks/tea/{hypothesis}/"),
        directory("tea/{hypothesis}/expansion_tibble/"),
        "tea/{hypothesis}/extended_BLAST_hits/extended_BLAST_hits.RDS",
    threads: 1
    script:
        "../scripts/expansion.R"


rule fasta_extraction:
    input:
        protein_lists = "tea/{hypothesis}/exp_OGs_proteinnames/{OG}.txt",
        hypothesis_fasta = "tea/{hypothesis}/ref_fasta/hypothesis_{hypothesis}_species.fa",
    output:
        "tea/{hypothesis}/fa_records/{OG}.fa"
    threads: 1
    shell:
        "faSomeRecords {input.hypothesis_fasta} {input.protein_lists} {output}"


rule muscle_MSA:
    input:
        "tea/{hypothesis}/fa_records/{OG}.fa"
    output:
        "tea/{hypothesis}/muscle/{OG}.afa"
    threads: 1
    shell:
        "muscle -in {input} -out {output}"



rule trimAl:
    input:
        "tea/{hypothesis}/muscle/{OG}.afa"
    output:
        "tea/{hypothesis}/trimAl/{OG}.afa"
    threads: 1
    shell:
        "trimal -automated1 -in {input} -out {output}"


rule FastTree:
    input:
        "tea/{hypothesis}/trimAl/{OG}.afa"
    output:
        "tea/{hypothesis}/trees/{OG}.tree"
    threads: 1
    shell:
        "FastTree {input} > {output}"


#CHECKPOINTS ARE SO AWESOME!
def solve_expansion(wildcards):
    checkpoint_output = checkpoints.expansion.get(**wildcards).output[0]
    file_names = expand("tea/{hypothesis}/trees/{OG}.tree", hypothesis=wildcards.hypothesis, OG=glob_wildcards(os.path.join(checkpoint_output, "{OG}.txt")).OG)
    return file_names


rule expansion_checkpoint_finish:
    input:
        solve_expansion
    output:
        "checks/expansion/{hypothesis}_finished.txt",
    shell:
        "touch {output}"

# necessary to be done here with the expression analysis;
#        expression = expand("R/deseq2/dea_final/dea_{species}", species=SPECIES),
