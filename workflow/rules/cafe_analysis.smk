# CAFE5 analysis
# requires orthofinder to be finished

# install & compile CAFE5
# by far the "unsafest" step in the workflow - if no conda integration in the future; should either be optional or perhaps s part of container..
rule cafe5_setup:
    input:
        ORTHOFINDER + "complete.check"
    output:
        "CAFE5/bin/cafe5"
    conda:
        "../envs/cafe5_compile.yaml"
    shell:
        """
        git clone https://github.com/tgstoecker/CAFE5 CAFE5_cloned/ &&
        mv CAFE5_cloned/* CAFE5/ &&
        rm -rf CAFE5_cloned/ &&
        cd CAFE5/ &&
        autoconf &&
        ./configure &&
        make
        """

# generate an ultrametric species tree from orthofinder output
rule ultrametric_species_tree:
    input:
        "CAFE5/bin/cafe5"
    output:
        "cafe/{hypothesis}/SpeciesTree_rooted_node_labels.txt.ultrametric.tre",
    params:
        species_tree = "orthofinder/final-results/Species_Tree/hypothesis_specific/{hypothesis}/SpeciesTree_rooted_node_labels.txt",
        ultrametric_tree = "orthofinder/final-results/Species_Tree/hypothesis_specific/{hypothesis}/SpeciesTree_rooted_node_labels.txt.ultrametric.tre"
    shell:
        "python workflow/scripts/make_ultrametric.py {params.species_tree} && mv {params.ultrametric_tree} cafe/{wildcards.hypothesis}/"
   

#if more than one species is compared to than we have to split the string based on ";"
#the goal in any scenarios to return a list (!) with just the species used in the hypothesis

def get_all_hypothesis_species_no_path(wildcards):
    """Get compared_to entries from hypotheses(.tsv) for each hypothesis. """
    exp = hypotheses.loc[ (wildcards.hypothesis), 'expanded_in']
    ct = hypotheses.loc[ (wildcards.hypothesis), 'compared_to']
    if ct.count(";") > 0:
        ct = str.split(ct, ";")
        ct.append(exp)       
        return ct
    else:
        output = []
        output.append(exp)
        output.append(ct)
        return output
    return

# also return number of hypothesis
def get_hypo_num(wildcards):
    """Get compared_to entries from hypotheses(.tsv) for each hypothesis. """
    num = hypotheses.loc[ (wildcards.hypothesis), 'hypothesis']
    return num

# includes creation of a genes per HOG per species table
# which is somewhat redundant since we also do this as part of the expansion.R script,
# however doing it also here leads to a parallelisation opportunity
rule reformat_HOG_table:
    input:
        "cafe/{hypothesis}/SpeciesTree_rooted_node_labels.txt.ultrametric.tre"
    output:
        "cafe/{hypothesis}/HOG_table_reformatted_filtered.tsv",
        "cafe/{hypothesis}/HOG_table_reformatted_complete.tsv",
    params:
        all_species = get_all_hypothesis_species_no_path,
        hypothesis_num = get_hypo_num,
    conda:
        "../envs/hypothesis_species_tree.yaml"
    script:
        "../scripts/HOG_table_reformat.R"


# we use CAFE5 twice: from the documentation:
# "Gene families that have large gene copy number variance can cause parameter estimates to be non-informative." 
# "You can remove gene families with large variance from your dataset," 
# "but we found that putting aside the gene families in which one or more species have ≥ 100 gene copies does the trick."
# first run - lambda value computation based on filtered set
rule cafe5_filtered_set:
    input:
        tree = "cafe/{hypothesis}/SpeciesTree_rooted_node_labels.txt.ultrametric.tre",
        table = "cafe/{hypothesis}/HOG_table_reformatted_filtered.tsv",
    output:
        directory("cafe/{hypothesis}/cafe_filtered_results"),
#        lambda_file = "filtered_results/Gamma_results.txt",
    shell:
        "CAFE5/bin/cafe5 --cores 32 -i {input.table} -t {input.tree} -o {output} -k 3"


# custom function to extract lambda value from first run of CAFE5 with filtered set
def get_lambda_value(wildcards):
    path = 'cafe/' + wildcards.hypothesis + '/cafe_filtered_results/Gamma_results.txt'
    # if condition protects from snakemake throwing an error because the file in question does not yet exist
    if not Path(path).exists():
        return -1
    else:
        # Open the file for reading
        with open(path) as fd:
            # Iterate over the lines
            for line in fd:
                # Capture one-or-more characters of non-whitespace after the initial match
                match = re.search(r'Lambda: (\S+)', line)
                # Did we find a match?
                if match:
                    # Yes, process it
                    value = match.group(1)
                    #print('{}'.format(value))
                    return(value)


 # second run - lambda value is used on complete set
rule cafe5_complete_set:
    input:
        tree = "cafe/{hypothesis}/SpeciesTree_rooted_node_labels.txt.ultrametric.tre",
        table = "cafe/{hypothesis}/HOG_table_reformatted_complete.tsv",
#        filtered_results = directory("cafe/cafe_filtered_results"),
        filtered_results = rules.cafe5_filtered_set.output,
    output:
        directory("cafe/{hypothesis}/cafe_complete_results"),
    params:
        lambda_value = get_lambda_value,
    shell:
        "CAFE5/bin/cafe5 --cores 64 -i {input.table} -t {input.tree} -o {output} -k 3 -l {params.lambda_value} -P 0.05"


#rule cafe_check:
#    input:
#        directory("cafe/{hypothesis}/cafe_complete_results"),
#    output:
#        touch("{hypothesis}_cafe_check")
    
