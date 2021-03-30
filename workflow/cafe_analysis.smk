# CAFE5 analysis
# requires orthofinder to be finished

# generate an ultrametric species tree from orthofinder output
rule ultrametric_species_tree:
    input:
        ORTHOFINDER + "complete.check"
    output:
        "cafe/SpeciesTree_rooted_node_labels.txt.ultrametric.tre"
    params:
        species_tree = "orthofinder/final-results/Species_Tree/SpeciesTree_rooted_node_labels.txt",
        ultrametric_tree = "orthofinder/final-results/Species_Tree/SpeciesTree_rooted_node_labels.txt.ultrametric.tre"
    shell:
        "python scripts/make_ultrametric.py {params.species_tree} && mv {params.ultrametric_tree} cafe/"
   

# includes creation of a genes per HOG per species table
# which is somewhat redundant since we also do this as part of the expansion.R script,
# however doing it also here leads to a parallelisation opportunity
rule reformat_HOG_table:
    input:
        "cafe/SpeciesTree_rooted_node_labels.txt.ultrametric.tre"
    output:
        "cafe/HOG_table_reformatted_filtered.tsv",
        "cafe/HOG_table_reformatted_complete.tsv",
    params:
        all_species = SPECIES,
    script:
        "scripts/HOG_table_reformat.R"


# we use CAFE5 twice: from the documentation:
# "Gene families that have large gene copy number variance can cause parameter estimates to be non-informative." 
# "You can remove gene families with large variance from your dataset," 
# "but we found that putting aside the gene families in which one or more species have ≥ 100 gene copies does the trick."
# first run - lambda value computation based on filtered set
rule cafe5_filtered_set:
    input:
        tree = "cafe/SpeciesTree_rooted_node_labels.txt.ultrametric.tre",
        table = "cafe/HOG_table_reformatted_filtered.tsv",
    output:
        directory("cafe/cafe_filtered_results"),
#        lambda_file = "filtered_results/Gamma_results.txt",
    shell:
        "CAFE5/bin/cafe5 --cores 64 -i {input.table} -t {input.tree} -o {output} -k 3"


# custom function to extract lambda value from first run of CAFE5 with filtered set
def get_lambda_value(wildcards):
    # if condition protects from snakemake throwing an error because the file in question does not yet exist
    if not Path('cafe/cafe_filtered_results/Gamma_results.txt').exists():
        return -1
    else:
        # Open the file for reading
        with open('cafe/cafe_filtered_results/Gamma_results.txt') as fd:
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
        tree = "cafe/SpeciesTree_rooted_node_labels.txt.ultrametric.tre",
        table = "cafe/HOG_table_reformatted_complete.tsv",
#        filtered_results = directory("cafe/cafe_filtered_results"),
        filtered_results = rules.cafe5_filtered_set.output,
    output:
        directory("cafe/cafe_complete_results"),
    params:
        lambda_value = get_lambda_value,
    shell:
        "CAFE5/bin/cafe5 --cores 64 -i {input.table} -t {input.tree} -o {output} -k 3 -l {params.lambda_value} -P 0.05"
