library(readr)
library(stringr)
library(ape)


# read-in hypothesis object
hypotheses <- read_delim("config/hypotheses.tsv", delim = "\t")

# loop through hypotheses
for (hypothesis in hypotheses$hypothesis) {

# assign or subset the orthofinder species_tree based on the current hypothesis
speciesTree <- ape::keep.tip(read.tree("orthofinder/final-results/Species_Tree/SpeciesTree_rooted_node_labels.txt"), 
                                 c(unlist(str_split(hypotheses$expanded_in[hypothesis], ";")), 
                                   unlist(str_split(hypotheses$compared_to[hypothesis], ";"))))

# create directories (handled by snakemake)
    
write.tree(speciesTree, 
           file = paste0("orthofinder/final-results/Species_Tree/", "hypothesis_specific/", hypothesis,
                         "/SpeciesTree_rooted_node_labels.txt"), 
           append = FALSE, 
           tree.names = FALSE)
}
