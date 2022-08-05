#backup R-based installation if conda didn't work or wasn't used
#we check if packages are installed first

#fix for: libicui18n.so.68: cannot open shared object file
#reinstall of stringi connects libraries correctly
#install.packages("stringi", repos="http://cran.us.r-project.org")

# list of cran packages
cran_packages = c("readr", "stringr", "ape")
# load or install&load all
package.check <- lapply(
  cran_packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE, repos = "http://cran.us.r-project.org")
      library(x, character.only = TRUE)
    }
  }
)


#load libraries
library(readr)
library(stringr)
library(ape)


# read-in hypothesis object
# we transpose the table since the following code was written for the old layout of hypotheses.tsv
hypotheses <- as.data.frame( 
                t(
                  read.table("config/hypotheses.tsv", 
                             header = FALSE,
                             sep = "\t", 
                             row.names = NULL)
                )              
              )
# first line as header/column names
names(hypotheses) <- hypotheses[1,]
# delete first line
hypotheses <- hypotheses[-1,]
# removal of row.names/numbering
row.names(hypotheses) <- NULL
#correct types
hypotheses$hypothesis <- as.numeric(hypotheses$hypothesis) 
hypotheses$Nmin_expanded_in <- as.numeric(hypotheses$hypothesis)
hypotheses$Nmin_compared_to <- as.numeric(hypotheses$hypothesis)
hypotheses$min_expansion_factor <- as.numeric(hypotheses$hypothesis)


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
