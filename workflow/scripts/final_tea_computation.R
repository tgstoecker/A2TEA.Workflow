## Install packages
# if they are not installed (check via the following neat lapply approach)
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

# list of bioconductor packages
bioc_packages = c("ggtree", "ggtreeExtra", "Biostrings", "DESeq2")

# load or install&load all
package.check <- lapply(
  bioc_packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
    BiocManager::install(x)
    library(x, character.only = TRUE)
    }
  }
)

# list of cran packages
cran_packages = c("UpSetR", "cowplot", "ggplotify", "seqinr", "tidyverse", "ape")
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


## Load all necessary packages
library(DESeq2)
#library(BiocGenerics) # comes with DESeq2
#library(apeglm) # currently difficulties with install # better DESseq2 calc; later 
library(tidyverse)
library(ggtree)
library(Biostrings)
library(seqinr)
#library(VennDiagram)
library(UpSetR)
library(cowplot)
library(ggplotify)
library(ape)


# get user choice for DEG FDR cutoff value
DEG_FDR = snakemake@params["DEG_FDR"]

#get functional annotation list object
SFA = readRDS(snakemake@input[["SFA"]])

## Find all DESeq2 differential expression RDS objects and load them
### Names are given based on the species/ecotype/etc. name
dea_list <- list.files(path = "R/deseq2/dea_final/", pattern = "dea_*", full.names=TRUE)
dea_list_short <- list.files(path = "R/deseq2/dea_final/", pattern = "dea_*", full.names=FALSE)

listRDS <- lapply(dea_list, readRDS, .GlobalEnv)

for (i in 1:length(listRDS)){
    assign(str_sub(dea_list_short[[i]], start=5), listRDS[[i]])
}

# remove unused combined list object
rm(listRDS)


## create dataframe of species and diff. exp. results
species_list <- vector()

for (i in 1:length(dea_list_short)){
    species_list <- c(species_list, str_sub(dea_list_short[[i]], start=5))
}

list_AllSpeciesDEResultsDataFrames <- list()

for (i in species_list) {
    # individual species file could be deleted but perhaps I'll need them later; kept for now...
    assign(paste0(i, "_DEresultsTable"), as.data.frame(results(get(i))))
    # create gene name vector here (when list is constructed rownames become unique species-gene combination)
    gene <- rownames(get(paste0(i, "_DEresultsTable")))
    df_with_genes <- add_column(get(paste0(i, "_DEresultsTable")), gene, .before = "baseMean")
    # also add a species column which will will also come in handy during the shiny steps
    species <- replicate(nrow(df_with_genes), i)
    df_with_genes <- add_column(df_with_genes, species, .before = "gene")
    # create list of dataframes, which will come in handy
    list_AllSpeciesDEResultsDataFrames[[i]] <- df_with_genes
}
#create combined "mega" dataframe of all species, which is going  to be used for shiny lookup
combined_AllSpeciesDEResultsDataFrames <- do.call("rbind", list_AllSpeciesDEResultsDataFrames)


## Create long format HOG-genes relation table
HOG_file_raw <- read_delim("orthofinder/final-results/Phylogenetic_Hierarchical_Orthogroups/N0.tsv", delim = "\t")

#drop the OG, Gene Tree Parent Clade
HOG_file_short <- HOG_file_raw %>% select (-c(OG, `Gene Tree Parent Clade`)) 


#merge all species columns into one and remove the solo ones
#use all colnames except first - thus all species
HOG_file_merged <- HOG_file_short %>% unite(x,
                                          c(colnames(HOG_file_short)[-1]),
                                          sep = ",", 
                                          na.rm = TRUE,
                                          remove = TRUE)

#transform into long format - each row 1 gene and it's corresponding HOG in seperate columns
HOG_file_long <- HOG_file_merged %>%
                    mutate(unpacked = str_split(x, ",")) %>%
                    unnest(cols = c(unpacked)) %>%
                    mutate(genes = str_trim(unpacked)) %>% 
                    select(-c(x, unpacked)) %>% 
                    rename(gene = genes)


## Adding specific HOG or singleton info as HOG column to DE tables
HOG_DE.a2tea <- full_join(combined_AllSpeciesDEResultsDataFrames, HOG_file_long, 
              by = c("gene")) %>% 
                replace_na(list(HOG = "singleton"))

#add column for significance - level is set by user in config.yaml ["DEG_FDR"]
#subset tidyverse with the right functions ;D
significant <- vector()
for (FDR in HOG_DE.a2tea$padj) {
    if (!is.na(FDR) && FDR < DEG_FDR) {
        significant <- c(significant, "yes")
    }
    else if (!is.na(FDR) && FDR > DEG_FDR) {
        significant <- c(significant, "no")
    }
    else if (is.na(FDR)) {
        significant <- c(significant, "no")
    }
}
HOG_DE.a2tea <- add_column(HOG_DE.a2tea, significant, .after = "padj")


## Load in the rest of the data - hypotheses, trees, fasta and msa for the start 
### VennDiagrammes included
### + toDo general stats, especially once Orthofinder calculates for HOG)
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
hypotheses$Nmin_expanded_in <- as.numeric(hypotheses$Nmin_expanded_in)
hypotheses$Nmin_compared_to <- as.numeric(hypotheses$Nmin_compared_to)
hypotheses$min_expansion_factor <- as.numeric(hypotheses$min_expansion_factor)

# create hypotheses object
# each object has list of exp. OGs
# iterate over list of hyptheses and associate expansions, fasta, msa and trees


# first goal - hypothesis object with hypothesis and all exp.HOGs
# parsing of hyptheses.tsv
hypotheses$hypothesis


# define three classes
# class for the expanded_OG - containing all different types of data we have on it
setClass("expanded_OG", slots=list(genes="spec_tbl_df",
                                   blast_table="tbl_df",
                                   nrow_table="numeric",
                                   num_genes_HOG="numeric",
                                   num_genes_extend="numeric",
                                   num_genes_complete="numeric",
                                   genes_HOG="tbl_df",
                                   genes_extend_hits="tbl_df",
                                   fasta_files="list", 
                                   msa="AAStringSet", 
                                   tree="phylo"))


# class for the hypotheses
# adding a prototype is essential here for ortho_intersect_plot since gg has to be defined as 
# register an old-style S3 class using setOldClass
#https://stackoverflow.com/questions/12636056/why-sometimes-i-cant-set-a-class-definition-as-slot-in-a-s4-class
#https://stackoverflow.com/questions/13841400/use-s3-virtual-class-as-slot-of-an-s4-class-got-error-got-class-s4-should-b
#had trouble getting this to work; thus solved differently
# -> seperate list object "Ortho_intersect_plots" hypothesis specific access via index
setClass("hypothesis", 
#         prototype=prototype(ortho_intersect_plot=structure(list(), class="gList")),
         slots=list(description="character", 
                                  number="character",
                                  expanded_in ="character", 
                                  compared_to="character", 
                                  expanded_OGs="list",
                                  species_tree="phylo"))
#                                  ortho_intersect_plot="gg"))

# class for extended BLAST hits info
setClass("extended_BLAST_hits", 
         slots=list(blast_table="tbl_df",
                    num_genes_HOG="numeric",
                    num_genes_extend="numeric",
                    num_genes_complete="numeric",
                    genes_HOG="tbl_df",
                    genes_extend_hits="tbl_df")
         )

#remove protein_names in the snakemake pipeline - directories clean enough
for (hypothesis in hypotheses$hypothesis) {  
    # read-in extended_BLAST_hits.RDS object of hypothesis
    extended_BLAST_hits <- readRDS(paste0("tea/", hypothesis, "/extended_BLAST_hits/extended_BLAST_hits.RDS"))

    # create empty list object for hypothesis
    assign(paste0("hypothesis_", hypothesis), list())
    # assign list of names
    expanded_OGs <- list.files(path = paste0("tea/", hypothesis, "/exp_OGs_proteinnames/"), 
                               pattern = "*", 
                               full.names=FALSE)
    expanded_OGs_short <- str_sub(expanded_OGs, end=-5)
    # assign or subset the orthofinder species_tree based on the current hypothesis
    speciesTree <- ape::keep.tip(read.tree("orthofinder/final-results/Species_Tree/SpeciesTree_rooted_node_labels.txt"), 
                                 c(unlist(str_split(hypotheses$expanded_in[hypothesis], ";")), 
                                   unlist(str_split(hypotheses$compared_to[hypothesis], ";"))))    
    #create empty list
    t <- list()

    for (exp_OG in expanded_OGs_short) {
        test <- new("expanded_OG", 
             genes=read_table(paste0("tea/", hypothesis, "/exp_OGs_proteinnames/", exp_OG, ".txt"), col_names = FALSE),
             blast_table=extended_BLAST_hits[[exp_OG]]@blast_table,
             num_genes_HOG=extended_BLAST_hits[[exp_OG]]@num_genes_HOG,
             num_genes_extend=extended_BLAST_hits[[exp_OG]]@num_genes_extend,
             num_genes_complete=extended_BLAST_hits[[exp_OG]]@num_genes_complete,
             genes_HOG=extended_BLAST_hits[[exp_OG]]@genes_HOG,
             genes_extend_hits=extended_BLAST_hits[[exp_OG]]@genes_extend_hits,
             fasta_files=read.fasta(paste0("tea/", hypothesis, "/fa_records/", exp_OG,".fa"), seqtype = "AA", as.string = TRUE), 
             msa=readAAStringSet(paste0("tea/", hypothesis, "/muscle/", exp_OG, ".afa")), 
             tree=read.tree(paste0("tea/", hypothesis, "/trees/", exp_OG, ".tree")))
        x <- list(test)
        names(x) <- paste0(exp_OG)
        t <- c(t, x)
    }
 
    # create for each hypothesis a complete hypothesis object with correct naming
    # order: hypothesis@info (name, number, expand, compared, exp OGs)$indiv. exp_OGs@objects (exp.OG list, fa, msa, tree)
    # "expanded_in" and "compared_to" are split on ";" in case their size is > 1
    assign(paste0("hypothesis_", hypothesis),
              new("hypothesis", 
                  description=hypotheses$name[hypothesis],
                  number=as.character(hypothesis),
                  expanded_in=unlist(str_split(hypotheses$expanded_in[hypothesis], ";")),
                  compared_to=unlist(str_split(hypotheses$compared_to[hypothesis], ";")), 
                  expanded_OGs=t,
                  species_tree=speciesTree))#,
#                  ortho_intersect_plot=NULL))
}


# create empty list for final complete hypothesis object
# create name list of all hypotheses
hypotheses_list <- ls(pattern = "hypothesis_")

# create empty list object to completely hold all hypotheses and associated data
HYPOTHESES.a2tea <- list()

for (hypothesis in hypotheses_list) {
    h <- list(get(hypothesis))
    names(h) <- paste0(hypothesis)
    HYPOTHESES.a2tea <- c(HYPOTHESES.a2tea, h)
}

rm(hypothesis_1)
rm(hypothesis_2)

# final object is called:
# HYPOTHESES.a2tea


## Creation of HOG level tables 
###  needs to be behind creation of hypotheses, because they are used here

# create hypothesis specfic HOG level file

HOG_level_list <- list()

for (i in 1:length(hypotheses$hypothesis)) {
    assign(paste0("hypothesis_", i, "_expansion_tibble"),
           readRDS(
               list.files(path = paste0("tea/", i, "/expansion_tibble"), 
                      pattern = "expansion_tibble.rds", 
                      full.names=TRUE)
                  )
           )
    # add to list and name list entries according to hypothesis
    hog_level <- list(get(paste0("hypothesis_", i, "_expansion_tibble")))
    names(hog_level) <- paste0("hypothesis_", i)
    HOG_level_list <- c(HOG_level_list, hog_level)
}

# e.g. the following displays the tibble
# hypothesis_2_expansion_tibble

# append "_gene_count" to all species in HOG tables
# get number of significantly regulated genes per HOG (per species and total) 
# get count of species, HOG, significant combination from HOG_DE.a2tea
# reduce it to HOG, species, count
sig_genes_per_species_and_HOG <- HOG_DE.a2tea %>% 
                                     filter(significant == c("yes")) %>%
                                     filter(HOG != c("singleton")) %>%
                                        group_by(HOG, species, significant) %>% 
                                        mutate(count = n()) %>%
                                        ungroup() %>%
                                        select(HOG, species, count) %>%
            # https://www.r-bloggers.com/2018/05/workaround-for-tidyrspread-with-duplicate-row-identifiers/
            # spread error when no indexing for data
                                            group_by(species) %>% 
                                            mutate(grouped_id = row_number()) %>%  
                                            spread(species, count) %>% 
                                            select(-grouped_id) %>%
            # easy workaround for duplicated rows
                                                distinct()%>% 
            # rename species columns containing now the counts of sig. DE genes
                                                    rename_at(vars(-HOG), ~ paste0(., '_sigDE'))
    # and also add a column summing the rowwise sig. DE counts for all species
    total_sigDE <- sig_genes_per_species_and_HOG %>% select(-HOG) %>% rowSums(na.rm = TRUE) 
    sig_genes_per_species_and_HOG <- add_column(sig_genes_per_species_and_HOG, total_sigDE)
    sig_genes_per_species_and_HOG <- sig_genes_per_species_and_HOG  %>%
                                         group_by(HOG) %>%
                                         # mutate all NAs to 0s
                                         mutate_at(vars(-group_cols()), ~replace(., is.na(.), 0)) %>%
                                         # merge rows per HOG - results in one line per HOG
                                         #summarise_all(funs(sum))
                                         summarise_all(list(sum))
# change the zeros back to NAs
# necessary for current implemntation of tea value computation
    sig_genes_per_species_and_HOG <- sig_genes_per_species_and_HOG  %>%
                                         na_if(0)

for (i in 1:length(HOG_level_list)) {
    HOG_level_list[[i]] <- HOG_level_list[[i]] %>% rename_at(vars(-HOG, -expansion), ~ paste0(., '_total'))
    HOG_level_list[[i]] <- full_join(HOG_level_list[[i]], sig_genes_per_species_and_HOG, by = c("HOG"))
}



##################
#### tea-value 
## after tests & talking with Heiko an enrichment Fisher approach might be better suited to tackle this problem..
## anyways in the following; commented out and not working with the newest code are the initial formula based ideas for the tea-value

##################################################################

## computation of tea value (first idea)
## expansion ratio (exp/com) (9/4) x expression ratio (com/exp) (2/6) 
## ---- don't need (n_species_expanded/n_species_compared); stays equal in a hypothesis; 
## similar to p-values not comparable between different experiments
## divided by HOG size (all genes in HOG for this particular hypothesis) x number of significantly regulated genes
## problem - how to deal with missing values???

# should I only have hypothesis species (total count and sigDE) HOG_level_list tables?

# new column tea-value
# if expansion yes continue; else NA value
# because we set expansion as criterium, species gene count has to be >=1
# sig DE counts can be NA ! - we could just add a pseudocount of 1 at these positions for the calculation

# so per line that has expansion = yes:
# get all expanded species; sum their gene counts; sum their sigDE counts
# get all compared species; sum their gene counts; sum their sigDE counts
# expansion_ratio <- divide summed gene count expanded by summed gene count compared
# expression_ratio <- divide summed SigDE count expanded by summed SigDE count compared
# e2_ratio <- divide expansion_ratio by expression_ratio
# tea-value <- divide e2_ratio by total number of genes in HOG and then by total number of SigDE in HOG

## outer most loop; going through hypotheses in HOG_level_list
#for (i in 1:length(HOG_level_list)) {
#    
#    # access the names of all species important for the hypothesis - expanded_in; compared_to
#    print(c(unlist(str_split(HYPOTHESES.a2tea[[i]]@expanded_in, ";")), 
#            unlist(str_split(HYPOTHESES.a2tea[[i]]@compared_to, ";"))))
#    
#    # creating a relevant species vector might make this more future proof;
#    # since I might change the completeness of the HOG_level_list elements
##    relevant_species <- c(unlist(str_split(HYPOTHESES.a2tea[[i]]@expanded_in, ";")), 
##                          unlist(str_split(HYPOTHESES.a2tea[[i]]@compared_to, ";")))
#    
#    expanded_species <- unlist(str_split(HYPOTHESES.a2tea[[i]]@expanded_in, ";"))
#    
#    compared_species <- unlist(str_split(HYPOTHESES.a2tea[[i]]@compared_to, ";"))                              
#                                  
#    ## pre-allocate empty vector for tea value with length of nrow of current hypothesis of HOG_level_list
#    # neat; because NAs already in there so I can just skip if the HOG shows no expansion ;D
#    tea_value <- rep(NA, nrow(HOG_level_list[[i]]))
#    
#    
#    ## create "geneCount_workset" containing HOGs and counts of all relevant species
#    geneCount_workset <-  HOG_level_list[[i]]
#
#    # create geneCount_EXP_workset & geneCount_COMP_workset + drop columns that are unnecessary                  
#    geneCount_EXP_workset <- geneCount_workset %>% 
#                             select(ends_with("_total")) %>%
#                             select(contains(expanded_species))
#                        # also remove the "_total" ending
#    geneCount_EXP_workset <- geneCount_EXP_workset %>% 
#                             setNames(names(geneCount_EXP_workset) %>% 
#                             stringr::str_replace("_total",""))
#                          
#    geneCount_COMP_workset <- geneCount_workset %>% 
#                             select(ends_with("_total")) %>%
#                             select(contains(compared_species))
#                        # also remove the "_total" ending
#    geneCount_COMP_workset <- geneCount_COMP_workset %>% 
#                             setNames(names(geneCount_COMP_workset) %>% 
#                             stringr::str_replace("_total",""))
#                    
#
#    ## create "geneSig_workset" containing HOGs and counts of all relevant species 
#    geneSig_workset <-  HOG_level_list[[i]]
#    
#    # create geneSig_EXP_workset & geneSig_COMP_workset + drop columns that are unnecessary                  
#    geneSig_EXP_workset <- geneSig_workset %>% 
#                             select(ends_with("_sigDE")) %>%
#                             select(contains(expanded_species))
#                        # also remove the "_sigDE" ending
#    geneSig_EXP_workset <- geneSig_EXP_workset %>% 
#                             setNames(names(geneSig_EXP_workset) %>% 
#                             stringr::str_replace("_sigDE",""))
#                          
#    geneSig_COMP_workset <- geneSig_workset %>% 
#                             select(ends_with("_sigDE")) %>%
#                             select(contains(compared_species))
#                        # also remove the "_sigDE" ending
#    geneSig_COMP_workset <- geneSig_COMP_workset %>% 
#                             setNames(names(geneSig_COMP_workset) %>% 
#                             stringr::str_replace("_sigDE",""))
#
#                              
#    # nested loop that goes through current hypothesis HOG_level_list line by line
#    for (j in 1:nrow(HOG_level_list[[i]])) {
#        # write value into index slot if current HOG does pass hard filter expansion criterium
#        # if expansion == "no" do nothing; since we prefilled all positions with NA
#        if (HOG_level_list[[i]][j,"expansion"] == "yes") {
#            
#            ## compute all the different sub-values for the tea_value
#            # count of all genes in expanded species in HOG
#            geneCount_EXP_HOG <- geneCount_EXP_workset[j,] %>%
#                                 dplyr::rowwise() %>%
#                                 sum()
#            
#            # count of all genes in compared species in HOG
#            geneCount_COMP_HOG <- geneCount_COMP_workset[j,] %>%
#                                      dplyr::rowwise() %>%
#                                      sum()
#            
#            # sum genes of expanded & compared species in HOG
#            geneCount_BOTH_HOG <- sum(geneCount_EXP_HOG, geneCount_COMP_HOG)
#            
#            
#            # count of all sig DE genes in expanded species in HOG
#            geneSig_EXP_HOG <- geneSig_EXP_workset[j,] %>%
#                                 dplyr::rowwise() %>%
#                                 sum()
#            
#            # count of all sig DE genes in compared species in HOG
#            geneSig_COMP_HOG <- geneSig_COMP_workset[j,] %>%
#                                      dplyr::rowwise() %>%
#                                      sum()
#            
#            ## sigDE columns can contain zero (NA in table) counts
#            ## we exchange this for a pseudocount 1
#            ## doesn't impact tea_value too much
#            if (is.na(geneSig_EXP_HOG)) {
#                geneSig_EXP_HOG <- 1
#            }
#            
#            if (is.na(geneSig_COMP_HOG)) {
#                geneSig_COMP_HOG <- 1
#            }
#            
#            # sum Sig DE genes of expanded & compared species in HOG
#            geneSig_BOTH_HOG <- sum(geneSig_EXP_HOG, geneSig_COMP_HOG)
#
#            
#            ## perform tea value calculation
#            int_tea_value <- (geneCount_EXP_HOG / geneCount_COMP_HOG) * (geneSig_COMP_HOG / geneSig_EXP_HOG) / 
#                                                  (geneCount_BOTH_HOG * geneSig_BOTH_HOG)
#            
#            # lastly, add computed value to tea_value at index of current j loop position
#            tea_value[j] <- int_tea_value
#        }
#    }
#    
#    ## add tea_value as new column to HOG_level_list table
#    HOG_level_list[[i]] <- HOG_level_list[[i]] %>% 
#                               add_column(tea_value, .after = "HOG")
#    
#}

######################################################################


## adding the computed CAFE p-values to the HOG_level_list(s)
for (hypothesis_num in 1:length(HOG_level_list)) {

  cafe_tibble <- read_tsv(paste0("cafe/", hypothesis_num, "/cafe_complete_results/Gamma_family_results.txt"))[1:2]

  HOG_level_list[[hypothesis_num]] <- HOG_level_list[[hypothesis_num]] %>%
                                        left_join(., 
                                                  cafe_tibble, 
                                                  by = c("HOG" = "#FamilyID")
                                                 ) %>%
                                        rename(cafe_pvalue = pvalue) %>%
                                        relocate(cafe_pvalue, .after = "HOG") %>%
                                        arrange(cafe_pvalue)
}


## Last step: saving everything to one file which is input for the A2TEA WebApp
save(hypotheses, 
     HYPOTHESES.a2tea, 
     HOG_DE.a2tea, 
     HOG_level_list,
     SFA,
     file = "tea/A2TEA_finished.RData")
