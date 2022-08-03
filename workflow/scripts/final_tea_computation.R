#set logging to file and terminal
log <- file(snakemake@log[[1]], open="wt")
sink(log, append=FALSE, split=FALSE, type="message")


#setting CRAN repository
r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)


## Install packages
# if they are not installed (check via the following neat lapply approach)
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

# list of bioconductor packages
bioc_packages = c("ggtree", "ggtreeExtra", "Biostrings", "DESeq2", "survcomp")

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
cran_packages = c("UpSetR", "cowplot", "ggplotify", "seqinr", "tidyverse", "ape", "stringr", "gtools")
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
library(UpSetR)
library(cowplot)
library(ggplotify)
library(ape)
library(stringr)
library(survcomp)
library(gtools)

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
setClass("expanded_OG", slots=list(blast_table="tbl_df",
                                   add_OG_analysis="list"))


# class for the hypotheses
# adding a prototype is essential here for ortho_intersect_plot since gg has to be defined as 
# register an old-style S3 class using setOldClass
#https://stackoverflow.com/questions/12636056/why-sometimes-i-cant-set-a-class-definition-as-slot-in-a-s4-class
#https://stackoverflow.com/questions/13841400/use-s3-virtual-class-as-slot-of-an-s4-class-got-error-got-class-s4-should-b
#had trouble getting this to work; thus solved differently
# -> seperate list object "Ortho_intersect_plots" hypothesis specific access via index
setClass("hypothesis", 
         slots=list(description="character", 
                                  number="character",
                                  expanded_in ="character", 
                                  compared_to="character", 
                                  expanded_OGs="list",
                                  species_tree="phylo"))

# class for extended BLAST hits info
setClass("extended_BLAST_hits", 
         slots=list(blast_table="tbl_df")
         )

#class for adding OGs analysis
setClass("add_OG_analysis",
         slots=list(add_OG_analysis="list")
         )

#another class for adding OGs analysis
setClass("add_OG_set",
         slots=list(genes="spec_tbl_df",
                    msa="AAStringSet", 
                    tree="phylo"
                   )
         )


#remove protein_names in the snakemake pipeline - directories clean enough
for (hypothesis in hypotheses$hypothesis) {  
   
    #to catch hypotheses making problems
    print("Working on:")
    print(hypothesis)

    # read-in extended_BLAST_hits.RDS object of hypothesis
    extended_BLAST_hits <- readRDS(paste0("tea/", hypothesis, "/extended_BLAST_hits/extended_BLAST_hits.RDS"))

    # add_OG_analysis_object.RDS object of hypothesis
    add_OG_analysis <- readRDS(paste0("tea/", hypothesis, "/add_OGs_object/add_OG_analysis_object.RDS"))
    
    ## adding addtional OGs analysis output as modular attachments
    #This needs to be done before the original steps since we add fasta, msa & trees 
    #to the sets of the add_OGs_analysis object first before adding it as a list object 
    #during the original extended BLAST hits approach

    for (exp_OG in names(add_OG_analysis)) {
      
      #to know faulty OGs
      print("Working on:")
      print(exp_OG)

      for (set in 1:length(add_OG_analysis[[exp_OG]]@add_OG_analysis)) {
      
        #adding msa info
        add_OG_analysis[[exp_OG]]@add_OG_analysis[[set]]@msa <- readAAStringSet(paste0("tea/", hypothesis, "/add_OGs_sets/muscle/", exp_OG, "/add_OGs_set_num-", set, ".afa.add"))
               
        #adding tree info
        add_OG_analysis[[exp_OG]]@add_OG_analysis[[set]]@tree <- read.tree(paste0("tea/", hypothesis, "/add_OGs_sets/trees/", exp_OG, "/add_OGs_set_num-", set, ".tree.add"))

      }
      #to not run into "too many open connections" problem
      closeAllConnections()
    }


    # create empty list object for hypothesis
    assign(paste0("hypothesis_", hypothesis), list())
    # assign list of names
    expanded_OGs <- list.files(path = paste0("tea/", hypothesis, "/expansion_cp_target_OGs/"), 
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
             #genes=read_table(paste0("tea/", hypothesis, "/add_OGs_sets/id_lists/", exp_OG, "/add_OGs_set_num-1.txt"), col_names = FALSE),
             blast_table=extended_BLAST_hits[[exp_OG]]@blast_table,
             #here we add the complete add_OG analysis results (x sets per OG) to the original object
             add_OG_analysis=add_OG_analysis[[exp_OG]]@add_OG_analysis)

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
                  species_tree=speciesTree))
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


#reorder with mixedsort so that "10" comes after 2, etc.
HYPOTHESES.a2tea <- HYPOTHESES.a2tea[mixedsort(names(HYPOTHESES.a2tea))]

print("All done with: HYPOTHESES.a2tea creation")

# final object is called:
# HYPOTHESES.a2tea
#saveRDS(HYPOTHESES.a2tea, "save_HYPOTHESES.a2tea")


## Creating a non-redundant fasta file containing all genes part of the final objects
#n_hypotheses <- length(HYPOTHESES.a2tea)

A2TEA.fa.seqs <- vector(mode = "list")

for (hypothesis in hypotheses$hypothesis) {
    
  print("Now, working on hypothesis:")
  print(hypothesis)

  exp_OGs <- names(HYPOTHESES.a2tea[[hypothesis]]@expanded_OGs)
    
  #add OGs analysis
  for (og in exp_OGs) {

      print("Now, working on OG:")
      print(og)
    
    n_sets <- length(HYPOTHESES.a2tea[[hypothesis]]@expanded_OGs[[og]]@add_OG_analysis)
      
    for (set in 1:n_sets) {
        
      A2TEA.fa.seqs <- c(A2TEA.fa.seqs, read.fasta(paste0("tea/", 
                                hypothesis, 
                                "/add_OGs_sets/fa_records/", 
                                og,
                                "/add_OGs_set_num-", set, ".fa.add"), 
                         seqtype = "AA", 
                         as.string = TRUE))

    }
    #to not run into "too many open connections" problem
    closeAllConnections()
  }  
    
}


A2TEA.fa.seqs <- unique(A2TEA.fa.seqs)

#naming list elements - with id of gene/transcript/protein
for (i in 1:length(A2TEA.fa.seqs)) {
    
  names(A2TEA.fa.seqs)[i] <- attr(A2TEA.fa.seqs[[i]], which = "name")
}


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

#catch edge case where a species posesses zero sig. regulated genes...
#if TRUE add column for this/these species since it won't have been created by the prev. step
all_species_HOG_DE <- HOG_DE.a2tea %>%
  group_by(species) %>%
  summarize(all_species = n()) %>%
  pull(species)

sigDE_species_HOG_DE <- HOG_DE.a2tea %>%
  filter(significant == c("yes")) %>%
  filter(HOG != c("singleton")) %>%
  group_by(species) %>%
  summarize(all_sig_species = n()) %>%
  pull(species)

sig_diff_species_check <- setdiff(all_species_HOG_DE, sigDE_species_HOG_DE)

if (length(sig_diff_species_check) > 0) {
    for (i in sig_diff_species_check) {
      assign("i_mod", paste0(i, "_sigDE"))
      sig_genes_per_species_and_HOG <- sig_genes_per_species_and_HOG %>%
        add_column("{i_mod}" := 0)
    }
} else {
    print("No species without any sig. diff. expressed genes ;D")
}


#handling per hypothesis as hypothesis specific total sig. diff. calculation and join with HOG_level_list 
for (i in 1:length(HOG_level_list)) {
    
    h_expanded_in <- unlist(str_split(hypotheses$expanded_in[i], ";"))
    h_compared_to <- unlist(str_split(hypotheses$compared_to[i], ";"))
    h_species <- c(h_expanded_in, h_compared_to)

    # adding a column summing the rowwise sig. DE counts for all species   
    #need to assign to a hypothesis specific name since this is a loop - duhh...
    h_sig_genes_per_species_and_HOG <- sig_genes_per_species_and_HOG %>% 
      mutate(
        total_sigDE = rowSums(select(., paste0(h_species, "_sigDE")), na.rm=TRUE)  
      )   

    h_sig_genes_per_species_and_HOG <- h_sig_genes_per_species_and_HOG  %>%
      group_by(HOG) %>%
      # mutate all NAs to 0s#
      mutate_at(vars(-group_cols()), ~replace(., is.na(.), 0)) %>%
      # merge rows per HOG - results in one line per HOG
      select(HOG, paste0(h_species, "_sigDE"), total_sigDE) %>%
      distinct() %>%
      summarise_all(list(sum))

    # change the zeros back to NAs
    h_sig_genes_per_species_and_HOG <- h_sig_genes_per_species_and_HOG  %>%
                                       na_if(0)
    
    HOG_level_list[[i]] <- HOG_level_list[[i]] %>% rename_at(vars(-HOG, -expansion), ~ paste0(., '_total'))
    HOG_level_list[[i]] <- full_join(HOG_level_list[[i]], h_sig_genes_per_species_and_HOG, by = c("HOG"))
}


#transforming all NA entries of the HOG tables to 0s - allows for correct filtering in the WebApp
for (i in 1:length(HOG_level_list)) {
    
    HOG_level_list[[i]] <- HOG_level_list[[i]] %>%
                             mutate_all(~replace(., is.na(.), 0))
}

#structure of final list object
str(HOG_level_list)



## Calculating enrichment/overrepresentation of sig. diff. exp. genes of expanded_in species per HOG
#normalization of expanded/compared sig. diff. genes with N species per group
for (i in 1:length(hypotheses$hypothesis)) {
    h_expanded_in <- unlist(str_split(hypotheses$expanded_in[i], ";"))
    h_compared_to <- unlist(str_split(hypotheses$compared_to[i], ";"))
    length_expanded_in <- length(h_expanded_in)
    length_compared_to <- length(h_compared_to)
    
    HOG_level_list[[i]] <- HOG_level_list[[i]] %>%
      mutate(
        norm_sum_expanded_sigDE = rowSums(select(., paste0(h_expanded_in, "_sigDE")), na.rm=TRUE) / length_expanded_in
      )
    
    HOG_level_list[[i]] <- HOG_level_list[[i]] %>%
      mutate(
        norm_sum_compared_sigDE = rowSums(select(., paste0(h_compared_to, "_sigDE")), na.rm=TRUE) / length_compared_to
      )
      
    #compute column sums for norm_sum_expanded_sigDE and norm_sum_compared_sigDE
    #variable nuames correspond to the classic urn problem
    sum_white <- HOG_level_list[[i]] %>% 
      summarise(summed = sum(norm_sum_expanded_sigDE)) %>%
      pull()
    
    sum_black <- HOG_level_list[[i]] %>% 
      summarise(summed = sum(norm_sum_compared_sigDE)) %>%
      pull()
    
    
    #calculate 
    HOG_level_list[[i]] <- HOG_level_list[[i]] %>%
      rowwise() %>%
       mutate(
        oaes_hyper = phyper(norm_sum_expanded_sigDE - 1, 
                             sum_white, 
                             sum_black, 
                             norm_sum_expanded_sigDE + norm_sum_compared_sigDE, 
                             lower.tail = FALSE)
        ) %>% ungroup() %>%      
      mutate(
        oaes_hyper = case_when(
                        oaes_hyper == 1 ~ NA_real_,
                        TRUE ~ oaes_hyper
                   )
      )

}


## adding the computed CAFE p-values to the HOG_level_list(s)
for (hypothesis_num in 1:length(HOG_level_list)) {

  cafe_tibble <- read_tsv(paste0("cafe/", hypothesis_num, "/cafe_complete_results/Gamma_family_results.txt"))[1:2]

  HOG_level_list[[hypothesis_num]] <- HOG_level_list[[hypothesis_num]] %>%
                                        left_join(., 
                                                  cafe_tibble, 
                                                  by = c("HOG" = "#FamilyID")
                                                 ) %>%
                                        rename(cafe_pvalue = pvalue) %>%
                                        relocate(cafe_pvalue, .after = oaes_hyper) %>%
                                        arrange(oaes_hyper)
}

#cafe only outputs p-values until 0.001 - so I will as a quick hack change all 0s to 0.001
#perhaps I can modify the C++ code of CAFE5 a bit so I get longer doubles?
#on the other hand this is probably not really necessary
for (hypothesis_num in 1:length(HOG_level_list)) {

  HOG_level_list[[hypothesis_num]] <- HOG_level_list[[hypothesis_num]] %>%
    mutate(
      cafe_pvalue = case_when(
                      cafe_pvalue == 0 ~ 0.001,
                      TRUE ~ cafe_pvalue
                   )
    )
}


## combine overrep. analysis of sig. exp. genes with cafe value of gene expansion
for (i in 1:length(hypotheses$hypothesis)) {
 
    #calculate the tea-value - trait-associated evolutionary adaption value
    HOG_level_list[[i]] <- HOG_level_list[[i]] %>%
      rowwise() %>%
      mutate(
        tea_value = case_when(
                      is.na(oaes_hyper) ~ NA_real_,
                      is.na(cafe_pvalue) ~ NA_real_,
                      is.na(oaes_hyper) & is.na(cafe_pvalue) ~ NA_real_,
                      TRUE ~ 1 # don't really get why I can't do it here:
                               # combine.test(c(oaes_hyper, cafe_pvalue), na.rm = TRUE)
                    )          # results in error as soon as both elements are NA
                               # seems like case_when calls everything anyways, even though the previous cases
                               # should make it skip the line..; anyways quick workaround with second mutate()
      ) %>%
      mutate(
        tea_value = case_when(
                      tea_value == 1 ~ combine.test(c(oaes_hyper, cafe_pvalue), na.rm = TRUE)
                    )
      ) %>%
      ungroup() %>%
      arrange(tea_value) 
}


#####################################################################
### Addition of HOG info to SFA tibbles
### Create per hypothesis additional functional annotation table on HOG level 
### (H/OG - non-redundant GO terms of all genes of all species in hypothesis)

#we can use the HOG_file_long containing all genes - two columms: HOG; genes
#to add HOG info to all rows in the SFA tables
for (j in 1:length(SFA)) {
  #removal of all singletons? - for now no
  SFA[[j]] <- left_join(SFA[[j]], HOG_file_long, by = c("Protein-Accession" = "gene")) %>% 
    relocate(HOG, .before = "Protein-Accession") #%>%
  }


#######
### ### Adding overview of all species (superset of all hypothesis) table - how many genes per OG + how many sigDE
##-> adding as as LAST element an overview table for all species to HOG_level_list
#this needs to happen AFTER all other HOG_level_list operations since this is not hypothesis specific

# why do we need this: only with such a superset/hypothesis independent table we have the possibility to create arbitrary sets we can compare the hypotheses to
#e.g. with this we can define any conserved set of species we want to compare subset hypotheses against (overrepresentation style analysis like Caro needs)
#by creating this superset, in the WebApp it is then possible to manually define any combination the user might be interested in
#this allows for the freedom necessary, in cases where one species might want/need to be removed without the need of re-running the pipeline, or
#more than one experiment is formulated in the hypotheses.tsv (e.g. truly seperated analysis, that were run in one analysis run)

#get complete orthofinder table
ph_orthogroups <- readr::read_delim(file = "orthofinder/final-results/Phylogenetic_Hierarchical_Orthogroups/N0.tsv",
                          delim = "\t"
                         )


#vector of all unique species names over all hypotheses
all_species <- unique(c(
               unlist(str_split(hypotheses$expanded_in, pattern = ";")), 
               unlist(str_split(hypotheses$compared_to, pattern = ";"))
                     ))

#create dataframe with numbers of genes per PHOG
#we combine expanded_in and compared_to vectors to easily compute for everything we need
all_HOG_df <- setNames(base::data.frame(matrix(ncol = 1 + length(all_species), nrow = 0)),
                  c("HOG", all_species))

for (i in 1:nrow(ph_orthogroups)) {
    
    row = 0
    
    row <- c(ph_orthogroups[i,]$HOG)
    
        for (j in c(all_species)) {
            if (is.na(ph_orthogroups[i,][[j]]))
            {
                test = 0
                row <- c(row, test)
            }
            else 
            {
                test = length(unlist(strsplit(ph_orthogroups[i,][[j]], ",")))
                row <- c(row, test)
            }
}
    all_HOG_df[i,] <- row
} 

all_HOG_tibble <- as_tibble(all_HOG_df)

for (k in 1:length(all_species)) {
    o = all_species[k]
    all_HOG_tibble[[o]] <- as.numeric(all_HOG_tibble[[o]])
}

# append "_gene_count" to all species in all_HOG_tibble
# get number of significantly regulated genes per HOG (per species and total) 
# get count of species, HOG, significant combination from HOG_DE.a2tea
# reduce it to HOG, species, count
sig_genes_per_species_and_HOG <- HOG_DE.a2tea %>%
                                     filter(significant == c("yes")) %>%
                                     filter(HOG != c("singleton")) %>%
                                        group_by(HOG, species, significant) %>%
                                        mutate(count = n()) %>%
                                        ungroup() %>%
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

#catch edge case where a species posesses zero sig. regulated genes...
#if TRUE add column for this/these species since it won't have been created by the prev. step
all_species_HOG_DE <- HOG_DE.a2tea %>%
  group_by(species) %>%
  summarize(all_species = n()) %>%
  pull(species)

sigDE_species_HOG_DE <- HOG_DE.a2tea %>%
  filter(significant == c("yes")) %>%
  filter(HOG != c("singleton")) %>%
  group_by(species) %>%
  summarize(all_sig_species = n()) %>%
  pull(species)

sig_diff_species_check <- setdiff(all_species_HOG_DE, sigDE_species_HOG_DE)

if (length(sig_diff_species_check) > 0) {
    for (i in sig_diff_species_check) {
      assign("i_mod", paste0(i, "_sigDE"))
      sig_genes_per_species_and_HOG <- sig_genes_per_species_and_HOG %>%
        add_column("{i_mod}" := 0)
    }
} else {
    print("No species without any sig. diff. expressed genes ;D")
}



h_sig_genes_per_species_and_HOG <- sig_genes_per_species_and_HOG  %>%
  group_by(HOG) %>%
  # mutate all NAs to 0s#
  mutate_at(vars(-group_cols()), ~replace(., is.na(.), 0)) %>%
  # merge rows per HOG - results in one line per HOG
  select(HOG, paste0(all_species, "_sigDE")) %>%
  distinct() %>%
  summarise_all(list(sum))

  
# change the zeros back to NAs
h_sig_genes_per_species_and_HOG <- h_sig_genes_per_species_and_HOG  %>%
                                       na_if(0)
#
all_HOG_tibble <- all_HOG_tibble %>% rename_at(vars(-HOG), ~ paste0(., '_total'))
all_HOG_tibble <- full_join(all_HOG_tibble, h_sig_genes_per_species_and_HOG, by = c("HOG"))
    
#substitute remaining NAs for 0s
all_HOG_tibble <- all_HOG_tibble %>% mutate_all(~replace(., is.na(.), 0))

#make it a list
all_species_overview <- list(all_HOG_tibble)

#name it a list
names(all_species_overview) <- "all_species_overview"

#add to HOG_level_list as last list element
HOG_level_list <- c(HOG_level_list, all_species_overview)


### add complete species tree to .RData since this is actually more sensible for the WebApp than the subset trees
all_speciesTree <- read.tree("orthofinder/final-results/Species_Tree/SpeciesTree_rooted_node_labels.txt")


#########################################################################################################

## Last step: saving everything to one file which is input for the A2TEA WebApp
save(hypotheses, 
     HYPOTHESES.a2tea, 
     A2TEA.fa.seqs,
     HOG_DE.a2tea, 
     HOG_level_list,
     SFA,
     all_speciesTree,
     file = "tea/A2TEA_finished.RData",
     compress = "xz", 
     compression_level = 9)
