#to reduce installation time; we could reduce down to readr, plyr, dplyr, stringr
#library(tidyverse)
library(readr)
library(plyr)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)

message("Acquiring hypothesis variables:")
num = snakemake@params[["num"]]
name = snakemake@params[["name"]]

# define the number of additional best blast hits to include in the follow-up analyses
add_blast_hits = snakemake@params[["add_blast_hits"]]

# define the required expansion factor (part of hypothesis.tsv)
expansion_factor = as.numeric(snakemake@params[["expansion_factor"]])

# get ploidy information of each species (supplied by user; part of species.tsv)
species_table <- read.delim("config/species.tsv", header = TRUE, sep = "\t", row.names = "species")

# using snakemake propagation + python strsplit() is really cool since the input is just a vector
## even if we have multiple species in expanded, compared or both ;D
expanded_in = snakemake@params[["expansion"]]
compared_to = snakemake@params[["comparison"]]
c_t_species <- compared_to
all_species <- unique(c(expanded_in, compared_to))

# defining:
# at least Nmin_expanded_in expanded species that are expanded in at least Nmin_compared_to compared_to species
Nmin_expanded_in <- as.numeric(snakemake@params[["Nmin_expanded_in"]])
Nmin_compared_to <- as.numeric(snakemake@params[["Nmin_compared_to"]])

# define whether or not only HOGs should be considered that have at least 1 gene from each expanded_in/compared_to species
expanded_in_all_found <- as.character(snakemake@params[["expanded_in_all_found"]])
compared_to_all_found <- as.character(snakemake@params[["compared_to_all_found"]])

#read-in Orthogroup-GeneCounts-table
message("Reading in Orthogroup-GeneCounts-table:")
OG.GC <- readr::read_tsv("orthofinder/final-results/Orthogroups/Orthogroups.GeneCount.tsv")


#concatenate all BLAST seaches to one file and create parseable tibble with read_delim
#could add that only BLAST searches of our species in this hypothesis are used - toDO for later
#easiest and efficient way is to bind outside of the loop
#do.call is a nice way rbind(bind_rows[[1]], bind_rows[[2]], ...)
message("Concatenating all BLAST seaches to one file and creating parseable tibble:")
datalist = list()
for (i in Sys.glob("orthofinder/search/*.txt")) {
    datalist[[i]] <- readr::read_delim(file = i,
                                delim = "\t",
                                col_names = c(
                                    "qseqid",
                                    "sseqid",
                                    "pident",
                                    "length",
                                    "mismatch",
                                    "gapopen",
                                    "qstart",
                                    "qend",
                                    "sstart",
                                    "send",
                                    "evalue",
                                    "bitscore"),
                                col_types = c(
                                    qseqid = col_character(),
                                    sseqid = col_character(),
                                    pident = col_double(),
                                    length = col_double(),
                                    mismatch = col_double(),
                                    gapopen = col_double(),
                                    qstart = col_double(),
                                    qend = col_double(),
                                    sstart = col_double(),
                                    send = col_double(),
                                    evalue = col_double(),
                                    bitscore = col_double()
                                                )
                                )
                                                }

all_BLAST <- base::do.call(bind_rows, datalist)


#read-in the conversion table: orthofinder sequence IDs and corresponding actual gene/protein names
#with "str_remove_all" I also easily get rid of the ":"
message("Reading in conversion table based on SequenceIDs.txt:")

conversion <- readr::read_delim(file = "orthofinder/SequenceIDs.txt",
                delim = " ",
                col_names = c("seqid", "name"),
                col_types = c(
                   seqid = col_character(),
                   name = col_character())               
               )
conversion <- conversion %>%  dplyr::mutate(seqid = str_remove_all(seqid, ":"))


message("Creating reformatted output with actual gene/protein names:")
#create new columns with the actual gene/protein names
all_BLAST_reformatted <- all_BLAST %>% dplyr::left_join(conversion, by = c("qseqid" = "seqid"))
all_BLAST_reformatted <- all_BLAST_reformatted %>% dplyr::left_join(conversion, by = c("sseqid" = "seqid"))
#position after qseqid and sseqid respectively
all_BLAST_reformatted <- all_BLAST_reformatted %>% dplyr::relocate(name.x, .after = qseqid)
all_BLAST_reformatted <- all_BLAST_reformatted %>% dplyr::relocate(name.y, .after = sseqid)
#rename to qseqid/sseqid_name
all_BLAST_reformatted <- dplyr::rename(all_BLAST_reformatted,  qseqid_name = name.x)
all_BLAST_reformatted <- dplyr::rename(all_BLAST_reformatted, sseqid_name = name.y)



message("Data read-in and reformat:")

ph_orthogroups <- readr::read_delim(file = "orthofinder/final-results/Phylogenetic_Hierarchical_Orthogroups/N0.tsv",
                          delim = "\t")

#create dataframe with numbers of genes per PHOG
#we combine expanded_in and compared_to vectors to easily compute for everything we need

HOG_df <- setNames(base::data.frame(matrix(ncol = 1 + length(all_species), nrow = 0)),
                  c("HOG", all_species))

for (i in 1:nrow(ph_orthogroups)) {
    
    row = 0
    
    #print(ph_orthogroups[i,]$HOG)
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
    HOG_df[i,] <- row
} 

HOG_tibble <- as_tibble(HOG_df)

for (k in 1:length(all_species)) {
    o = all_species[k]
    HOG_tibble[[o]] <- as.numeric(HOG_tibble[[o]])
}

#create "copy" as HOG_tibble is currently modified in the upcoming step
HOG_tibble_complete <- HOG_tibble


#here, we apply the expansion rule/s - could be also be parsed from config? (toDO!)
#get() function solves my problems of using the character functions inside dplyr
message("Applying expansion rule/s per hypothesis:")

# function has to reduce input dataset; this way each iteration (compared species) is included
# for all compared species first check for each HOG if expanded species is >= #* 
# then keep rows in which expanded species has at least # (default: 2)


for (e in expanded_in) {
    # create expansion counter column for each expanded species
    exp_species_name = paste0("exp_counter_", e)
    HOG_tibble <- HOG_tibble %>% tibble::add_column(!!(exp_species_name) := 0)

    for (c in c_t_species) {
        HOG_tibble <- HOG_tibble %>% 
                       dplyr::mutate(!!(exp_species_name) := 
                              dplyr::case_when(
                                # although cases in which the ct species has no genes are ignored via 
                                # the multiplication in the 2nd step, this underlines that we do so
                                get(c) > 0 &
                                get(e)/species_table[e, "ploidy"] >= expansion_factor*(get(c))/species_table[c, "ploidy"] ~ get(exp_species_name) + 1,
                                TRUE ~ 0,
                              )
                        ) 
        
    }
}

# then we perform summing over all exp_counter_ columns based on user choices
# at least Nmin_expanded_in expanded species that are expanded in at least Nmin_compared_to compared_to species
# create new column
# for each exp_counter_* species check if value >= y
# all passing exp_counter_* species are counted; if sum >= x retain HOG

HOG_tibble <- HOG_tibble %>%
  mutate(pass =
    # selecting columns that we want to use use for subsequent function
    select(., contains("exp_counter_")) %>%
    # use pmap on this subset to get a vector of min from each row
    # dataframe is a list so pmap works on each element of the list; here, each row  
    # we sum the occurences (per row!) of exp_counter_* cols that greater/equal to user cutoff 
    purrr::pmap_dbl(., ~ sum(c(...) >= Nmin_compared_to))
  )

# lastly, retain rows/HOGs that pass x cutoff = number of expanded species
# & drop unnecessary columns
HOG_tibble <- HOG_tibble %>%
  filter(pass >= Nmin_expanded_in) %>%
  select(-contains("exp_counter_"), -pass)


# additional optional filter1 - expanded_in gene family complete criterium
if (expanded_in_all_found == "YES") {
  HOG_tibble <- HOG_tibble %>%
  # create column with per-row sum of expanded_in species that have at least 1 gene in the HOG
  mutate(expanded_in_pass =
    select(., contains(expanded_in)) %>%
    purrr::pmap_dbl(., ~ sum(c(...) >= 1))
  ) %>%
  # only keep rows/HOGs in which at least 1 gene of each expanded_in species occurs
  filter(expanded_in_pass >= length(expanded_in)) %>%
  # remove expanded_in_pass column
  select(-expanded_in_pass)
}


# additional optional filter2 - compared_to gene family complete criterium
if (compared_to_all_found == "YES") {
  HOG_tibble <- HOG_tibble %>%
  # create column with per-row sum of compared_to species that have at least 1 gene in the HOG
  mutate(compared_to_pass =
    select(., contains(compared_to)) %>%
    purrr::pmap_dbl(., ~ sum(c(...) >= 1))
  ) %>%
  # only keep rows/HOGs in which at least 1 gene of each compared_to species occurs
  filter(compared_to_pass >= length(compared_to)) %>%
  # remove compared_to_pass column
  select(-compared_to_pass)
}

# new object: expanded_HOGs
expanded_HOGs <- HOG_tibble


# based on filtering criteria create per hypothesis table with:
# all HOGs and gene counts + additional column expansion (yes/no)
# this is later used in the creation of final tea outputs to create a HOG level table per hypothesis

#create expansion vector with length of expandsion HOGs
expansion <- replicate(nrow(expanded_HOGs), "yes")
expansion_tibble <- full_join(HOG_tibble_complete, tibble::add_column(expanded_HOGs, expansion, .after = "HOG"), 
                       by = c("HOG"), suffix = c("", ".remove")) %>% 
                           tidyr::replace_na(list(expansion = "no")) %>%
                               select(-c(ends_with(".remove"))) 

dir.create(paste("tea/", num, "/expansion_tibble/", sep = ""))
saveRDS(expansion_tibble, paste("tea/", num, "/expansion_tibble/expansion_tibble.rds", sep = ""))


# create genes column in ph_orthogroups file
# row merge? - no unite function, really handy ;D
ref_ph_orthogroups <- ph_orthogroups %>% tidyr::unite("genes", all_of(all_species), sep =", ", na.rm = TRUE, remove = TRUE)


message("Creating .txt files for all expanded OGs with reciprocal best BLAST hits of species in respective hypothesis:")

#for each gene/protein name in an interesting OG do (is there a file just with names per OG?):
#check "all_BLAST_reformatted" for all entries including these names and create new dataframe/tibble
# then, perform filtering and retain certain set of genes/proteins per OG analysis
#  then, create .txt file per OG with these gene names
## come-up with filter criteria to have better trees?
## I could of course just keep everything and save the evalues, etc.; well, problem for later.. ;D
####> output for snakemake? what about inidividual OG txt files, because starting here parallelisation can really impact
dir.create(paste("tea/", num, "/exp_OGs_proteinnames/", sep = ""))
#dir.create(paste("tea/", num, "/extended_BLAST_hits/", sep = ""))

## define custom class for extended blast hits
# need a list object to hold all data of this class
extended_BLAST_hits <- list()

# class for extended BLAST hits info
setClass("extended_BLAST_hits", 
         slots=list(blast_table="tbl_df",
                    num_genes_HOG="numeric",
                    num_genes_extend="numeric",
                    num_genes_complete="numeric",
                    genes_HOG="tbl_df",
                    genes_extend_hits="tbl_df")
         )


for (i in expanded_HOGs$HOG) {
    exp_og_genes <- unlist(strsplit(ref_ph_orthogroups[ref_ph_orthogroups$HOG == i,]$genes, split = ", "))
    BLAST_hits_exp_og_genes <- dplyr::filter(all_BLAST_reformatted, 
                                             qseqid_name %in% exp_og_genes | sseqid_name %in% exp_og_genes)
    sorted_BLAST_hits_exp_og_genes <- arrange(BLAST_hits_exp_og_genes, evalue, -bitscore, -pident)
    # add number of chosen additional best blast hits to size of HOG 
    HOG_set_extended <- length(exp_og_genes) + add_blast_hits
    
        # get gene name of last gene to be added based on number of add_blast_hits
    all_blast_genes <- na.omit(
       unique(
         c(
           rbind(
             sorted_BLAST_hits_exp_og_genes$qseqid_name, 
             sorted_BLAST_hits_exp_og_genes$sseqid_name
           )
         )
       )
     ) 
    
    # set of all extended blast hits (based on threshold) - vector of gene names (ordered!)
    # also nice: don't need a conditional since `%>% head(n = add_blast_hits)` will work,
    # even if add_blast_hits param is > setdiff(all_blast_genes, exp_og_genes) 
    extended_blast_hits_genes <- setdiff(all_blast_genes, exp_og_genes) %>% head(n = add_blast_hits)
    # length(extended_blast_hits_genes)
  
    # non redundant set of gene names of HOG + n additional blast hits as defined in the user threshold
    HOG_and_ext_blast_hits_genes <- c(exp_og_genes, extended_blast_hits_genes)

    #create subset of sorted_BLAST_hits_exp_og_genes table in which only:
    # exp_og_genes & extended_blast_hits_genes are allowed to appear
    # this way we have cutoff for the nth best blast hit/gene but also keep all secondary hits
    HOG_and_ext_blast_hits_table <- sorted_BLAST_hits_exp_og_genes %>%
                                      filter(qseqid_name %in% HOG_and_ext_blast_hits_genes) %>%
                                      filter(sseqid_name %in% HOG_and_ext_blast_hits_genes)

    write_lines(HOG_and_ext_blast_hits_genes,
                paste("tea/", num, "/exp_OGs_proteinnames/", i, ".txt", sep = "")
    )
 
    # for each exp. HOG create an extended_BLAST_hits S4 object and collect as part of list
    ext_B_hits <- new("extended_BLAST_hits",
      blast_table=HOG_and_ext_blast_hits_table,
      num_genes_HOG=length(exp_og_genes),
      num_genes_extend = length(extended_blast_hits_genes),
      num_genes_complete=length(HOG_and_ext_blast_hits_genes),
      genes_HOG=as_tibble(exp_og_genes),
      genes_extend_hits=as_tibble(extended_blast_hits_genes)
            )
    # assign name based on name of the underlying expanded HOG
    ext_B_hits <- list(ext_B_hits)
    names(ext_B_hits) <- paste0(i)
    
    # append to list object
    extended_BLAST_hits <- c(extended_BLAST_hits, ext_B_hits)
}

# save extended BLAST hits to hypothesis specific ("num") RDS file 
#-> to be read and used in final_tea_computation.R script
saveRDS(extended_BLAST_hits, paste("tea/", num, "/extended_BLAST_hits/extended_BLAST_hits.RDS", sep = ""))

#lastly create .check to know it's done
message("Creating .check - Expansions successfully computed for hypothesis ", num)
exp_OGs_proteinnames.check <- "check"
dir.create(paste("checks/tea/", num, "/", sep=""))
write_file(exp_OGs_proteinnames.check, paste("checks/tea/", num, "/exp_OGs_proteinnames.check", sep = ""))
