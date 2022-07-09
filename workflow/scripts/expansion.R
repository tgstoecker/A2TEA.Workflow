#backup R-based installation if conda didn't work or wasn't used
#we check if packages are installed first

# list of cran packages
cran_packages = c("readr", "plyr", "dplyr", "stringr", "tidyr", "tibble", "reshape2", "foreach", "doParallel")
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
library(plyr)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(reshape2)
library(foreach)
library(doParallel)

message("Acquiring hypothesis variables:")
num = snakemake@params[["num"]]
name = snakemake@params[["name"]]

#threads info for parallel FORK cluster
threads = as.numeric(snakemake@threads)

# define the number of additional best blast hits to include in the follow-up analyses
add_OGs = snakemake@params[["add_OGs"]]

# define the minimum expansion factor & expansion difference to call an expansion (part of hypothesis.tsv)
# + minimum number of genes of expanded species to even consider an OG
expansion_factor = as.numeric(snakemake@params[["expansion_factor"]])
expansion_difference = as.numeric(snakemake@params[["expansion_difference"]])
expanded_genes_min = as.numeric(snakemake@params[["expanded_genes_min"]])

# get user choice whether to perform ploidy normalization or not
ploidy_normalization = as.character(snakemake@params[["ploidy_normalization"]])

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

    if (ploidy_normalization == "YES") {
        for (c in c_t_species) {
            HOG_tibble <- HOG_tibble %>% 
                            mutate(!!(exp_species_name) := 
                            dplyr::case_when(
                                # although cases in which the ct species has no genes are ignored via 
                                # the multiplication in the 2nd step, this underlines that we do so
                                # now with the second expansion distance criterium we should have it anyway
                                get(c) > 0 &
                                get(e)/species_table[e, "ploidy"]*2 >= expansion_factor*(get(c))/species_table[c, "ploidy"]*2 ~ get(exp_species_name) + 1,
                                get(c) > 0 &
                                get(e)/species_table[e, "ploidy"]*2 >= expansion_difference + (get(c))/species_table[c, "ploidy"]*2 ~ get(exp_species_name) + 1,
                                TRUE ~ 0,
                              )
                        )  
        }
    }
    #if the user did not set ploidy_normalization to "YES"
    #here, no ploidy information is used to divide the number of genes per OG
    else {
        for (c in c_t_species) {
            HOG_tibble <- HOG_tibble %>% 
                            mutate(!!(exp_species_name) := 
                            dplyr::case_when(
                                get(c) > 0 &
                                get(e) >= expansion_factor*(get(c)) ~ get(exp_species_name) + 1,
                                get(c) > 0 &
                                get(e) >= expansion_difference + get(c) ~ get(exp_species_name) + 1,
                                TRUE ~ 0,
                              )
                        )
        }
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

# optional hard filter for at least # genes in all expanded species
# this is useful in difficult ploidy cases and solves downstream issues in small OGs
HOG_tibble <- HOG_tibble %>%
  filter(if_all(contains(expanded_in), ~ . > expanded_genes_min))

# new object: expanded_HOGs
expanded_HOGs <- HOG_tibble


# based on filtering criteria create per hypothesis table with:
# all HOGs and gene counts + additional column expansion (yes/no)
# this is later used in the creation of final tea outputs to create a HOG level table per hypothesis

#need to check if there is at least one expanded (H)OG
#if not replace_na woud throw error since we are changing types
#all of the downstream stuff in this script only makes sense if we in fact do have expanded groups
if (nrow(expanded_HOGs) > 0) {

  #create expansion vector with length of expandsion HOGs
  expansion <- replicate(nrow(expanded_HOGs), "yes")
  #print(expansion)
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
  dir.create(paste("tea/", num, "/expansion_cp_target_OGs/", sep = ""))

  ## define custom class for extended blast hits
  # need a list object to hold all data of this class
  extended_BLAST_hits <- list()

  # class for extended BLAST hits info
  setClass("extended_BLAST_hits", 
           slots=list(blast_table="tbl_df")
          )


  for (i in expanded_HOGs$HOG) {
    exp_og_genes <- unlist(strsplit(ref_ph_orthogroups[ref_ph_orthogroups$HOG == i,]$genes, split = ", "))
    BLAST_hits_exp_og_genes <- dplyr::filter(all_BLAST_reformatted, 
                                             qseqid_name %in% exp_og_genes | sseqid_name %in% exp_og_genes)
    sorted_BLAST_hits_exp_og_genes <- arrange(BLAST_hits_exp_og_genes, evalue, -bitscore, -pident)
    
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
    extended_blast_hits_genes <- setdiff(all_blast_genes, exp_og_genes) %>% head(n = add_OGs)
  
    # non redundant set of gene names of HOG + n additional blast hits as defined in the user threshold
    HOG_and_ext_blast_hits_genes <- c(exp_og_genes, extended_blast_hits_genes)

    #create subset of sorted_BLAST_hits_exp_og_genes table in which only:
    # exp_og_genes & extended_blast_hits_genes are allowed to appear
    # this way we have cutoff for the nth best blast hit/gene but also keep all secondary hits
    HOG_and_ext_blast_hits_table <- sorted_BLAST_hits_exp_og_genes %>%
                                      filter(qseqid_name %in% HOG_and_ext_blast_hits_genes) %>%
                                      filter(sseqid_name %in% HOG_and_ext_blast_hits_genes)

#tt    write_lines("",
#tt                paste("tea/", num, "/expansion_cp_target_OGs/", i, ".txt", sep = "")
#tt    )
 
    # for each exp. HOG create an extended_BLAST_hits S4 object and collect as part of list
    ext_B_hits <- new("extended_BLAST_hits",
      blast_table=HOG_and_ext_blast_hits_table
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


  ### Adding OGs instead of BLAST hits ###
  message("Adding OGs instead of BLAST hits")

  #transforming ph_orthogroups to long format - nice & neat lookup table ;D
  long_ph_orthogroups <- ph_orthogroups %>%
    select(-OG, -`Gene Tree Parent Clade`) %>%
    melt(., id.vars = c("HOG")) %>%
    rename(species=variable, id=value) %>%
    mutate(species = as.character(species)) %>%
    separate_rows(id, sep = ", ") %>%
    drop_na()

  #create final summary list - per OG (name) the cumulative steps of additional OGs
  summary_add_OG_analysis_list <- vector(mode = "list", length = length(expanded_HOGs$HOG))

  length(summary_add_OG_analysis_list)

  #create classes for nested S3-list structure holding all additional OG sets per OG analysis
  setClass("add_OG_analysis", 
           slots=list(add_OG_analysis="list")
          )

  setClass("add_OG_set", 
           slots=list(genes="tbl_df")
          )


  dir.create(paste("tea/", num, "/add_OGs_sets/id_lists/", sep = ""))

  #removing all big files to minimize mem impact of FORK cluster
  rm(datalist)
  rm(all_BLAST)
  rm(all_blast_genes)
  rm(HOG_and_ext_blast_hits_table)
  rm(HOG_and_ext_blast_hits_genes)
  rm(sorted_BLAST_hits_exp_og_genes)
  rm(extended_BLAST_hits)
  rm(ext_B_hits)

  #### FORK cluster since I expect a Linux machine
  #### autostop=TRUE since I don't want to handle this manually
  #with my.cluster & stopCluster(my.cluster) I could check the status

  setup_cluster <- function(){

    #define cluster
    parallel::detectCores()
    n.cores <- threads
    n.cores

    #create the cluster - FORK because this way libraries, variables etc. are copied to the clusters!
    my.cluster <- parallel::makeForkCluster(
      n.cores, 
      type = "FORK",
      autostop=TRUE
    )
   
    #check cluster definition (optional)
    print(my.cluster)

    #register it to be used by %dopar%
    doParallel::registerDoParallel(cl = my.cluster)

    #check if it is registered (optional)
    print(
      foreach::getDoParRegistered()
    )
    #how many workers are available? (optional)
    print(
      foreach::getDoParWorkers()
    )

  }

  #function to completely remove a fork cluster
  burn_socks <- function(x){
    close.connection(getConnection(x))
  }

  #function to truly get rid of old Cluster/sockets
  rm_cluster <- function(){
    stopImplicitCluster()

    connections <- showConnections(all = FALSE) 
    socket_connections <- as.data.frame(connections) %>%
      filter(class == "sockconn") %>%
      rownames()

    message("Removing all unwanted FORK connections - purging closed cluster sockets")
    message("This will kill zombie proesses and free up RAM")

    lapply(X = socket_connections, 
           FUN = burn_socks)

  }

  #setup & start FORK cluster
  setup_cluster()

  #we iterate over the expanded OGs
  pre_summary_add_OG_analysis_list <- foreach(i = expanded_HOGs$HOG) %dopar% {
    
    exp_og_genes <- unlist(strsplit(ref_ph_orthogroups[ref_ph_orthogroups$HOG == i,]$genes, split = ", "))
    BLAST_hits_exp_og_genes <- dplyr::filter(all_BLAST_reformatted, 
                                             qseqid_name %in% exp_og_genes | sseqid_name %in% exp_og_genes)
    sorted_BLAST_hits_exp_og_genes <- arrange(BLAST_hits_exp_og_genes, evalue, -bitscore, -pident)

    #after the sorting of BLAST hits we move on to merging the OG information into the current exp. OG BLAST hits results

    #we can create a merged dataframe from the blast table (HOG specific) and the long format HOG table (general)
    #difficulty is that in case of singletons we have NA under HOG
    #easiest way I can think of is a named list
    self_and_closest_ogs <- left_join(sorted_BLAST_hits_exp_og_genes, long_ph_orthogroups, by = c("sseqid_name" = "id")) %>%
      group_by(sseqid_name) %>%
      arrange(evalue, -bitscore, -pident) %>%
      slice(1) %>%
      ungroup() %>%
      arrange(evalue, -bitscore, -pident) %>%
      mutate(HOG = as.character(HOG),
             HOG = ifelse(is.na(HOG), paste("singleton", sep = "-", cumsum(is.na(HOG))),
                          as.character(HOG))) %>%
      group_by(HOG) %>%
      arrange(evalue, -bitscore, -pident) %>%
      slice(1) %>%
      ungroup() %>%
      arrange(evalue, -bitscore, -pident)
     
    #first row will should correspond to self match - the OG itself
    #self_and_closest_ogs

    #here the information provided by the user regarding max. additional OGs is used
    #we copy the value to an "internal" object since in the following commands we modify in each iteration depending on the underlying OG/s
    #the initial add_OGs value chosen by the user is used later and must not be modified
    add_OGs_internal = add_OGs

    #need to add check for number of additonal OGs/singletons
    #e.g. user sets max. add. OGs/singletons to 5 but we only can provide 4
    #+ have to consider the expanded OG itself, so:
    available_add_OGs <- nrow(self_and_closest_ogs) - 1 
    #available_add_OGs

    #in this case set add_OGs to max available
    if (available_add_OGs < add_OGs_internal) {
      add_OGs_internal <- available_add_OGs
    } 

    #empty list of precomputed size - add OGs plus the expanded OG itself
    add_OG_analysis_list <- vector(mode = "list", length = add_OGs_internal + 1)


    for (j in 1:(add_OGs_internal + 1)) {
      og_name <- self_and_closest_ogs[j,] %>%
      pull(HOG)
  
      #differnetiate between OG case and a singleton case - different handling:
      #for HOG get all associated genes/proteins from large long format table
      #for singletons get singelton gene/protein from the table itself
      if (!str_detect(og_name, "singleton") == TRUE) {
  
        name_curr_close_HOG <- self_and_closest_ogs[j,] %>%
                                 pull(HOG)
  
        ids_curr_close_HOG <- long_ph_orthogroups %>%
                                filter(HOG %in% name_curr_close_HOG) %>%
                                pull(id)
  
        add_OG_analysis_list[[j]] <- ids_curr_close_HOG
        names(add_OG_analysis_list)[j] <- name_curr_close_HOG
  
      } else {
        name_curr_close_singleton <- self_and_closest_ogs[j,] %>%
                                     pull(HOG)
  
        id_curr_close_HOG <- self_and_closest_ogs[j,] %>%
                               pull(sseqid_name)
  
        add_OG_analysis_list[[j]] <- id_curr_close_HOG
        names(add_OG_analysis_list)[j] <- name_curr_close_singleton
      }
    }
    
    #create copy of the list and remove the names
    #add all previous gene/proteins to the next OG
    # -> e.g. the former seventh OG list element will contain all genes/proteins of all lower numbered OGs
    #the final "cum_add_OG_analysis_list" will be a list with each next element having the cumulative set of all previous and own gene/protein ids 
    cum_add_OG_analysis_list <- add_OG_analysis_list
    names(cum_add_OG_analysis_list) <- NULL
    #copy_add_OG_analysis_list

    for (k in 0:(length(cum_add_OG_analysis_list)-1)) {
      cum_add_OG_analysis_list[[k+1]] <- unlist(c(cum_add_OG_analysis_list[k], cum_add_OG_analysis_list[k+1]))
    }

    #create directory/ies - for each HOG under "add_OGs_sets/id_lists/"
    dir.create(paste("tea/", num, "/add_OGs_sets/id_lists/", i, sep = ""))
    
    #iterate over cum. list and write cumulative sets into seperate files
    for (l in 1:length(cum_add_OG_analysis_list)) {
      write_lines(cum_add_OG_analysis_list[[l]],
                  paste("tea/", num, "/add_OGs_sets/id_lists/", i, "/add_OGs_set_num-", l, ".txt", sep = "")
                  #paste("../tea/", num, "/add_OGs_sets/", i, "/add_OGs_set_num-", l, ".txt", sep = "")
      )

    }

    #append current cum. list to summary list & name element after OG
    #actually we need to do this again and this is redundant;toDo
    names(cum_add_OG_analysis_list) <- i
    return(cum_add_OG_analysis_list)
    
  }

  #after cluster is done - remove it
  rm_cluster()

  #some final formatting
  for (i in 1:length(pre_summary_add_OG_analysis_list)) {
    names(pre_summary_add_OG_analysis_list[[i]]) <- NULL
  }

  names(pre_summary_add_OG_analysis_list) <- expanded_HOGs$HOG

  summary_add_OG_analysis_list <- pre_summary_add_OG_analysis_list



  ## in this final step we use the summary list to create S3-list nest structure which contains all info 
  ## of the summary but in the format we will embed in the final results object as part of the final_tea computation

  #length of add OGs analysis summary list
  summary_length <- length(summary_add_OG_analysis_list)

  #create empty list of length of the summary list
  add_og_complete_object <- vector(mode = "list", length = summary_length)

  #iterate over all OGs in summary list
  for (og in 1:summary_length) {
 
    #amount of sets for current OG - (necessary since less addtional OGs than max wished by user is possible)
    curr_og_n_sets <- length(summary_add_OG_analysis_list[[og]])
    
    #empty list with n elements = n sets for currently analyzed OG
    og_all_sets <- vector(mode = "list", length = curr_og_n_sets)
    
    for (set in 1:curr_og_n_sets) {

      curr_set <- new("add_OG_set",
                      genes=as_tibble(
                        unlist(
                          summary_add_OG_analysis_list[[og]][set],
                        )  
                      ) 
                  )
      curr_set <- list(curr_set)
  
      og_all_sets[set] <- curr_set
      names(og_all_sets)[set] <- paste0("set_", set)
      
    }
    
    curr_add_og <- new("add_OG_analysis",
                       add_OG_analysis=og_all_sets
                      )

    curr_add_og <- list(curr_add_og) 
 
    add_og_complete_object[og] <- curr_add_og

    names(add_og_complete_object)[og]  <- names(summary_add_OG_analysis_list[og])
  }

  # save summary table for aditional OG analysis to hypothesis specific ("num") RDS file
  saveRDS(add_og_complete_object, paste("tea/", num, "/add_OGs_object/add_OG_analysis_object.RDS", sep = ""))

  # nested snakemake checkpoints are annoying at the moment 
  # quick fix - create empty addtional set_num files which we can ignore but nonetheless exist
  message("adding empty missing files - current workaround to avoid nested snakemake checkpoints")

  #max_n_sets + 1 since user choice is the number of ADDITIONAL sets but we have ti remember the expanded OG itself (always set_num = 1)
  max_n_sets = add_OGs + 1

  for (n in names(add_og_complete_object)) {
    
    og_n_sets <- length(
         list.files(path = paste0("tea/", num, "/add_OGs_sets/id_lists/", n))
        )
    
    og_sets <- list.files(path = paste0("tea/", num, "/add_OGs_sets/id_lists/", n))
    
    if (og_n_sets < max_n_sets) {
        n_missing_sets <- max_n_sets - og_n_sets
        
        for (m in 1:n_missing_sets) {
          value <- max_n_sets - m + 1
          missing_file <- paste0("tea/", num, "/add_OGs_sets/id_lists/", n, "/", "add_OGs_set_num-", value, ".txt")
            
          file.create(missing_file)
        }
    }
    
  }


#all of the previous only runs in case we actually find an expanded group in this particular hypothesis
#if this not the case we have to perform considerably less work (although of course the hypothesis is rather uninformative)
} else {
  #just add expansion column with "no" for all rows
  expansion_tibble <- HOG_tibble_complete %>%
                        mutate(expansion = "no")

  dir.create(paste("tea/", num, "/expansion_tibble/", sep = ""))
  saveRDS(expansion_tibble, paste("tea/", num, "/expansion_tibble/expansion_tibble.rds", sep = ""))

  message("No expanded OGs for this hypothesis - creating empty outputs")
  extended_BLAST_hits <- "empty"
  # save extended BLAST hits to hypothesis specific ("num") RDS file 
  #-> to be read and used in final_tea_computation.R script
  saveRDS(extended_BLAST_hits, paste("tea/", num, "/extended_BLAST_hits/extended_BLAST_hits.RDS", sep = ""))

  add_og_complete_object <- "empty"
  # save summary table for aditional OG analysis to hypothesis specific ("num") RDS file
  saveRDS(add_og_complete_object, paste("tea/", num, "/add_OGs_object/add_OG_analysis_object.RDS", sep = ""))
 
  dir.create(paste("tea/", num, "/expansion_cp_target_OGs/", sep = ""))

  #create empty file under add_OG id_lists
  dir.create(paste0("tea/", num, "/add_OGs_sets/id_lists/"))
  file.create(paste0("tea/", num, "/add_OGs_sets/id_lists/", "empty.txt"))

}

#####

#### Lastly, create .check to know everything is done
#message("Creating .check - Expansions successfully computed for hypothesis ", num)
#expansion_cp_target_OGs.check <- "check"
#dir.create(paste("checks/tea/", num, "/", sep=""))
#write_file(expansion_cp_target_OGs.check, paste("checks/tea/", num, "/expansion_cp_target_OGs.check", sep = ""))
