# we are using readr 1.4.0 from conda - doesn't support some of the cool newer features like col_select
library(readr)
library(stringr)
library(dplyr)


ahrd_files <- snakemake@params[["ahrd_files"]]
renamed_user_files <- snakemake@input[["renamed_user_files"]]

#compute length of list - could also get this via Snakemake and the N_SPECIES constant
n_species <- length(ahrd_files) + length(renamed_user_files)

#"SFA" stands for Species Functional Annotation 
SFA <- vector(mode = "list", length = n_species)


###we start woth the AHRD annotated species
#here we just add to the list starting at the the index positon 1

#To Do - make this nice xD
#I know there are better ways to do this in R - #never use for loops 
for (i in 1:length(ahrd_files)) {

  #drop paths and assign species name to "species_name"
  species_name <- str_replace(ahrd_files[i], "results/", "")
  species_name <- str_replace(species_name, ".ahrd_output.tsv", "")
  #assign tibble to species_name - skipping two rows since current AHRD output comes with two useless header rows
  assign(species_name, read_tsv(ahrd_files[i], skip=2, col_names=TRUE))
  #drop AHRD-Quality-Code (immediately with col_select - only suppported by readr >2.0.0 - conda currently 1.4.0)
  assign(species_name, get(species_name) %>% select(-"AHRD-Quality-Code"))

  #assign name of species to the ith element in list
  names(SFA)[[i]] <- species_name
  #assign actual tibble with functional annotation of species to ith element of list
  SFA[[i]] <- get(species_name)

}


###in a second step we deal with user supplied species functional annotation files (prev. validity checked with python)
#here again: To Do change this from for loop to sth. nicer
for (i in 1:length(renamed_user_files)) {

  #drop paths and assign species name to "species_name"
  species_name <- str_replace(renamed_user_files[i], "resources/functional_annotation/", "")
  species_name <- str_replace(species_name, ".func_annotation.tsv", "")

  #assign tibble to species_name
  assign(species_name, read_tsv(renamed_user_files[i], col_names=TRUE))

  #assign name of species to the length(ahrd_files)+ith element in list
  names(SFA)[[length(ahrd_files)+i]] <- species_name
  #assign actual tibble with functional annotation of species to the length(ahrd_files)+ith element of list
  SFA[[length(ahrd_files)+i]] <- get(species_name)

}

#save SFA
saveRDS(SFA, "results/functional_annotation/species_functional_annotation.rds")
