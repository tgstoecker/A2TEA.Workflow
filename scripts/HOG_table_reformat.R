#to reduce loading time; we could reduce down to readr, plyr, dplyr, stringr
library(tidyverse)

# using snakemake propagation + python strsplit() is really cool since the input is just a vector
## even if we have multiple species in expanded, compared or both ;D
all_species <- snakemake@params[["all_species"]]
hypothesis_num <- snakemake@params[["hypothesis_num"]]

message("Started work on hypothesis:")
print(hypothesis_num)
message("With species:")
print(all_species)

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

HOG_tibble <- add_column(HOG_tibble, Desc = "(null)" , .before = "HOG")

# create "copy" for the filtered set
# also here a much better tidyverse solution likely - problem for later
HOG_tibble_filtered <- HOG_tibble

for (species in all_species) {
  HOG_tibble_filtered <- HOG_tibble_filtered %>%
    filter(get(species) < 100)
}


message("Writing filtered set (<100 genes per species per HOG) to file - HOG_table_reformatted_filtered.tsv")
write_tsv(HOG_tibble_filtered, paste0("cafe/", hypothesis_num, "/HOG_table_reformatted_filtered.tsv"))

message("Writing complete set to file - HOG_table_reformatted_complete.tsv")
write_tsv(HOG_tibble, paste0("cafe/", hypothesis_num, "/HOG_table_reformatted_complete.tsv"))
