#to reduce installation time; we could reduce down to readr, plyr, dplyr, stringr
library(readr)
library(dplyr)
library(stringr)


message("Acquiring hypothesis variables:")
num = snakemake@params[["num"]]
name = snakemake@params[["name"]]
expanded_in = snakemake@params[["expansion"]]
compared_to = snakemake@params[["comparison"]]

#read-in Orthogroup-GeneCounts-table
message("Reading in Orthogroup-GeneCounts-table:")
OG.GC <- readr::read_tsv("orthofinder/final-results/Orthogroups/Orthogroups.GeneCount.tsv")


#here, we apply the expansion rule/s - could be also be parsed from config? (toDO!)
#get() function solves my problems of using the character functions inside dplyr
message("Applying expansion rule/s per hypothesis:")
expanded_OGs <- OG.GC
for (i in compared_to) {
    expanded_OGs <- expanded_OGs %>% dplyr::filter(get(expanded_in) > 3*(get(i)))
    }


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


#read-in Orthogroups.txt
message("Reading in Orthogroups.txt:")
orthogroups <- readr::read_delim(file = "orthofinder/final-results/Orthogroups/Orthogroups.txt",
                          delim = ":",
                          col_names = c("OG", "genes"),
                          col_types = c(
                              OG = col_character(),
                              genes = col_character())
                         )


#for each gene/protein name in an interesting OG do:
#check "all_BLAST_reformatted" for all entries including these names and create new dataframe/tibble
# then, perform filtering and retain certain set of genes/proteins per OG analysis
#  then, create .txt file per OG with these gene names
## come-up with filter criteria to have better trees?
## I could of course just keep everything and save the evalues, etc.; well, problem for later.. ;D
####> we create the actual output files here; incl. a genereal ".check" - the sub-/directories in the path are created in the Snakefile

#additional quotation mark as first element is introduced for some reason
#(delim space after : isn't handled properly)
#workaround, we get rid with trailing [-1]

message("Creating .txt files for all expanded OGs with reciprocal best BLAST hits of species in respective hypothesis:")

for (i in expanded_OGs$Orthogroup) {
    exp_og_genes <- unlist(strsplit(orthogroups[orthogroups$OG == i,]$genes, split = " "))[-1]
    BLAST_hits_exp_og_genes <- dplyr::filter(all_BLAST_reformatted, 
                                             qseqid_name %in% exp_og_genes | sseqid_name %in% exp_og_genes)
    sorted_BLAST_hits_exp_og_genes <- arrange(BLAST_hits_exp_og_genes, evalue, -bitscore, -pident)
    list_qseqid <- as.character(sorted_BLAST_hits_exp_og_genes$qseqid_name)
    list_sseqid <- as.character(sorted_BLAST_hits_exp_og_genes$sseqid_name)
    list_merged <- unique(c(list_qseqid, list_sseqid))
    write_lines(list_merged,
           paste("tea/", num, "/exp_OGs_proteinnames/proteinnames_", i, ".txt", sep = ""))
}

#lastly create .check to know it's done
#message("Expansion successfully computed - creating .check:")
#exp_OGs_proteinnames.check <- "check"
#write_file(exp_OGs_proteinnames.check, paste("checks/tea/", num, "/exp_OGs_proteinnames.check", sep = ""))
