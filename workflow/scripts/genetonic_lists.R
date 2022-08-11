library(GeneTonic)
library(topGO)
library(dplyr)
library(stringr)
library(tidyr)
library(DESeq2)


#get dea object
dea_species_raw = readRDS(snakemake@input[["dea"]])

#get functional annotation
SFA = readRDS(snakemake@input[["SFA"]])

#get current species
species = snakemake@params[["species"]]

#get current ontology
ontology = snakemake@params[["ontology"]]


#Part1 - dds
dea_species <- dea_species_raw

rowData(dea_species) <- rowData(dea_species) %>%
  as.data.frame() %>%
  mutate(gene_id=rownames(.), SYMBOL=rownames(.)) %>%
  dplyr::select(gene_id, SYMBOL)


#Part2 - dea
dea_species_filtered <- dea_species_raw
keep <- rowSums(counts(dea_species_filtered) >= 10) >= 6
dea_species_filtered <- dea_species_filtered[keep, ]

results_dea_species <- results(dea_species_filtered)

results_dea_species$SYMBOL <- rownames(rowData(dea_species_filtered))


#Part 3 - res_enrich - need to create seperate for all ontologies

#quick solution to generate an all OGs set (no singletons! for now)
SFA_species <- SFA[[species]]

all_set <- bind_rows(SFA_species) %>%
  dplyr::select(`Protein-Accession`, `Gene-Ontology-Term`) %>%
  #we should remove all rows that have NO GO terms associated
  #and remove all NAs from the Gene-Ontology-Term column if they have at least 1 GO term
  #remove inline NAs in Gene-Ontology-Term column
  mutate(
    `Gene-Ontology-Term` = str_remove_all(`Gene-Ontology-Term`, "NA, |, NA") 
  ) %>% 
  #remove lines with only NA in Gene-Ontology-Term column
  filter(!str_detect(`Gene-Ontology-Term`, 'NA'))

universe <- all_set %>%
      mutate(universe_list = setNames(.[["Gene-Ontology-Term"]], pull(.["Protein-Accession"]))) %>%
      pull(universe_list)

#important to create true gene 2 GO mapping; GO column char vector!
universe <- universe %>%
  strsplit(split = ", ")

#reverse - all genes per GO term
genes_per_GO <- all_set %>%
      #reverse set for enrich table should only consist of genes with actual expression
      filter(`Protein-Accession` %in% rownames(results_dea_species)) %>%
      #normal formatting
      separate_rows(`Gene-Ontology-Term`, sep = ", ") %>%
      group_by(`Gene-Ontology-Term`) %>% 
      summarise(`Protein-Accession` = paste(`Protein-Accession`, collapse = ",")) %>%
      mutate(reverse_all_set = setNames(.[["Protein-Accession"]], pull(.["Gene-Ontology-Term"]))) %>%
      pull(reverse_all_set)

#whatever is compared to the universe it is always a subset of it!
int_set_df <- results_dea_species %>%
  as.data.frame() %>%
  mutate(id = rownames(.)) %>%
  filter(padj <= 0.1)

int_set <- int_set_df %>% pull(id)

inGenes <- factor(as.integer(names(universe) %in% int_set))
names(inGenes) <- names(universe)


GOdata <- new("topGOdata", 
              ontology=ontology, 
              allGenes=inGenes,
              annot=annFUN.gene2GO, 
              gene2GO=universe)

#Perform Fisher's exact test:
resultFisher <- runTest(GOdata, algorithm="weight01", statistic = "fisher")
#We can list the top significant results found:
allRes <- GenTable(GOdata, 
                   classicFisher = resultFisher, 
                   orderBy = "resultFisher", 
                   topNodes = 500, 
                   numChar=1000)

allRes <- allRes %>%
  dplyr::rename(p.value_classic = classicFisher) %>%
  mutate(p.value_classic = as.numeric(p.value_classic)) %>%
  #safety for unwanted GO terms
  filter(GO.ID %in% names(genes_per_GO)) %>%
  rowwise() %>%
  #if elim or another algorithm that adds related GO terms were used, 
  #we would need to catch these instances with case_when...
  mutate(genes = genes_per_GO[[GO.ID]])


#final step for part3 - get this into GeneTonic format
res_enrich_allRes <- shake_topGOtableResult(obj = allRes, p_value_column = "p.value_classic")


#Part 4
anno_df <- SFA_species %>%
  dplyr::select(`Protein-Accession`) %>%
  dplyr::rename(gene_id=`Protein-Accession`) %>%
  mutate(gene_name = gene_id) %>%
  as.data.frame()

rownames(anno_df) <- anno_df$gene_id


#enrich scores
res_enrich_scores <- get_aggrscores(res_enrich = res_enrich_allRes,
                                    res_de = results_dea_species,
                                    annotation_obj = anno_df,
                                    aggrfun = mean)

#remove cases of single genes for a GO term
need_removed <- c()

for (i in seq_len(nrow(res_enrich_scores))) {

  genes <- unlist(strsplit(res_enrich_scores[i, "gs_genes"], ","))
  if (length(genes) < 2){
    need_removed <- c(need_removed, i)
  }
}

res_enrich_scores_filtered <- res_enrich_scores[-need_removed,]


gtl <- GeneTonic_list(
  dds = DESeq2::estimateSizeFactors(dea_species),
  res_de = results_dea_species,
  res_enrich = res_enrich_scores_filtered,
  annotation_obj = anno_df
)

#save the object
saveRDS(object=gtl, file=snakemake@output[[1]])
