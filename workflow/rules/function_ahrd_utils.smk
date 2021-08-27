import os
import pandas
import csv


#create a subset of the species - those that should be annotated using AHRD
species_df = pd.read_csv('config/species.tsv', sep='\t')
ahrd_species_df = species_df.query("function=='AHRD'")
AHRD_SPECIES = ahrd_species_df["species"].tolist()


#function to write a custom species.tsv file for running the AHRD snakemake wrapper
def create_ahrd_species_tsv(ahrd_species):
    with open('AHRD_Snakemake/resources/species.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['species', 'pep_fasta'])
        
        for s in ahrd_species:
            s_path = os.path.abspath(s)
            s_name = os.path.basename(s)
            s_name_noext = os.path.splitext(s_name)[0]
            tsv_writer.writerow([s_name_noext, s_path])
