import os
import pandas as pd
import csv
import posixpath
import gzip
import shutil

from urllib.request import urlretrieve
from urllib.request import urlcleanup


#create a subset of the species - those that should be annotated using AHRD
species_df = pd.read_csv('config/species.tsv', sep='\t')
ahrd_species_df = species_df.query("function=='AHRD'")
non_ahrd_species_df = species_df.query("function!='AHRD'")
AHRD_SPECIES = ahrd_species_df["species"].tolist()
NON_AHRD_SPECIES= non_ahrd_species_df["species"].tolist()

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


#function handle supplied function annotation tables
#checking if URL is supplied and if gzipped based on functions currently part of ref_utils.smk
def handle_user_func_annotation_table(file, species):
    urlcleanup()
    if is_url(file) and not is_gzipped(file):
        print("Downloading user supplied functional annotation:")
        urlretrieve(annotation, 'resources/functional_annotation/' + species + '.func_annotation.tsv')
    if is_url(file) and is_gzipped(file):
        print("Downloading user supplied functional annotation (gzipped):")
        urlretrieve(annotation, 'resources/functional_annotation/' + species + '.func_annotation.tsv.gz')
    if not is_url(file) and not is_gzipped(file):
        print("Linking annotation")
        os.symlink(os.path.abspath(file), 'resources/functional_annotation/' + species + '.func_annotation.tsv')
    if not is_url(file) and is_gzipped(file):
        print("Linking annotation")
        os.symlink(os.path.abspath(file), 'resources/functional_annotation/' + species + '.func_annotation.tsv.gz')
    if is_gzipped(file):
        gunzip_annotation('resources/functional_annotation/' + species + '.func_annotation.tsv.gz')


#functions to check validity of user supplied function annotation table
#we require a tab seperated table with at least one column "Protein-Accession" & another column "Gene-Ontology-Term"
#  contained gene/protein identifier & corresponding gene ontology terms (seperated by ", ") respectively
def check_user_func_annotation_table_validity(file, species):
    check_tsv = pd.read_csv(os.path.abspath(file),  sep='\t')
    #create list of column headers
    check_tsv_columns = list(check_tsv.columns)
    #create list of the two mandatory column headers
    mandatory_list = ["Protein-Accession", "Gene-Ontology-Term"]
    #check whether or not all mandatory headers are found in the user function annotation table
    if not set(mandatory_list).issubset(check_tsv_columns):
        print("Necessary headers Protein-Accession & Gene-Ontology-Term not found in function annotation table for species:", species, "!!!")
        print("Check the file - perhaps it is not tab seperated?")
        exit()
    else:
        print("User supplied functional annotation file for species - ", species, "- checked")

