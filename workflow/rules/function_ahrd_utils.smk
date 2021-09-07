import os
import pandas as pd
import csv
import posixpath
import gzip
import shutil

from urllib.request import urlretrieve
from urllib.request import urlcleanup


##PART 1 - dealing with running AHRD for species the user did not supply their own functional annotation file

#create a subset of the species - those that should be annotated using AHRD
species_df = pd.read_csv('config/species.tsv', sep='\t')
ahrd_species_df = species_df.query("function=='AHRD'")
non_ahrd_species_df = species_df.query("function!='AHRD'")
AHRD_SPECIES = ahrd_species_df["species"].tolist()
NON_AHRD_SPECIES= non_ahrd_species_df["species"].tolist()

#function to write a custom species.tsv file for running the AHRD snakemake wrapper
def create_ahrd_species_tsv(ahrd_species):
    with open('workflow/rules/AHRD_Snakemake/resources/species.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['species', 'pep_fasta'])
        
        for s in ahrd_species:
            s_path = os.path.abspath(s)
            s_name = os.path.basename(s)
            s_name_noext = os.path.splitext(s_name)[0]
            tsv_writer.writerow([s_name_noext, s_path])


#read-in of AHRD species.tsv in function that is called by the next two functions (modified from AHRD_snakemake common.smk)
#the trick is that since the AHRD_Snakemake example species.tsv exists at the first evaluation no error is thrown;
#then after we have created the AHRD species.tsv based on the A2TEA species.tsv these functions are called by the AHRD_wrapper
#the new/actual species.tsv (updated in AHRD_Snakemake/resources/) is being parsed and used
def get_ahrd_species_table():
    species_table = pd.read_table("workflow/rules/AHRD_Snakemake/resources/species.tsv").set_index("species")
    return species_table

def get_species_fasta(wildcards):
    species_table = get_ahrd_species_table()
    return species_table.loc[wildcards.species]["pep_fasta"]

def calc_mem_usage_in_mb(wildcards):
    species_table = get_ahrd_species_table()
    fastaPath = species_table.loc[wildcards.species]["pep_fasta"]
    fastaSizeInByte = Path(fastaPath).stat().st_size
    fastaSizeInKB = fastaSizeInByte/1024
    fastaSizeInMB = fastaSizeInKB/1024
    projectedMemUsage = 11000+fastaSizeInMB*540
    return int(projectedMemUsage)

def species_tsv_checkpoint_end(wildcards):
    checkpoint_output = checkpoints.create_ahrd_species_tsv.get(**wildcards).output
    ahrd_results = expand("results/{species}.ahrd_output.tsv", species=AHRD_SPECIES)
    return ahrd_results


##PART 2 - dealing with user supplied functional annotation input

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
