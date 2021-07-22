#neat placeholder variable for the PATH
ORTHOFINDER = "orthofinder/"

# added choice of chunk usage or not ;D
with open('config/config.yaml', 'r') as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

CHUNKS_ORTHO = config["chunks_number"]

def get_species_fasta(wildcards):
    return species_table.loc[wildcards.species]["pep_fasta"]

def get_longest_isoforms(wildcards):
    return os.path.join("FS/longest_isoforms/", os.path.split(get_species_fasta(wildcards))[1])
