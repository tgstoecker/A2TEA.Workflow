#Expansion Analysis
# can handle both multi expansion and multi comparison species, now ;D
# the strplit command works really well in this context
## this way, multiple species are propagated into R as a simple vector  

def get_hypo_num(wildcards):
#    """Get compared_to entries from hypotheses(.tsv) for each hypothesis. """
#    num = hypotheses.loc[ (wildcards.hypothesis), 'hypothesis']
    num = str(wildcards.hypothesis)
    return num

def get_hypo_name(wildcards):
    """Get compared_to entries from hypotheses(.tsv) for each hypothesis. """
    name = hypotheses.loc[ 'name', (wildcards.hypothesis) ]
    return name

# for both the species we want to check for expansion as well as those being compared to:
# if more than one species is compared to than we have to split the string based on ";" 

def get_exp_species(wildcards):
    """Get compared_to entries from hypotheses(.tsv) for each hypothesis. """
    exps_1 = hypotheses.loc[ 'expanded_in', (wildcards.hypothesis) ]
    if exps_1.count(";") > 0:
        exps_2 = str.split(exps_1, ";")
        return exps_2
    else:
        return exps_1
    return


#if more than one species is compared to than we have to split the string based on ";"
def get_com_species(wildcards):
    """Get compared_to entries from hypotheses(.tsv) for each hypothesis. """
    ct_1 = hypotheses.loc[ 'compared_to', (wildcards.hypothesis) ]
    if ct_1.count(";") > 0:
        ct_2 = str.split(ct_1, ";")
        return ct_2
    else:
        return ct_1
    return


#if more than one species is compared to than we have to split the string based on ";"
#the goal in any scenarios to return a list (!) with just the species used in the hypothesis
#we also need to add the path to the longest isoform peptide fasta files and add the .fa suffix
#another trick is that we can simply use the species base name in any case,
#since the longest isoform output is always named after it!
# this also applies if the user did the isoform filtering on their own since the renamed pep fastas are linked to isoform directory in that case ;D
#+the filtering for longest isoform also only retains the base name of the gene/protein
def get_all_hypothesis_species(wildcards):
    """Get compared_to entries from hypotheses(.tsv) for each hypothesis. """
    path_prefix = 'resources/longest_isoforms/'
    suffix = '.fa'
    exp = hypotheses.loc[ 'expanded_in', (wildcards.hypothesis) ]
    ct = hypotheses.loc[ 'compared_to', (wildcards.hypothesis) ]
    # split by ";", if no ";" then transform string to single-element list (so concatenation works)
    if exp.count(";") > 0:
        exp = str.split(exp, ";")
    else:
        exp = [exp]
    if ct.count(";") > 0:
        ct = str.split(ct, ";")
    else:
        ct = [ct]
    # concatenate both lists
    output = exp + ct
    # removing dups - (complex hypotheses?)
    output = list( dict.fromkeys(output) )
    # add .fa suffix
    output = [path_prefix + x + suffix for x in output]
    return output


# get user defined expansion_factor for each hypothesis by parsing hypotheses.tsv
def get_expansion_factor(wildcards):
    return hypotheses.loc["min_expansion_factor"][wildcards.hypothesis]


# get user defined Nmin_expanded_in & Nmin_compared_to for each hypothesis by parsing hypotheses.tsv
# explained: at least Nmin_expanded_in expanded species that are expanded in at least Nmin_compared_to compared_to species
def get_Nmin_expanded_in(wildcards):
    return hypotheses.loc["Nmin_expanded_in"][wildcards.hypothesis]

def get_Nmin_compared_to(wildcards):
    return hypotheses.loc["Nmin_compared_to"][wildcards.hypothesis]


# get user defined expanded_in_all_found & compared_to_all_found for each hypothesis by parsing hypotheses.tsv
# explained: at least Nmin_expanded_in expanded species that are expanded in at least Nmin_compared_to compared_to species
def get_expanded_in_all_found(wildcards):
    return hypotheses.loc["expanded_in_all_found"][wildcards.hypothesis]

def get_compared_to_all_found(wildcards):
    return hypotheses.loc["compared_to_all_found"][wildcards.hypothesis]


#SOLVING expansion checkpoint HERE!
#after rule FastTree in expansion_analysis.smk
def solve_expansion(wildcards):
    checkpoint_output = checkpoints.expansion.get(**wildcards).output[0]
    file_names = expand("tea/{hypothesis}/trees/{OG}.tree", hypothesis=wildcards.hypothesis, OG=glob_wildcards(os.path.join(checkpoint_output, "{OG}.txt")).OG)
    return file_names

