#Expansion Analysis
# can handle both multi expansion and multi comparison species, now ;D
# the strplit command works really well in this context
## this way, multiple species are propagated into R as a simple vector  

def get_hypo_num(wildcards):
    """Get compared_to entries from hypotheses(.tsv) for each hypothesis. """
    num = hypotheses.loc[ (wildcards.hypothesis), 'hypothesis']
    return num

def get_hypo_name(wildcards):
    """Get compared_to entries from hypotheses(.tsv) for each hypothesis. """
    name = hypotheses.loc[ (wildcards.hypothesis), 'name']
    return name

# for both the species we want to check for expansion as well as those being compared to:
# if more than one species is compared to than we have to split the string based on ";" 

def get_exp_species(wildcards):
    """Get compared_to entries from hypotheses(.tsv) for each hypothesis. """
    exps_1 = hypotheses.loc[ (wildcards.hypothesis), 'expanded_in']
    if exps_1.count(";") > 0:
        exps_2 = str.split(exps_1, ";")
        return exps_2
    else:
        return exps_1
    return


#if more than one species is compared to than we have to split the string based on ";"
def get_com_species(wildcards):
    """Get compared_to entries from hypotheses(.tsv) for each hypothesis. """
    ct_1 = hypotheses.loc[ (wildcards.hypothesis), 'compared_to']
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
    path_prefix = 'FS/longest_isoforms/'
    suffix = '.fa'
    exp = hypotheses.loc[ (wildcards.hypothesis), 'expanded_in']
    ct = hypotheses.loc[ (wildcards.hypothesis), 'compared_to']
    if ct.count(";") > 0:
        ct = str.split(ct, ";")
        ct.append(exp)
        ct = [path_prefix + x + suffix for x in ct]
        return ct
    else:
        output = []
        output.append(exp)
        output.append(ct)
        output = [path_prefix + x + suffix for x in output]
        return output
    return


# get user defined expansion_factor for each hypothesis by parsing hypotheses.tsv
def get_expansion_factor(wildcards):
    return hypotheses.loc[wildcards.hypothesis]["min_expansion_factor"]


# get user defined Nmin_expanded_in & Nmin_compared_to for each hypothesis by parsing hypotheses.tsv
# explained: at least Nmin_expanded_in expanded species that are expanded in at least Nmin_compared_to compared_to species
def get_Nmin_expanded_in(wildcards):
    return hypotheses.loc[wildcards.hypothesis]["Nmin_expanded_in"]

def get_Nmin_compared_to(wildcards):
    return hypotheses.loc[wildcards.hypothesis]["Nmin_compared_to"]


# get user defined expanded_in_all_found & compared_to_all_found for each hypothesis by parsing hypotheses.tsv
# explained: at least Nmin_expanded_in expanded species that are expanded in at least Nmin_compared_to compared_to species
def get_expanded_in_all_found(wildcards):
    return hypotheses.loc[wildcards.hypothesis]["expanded_in_all_found"]

def get_compared_to_all_found(wildcards):
    return hypotheses.loc[wildcards.hypothesis]["compared_to_all_found"]


#SOLVING expansion checkpoint HERE!
#after rule FastTree in expansion_analysis.smk
def solve_expansion(wildcards):
    checkpoint_output = checkpoints.expansion.get(**wildcards).output[0]
    file_names = expand("tea/{hypothesis}/trees/{OG}.tree", hypothesis=wildcards.hypothesis, OG=glob_wildcards(os.path.join(checkpoint_output, "{OG}.txt")).OG)
    return file_names

