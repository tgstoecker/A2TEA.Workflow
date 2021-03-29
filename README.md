# A2TEA
Automated Assessment of Trait specific Evolutionary Adaptations

This workflow combines RNA-seq analyses, differential gene expression, with evolutionary analyses - most notably gene family expansion events.
We use Orthofinder2 to infer gene duplication events and correlate these with significant physiological reaction patterns in the compared species.

At the moment for each species A2TEA requires as input RNA-Seq reads (both PE/SE possible) suitable for a differential expression experiments (control vs. treatment), either a genomic or transcriptomic fasta file + annotation (.gtf) as well as a peptide fasta.  


# Setup:
Install the Python 3 version of Miniconda.
you can get it here: https://docs.conda.io/en/latest/miniconda.html

Answer yes to the question whether conda shall be initialized and put into your PATH.

Then, you can install [mamba](https://github.com/QuantStack/mamba) (a faster replacement of conda in C++) with:

`conda install -c conda-forge mamba`

Download/Clone the current release of the A2TEA workflow into the directory.

The included environment.yaml file can be used to install all required software into an isolated Conda environment with a name of your choice - in the following we will call it "A2TEA":

`mamba env create --name A2TEA --file environment.yaml`

Activating the environment
To activate the snakemake-tutorial environment, execute

`conda activate A2TEA`

Now you can use the installed tools and our workflow without any software dependency issues.
For detailed options of snakemake see: https://snakemake.readthedocs.io/en/v5.5.1/executable.html

Should you want to remove the conda environment, execute
`conda env remove -n A2TEA`  


# General usage
## Recommended steps
1) Add fasta, annotation and fastq files to the input directories:  
  FASTAs and GTFs -> FS/  
  FASTQ files -> rawreads/  

2) Modify species.tsv and samples.tsv files
-> only when using cDNA FASTA and single-end reads for a species you NEED to add information to the fragment_length_mean column (single-end read length) as well as the standard deviation to the samples.tsv file  

3) Using the activated environment perform a dry-run and check for problems with:    
`snakemake --snakefile Snakefile_complete -np`  

4) Configure the config.yaml file to your needs  

5) Run A2TEA with (exchange XX for the amount of cores you can offer):  
`snakemake --snakefile Snakefile_complete --cores XX`  

# Important note on cDNA vs genomic fasta as choice for a species/ecotype/etc.:
cDNA input leads to kallisto as quantification software. This is much faster than using STAR and also requires much less resources.  
However, since our approach focuses on gene loci, the transcript-level quantification of kallisto needs to be aggregated to gene level as part of the differential expression analysis.  
This is done via the "makeTxDbFromGFF" function of the "GenomicFeatures" package in R (requires as input gff3 of gtf file).  
It works really well for the annotation files I have tested so far but this is i.m.O. a source of potential errors if e.g. non-standard annotations are used.  
In such cases, changes to the tximport.R script in scripts/ might be necessary - or one switches to the genomic FASTA/ STAR-based approach which directly quantifies at gene-level.  
Also: Currently working with genomic FASTA requires the use of a GTF annotation, since automatic changes to STAR and featureCounts options necessary for GFF usage is not included at the moment (will come in an upcoming update).  If you only possess a .gff/3 file you can solve this by converting it to .gtf with e.g. [gffread] (https://github.com/gpertea/gffread).  
  
  
# Some additional important pointers on usage:
1) Keep/add the "FS/" before the files in the species.tsv table 
2) Do NOT provide both a cDNA and genome fasta for a given species in the species.tsv file!  
  However using cDNA fasta for one species and genome fasta for another is totally fine.  
3) If you are using a genome fasta please also provide annotation .gtf file.    
If you are using cDNA fasta then also URL to the the annotation file suffices.  
4) Peptide & genome FASTA as well as GTF files shouldn't be compressed; cDNA FASTA should be gzipped
5) Fastq files should also be gzipped
6) The amount of cores specified on the command-line sets the maximum that snakemake will be able to use. If rule threads set in the Snakefile exceed this limit, they will be automatically scaled down. This means that if you diverge from my standard (= 24 cores) A2TEA will still run, however by modifying the threads for individual rules (in config.yaml / the Snakefile itself) you can improve performance for your particular computational setup.  
7) With "auto_isoform_filtering" you can choose whether to try an automated approach for filtering the peptide fastas for their longest isoform or doing this yourself before starting the workflow. If the option is not set to "YES" the filtering is skipped.
8) Although "all" software is installed old gnu installs might have problems;  
  the efficient usage of the split command as part of the gnu coreutils package (used during Orthofinder_split) requires a "newish" version,  
  namely one with the `-n, --number=chunks` option flag; currently conda is not able to install a proper version (this could be circumvented by distributing the workflow as a docker or singularity image)  
  -> if you can't install a newer version of coreutils split you can change the option "chunk_usage" to not "ON" and this (to be fair minor optimization) is not used.  
9) With "add_blast_hits" you can define the max number of additional best blast hits to include in the follow-up analyses. Tree visualizations often are more informative if we use more than an individual (H)OG.  

# Common resons for errors: 
- falsely formatted annotations; e.g. gene_id field is called different in some lines geneID  
- format of fasta files -> same lengths of lines and shorter; otherwise samtools faidx etc. won't work  


# To do:
- expansion analysis in R; currently only expansion in ONE species compared to ONE OR MANY is implemented  
- new release of Orthofinder (2.4.0), some changes to orthogroups output; should use this to overhaul the orthofinder parts and remove the unncessary bits  
- option for installation of all software and dependencies during runtime  - moved up in priorities because has become important to use multiple yaml files  
- related to prior point: stringr doesn't install in my environment
- checkpoints a new feature that is also going to be without alternative starting with snakemake v6 is the solution for using all the generated fasta records in a parallel manner in the subsequent steps; works now ;D    
- once that runs I saw that the new version of snakemake doesn't like the way I implemented the linking of the isoform filtering (currently commented out)..  
  
- muscle install and use instead of maaft; also mamba install -c anaconda gmp (The GNU multiprecision library)  
- remove usage of chunks? with coreutils split - e.g. on the cluster this is/can not be installed  
- ? picard/samtools for duplicate removal
- combination analyses of diff. exp. and orthologous groups  
-> R shiny overlay(?) integrating the data (+ GO analysis) and making them explorable (trees, etc.)  
(-> provide link to AHRD for researchers not possessing GO-/Annotation for their species of interest  

- more cleanup:  
-> add more options to the config.yaml files (e.g. trimmomatic options) so that all can be changed modified there  
-> add snakemake internal report  
-> restructering of A2TEA to modular layout  
-> this incl. seperate yaml files for the softwares and tools   
-> move the STAR index log to the logs/ directory  

  
  
### The workflow in its current form:
![Alt text](./latest_rulegraph.svg)

