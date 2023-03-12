FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="ad40e71b7f00dc1985a36a3b8e0a29fa4a3b3cc8522abe7ebe48b068ae2ba056"

RUN mamba install -c anaconda -c bioconda numpy python-newick --yes

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/cafe5_compile.yaml
#   prefix: /conda-envs/3ad9877d8ed81c644892fb2075504468
#   channels:
#     - conda-forge
#     - defaults
#     - anaconda
#   dependencies:
#     - _libgcc_mutex=0.1=conda_forge
#     - _openmp_mutex=4.5=1_gnu
#     - autoconf
#     - binutils=2.36.1=hdd6e379_1
#     - binutils_impl_linux-64=2.36.1=h193b22a_1
#     - binutils_linux-64=2.36=hf3e587d_32
#     - c-compiler=1.2.0=h7f98852_0
#     - compilers=1.2.0=ha770c72_0
#     - cxx-compiler=1.2.0=h4bd325d_0
#     - fortran-compiler=1.2.0=h1990efc_0
#     - gcc_impl_linux-64=9.3.0=h70c0ae5_19
#     - gcc_linux-64=9.3.0=hf25ea35_32
#     - gfortran_impl_linux-64=9.3.0=hc4a2995_19
#     - gfortran_linux-64=9.3.0=hdc58fab_32
#     - gxx_impl_linux-64=9.3.0=hd87eabc_19
#     - gxx_linux-64=9.3.0=h3fbe746_32
#     - kernel-headers_linux-64=2.6.32=he073ed8_14
#     - ld_impl_linux-64=2.36.1=hea4e1c9_1
#     - libgcc-devel_linux-64=9.3.0=h7864c58_19
#     - libgcc-ng=9.3.0=h2828fa1_19
#     - libgfortran-ng=9.3.0=hff62375_19
#     - libgfortran5=9.3.0=hff62375_19
#     - libgomp=9.3.0=h2828fa1_19
#     - libstdcxx-devel_linux-64=9.3.0=hb016644_19
#     - libstdcxx-ng=9.3.0=h6de172a_19
#     - sysroot_linux-64=2.12=he073ed8_14
#     - make
RUN mkdir -p /conda-envs/3ad9877d8ed81c644892fb2075504468
COPY workflow/envs/cafe5_compile.yaml /conda-envs/3ad9877d8ed81c644892fb2075504468/environment.yaml

# Conda environment:
#   source: workflow/envs/coreutils.yaml
#   prefix: /conda-envs/1301bba4ef77c8a38bbd2bcc61bacdcf
#   channels:
#     - conda-forge
#   dependencies:
#     - coreutils ==8.31
RUN mkdir -p /conda-envs/1301bba4ef77c8a38bbd2bcc61bacdcf
COPY workflow/envs/coreutils.yaml /conda-envs/1301bba4ef77c8a38bbd2bcc61bacdcf/environment.yaml

# Conda environment:
#   source: workflow/envs/deseq2_tximport.yaml
#   prefix: /conda-envs/8d2fb7bba52392c02883cdfdd167504c
#   channels:
#     - conda-forge
#     - bioconda
#     - r
#   dependencies:
#     - r-base >=3.6.0
#     - bioconductor-tximport >=1.14.0
#     - r-readr >=1.3.1
#     - bioconductor-genomeinfodbdata==1.2.6
#     - bioconductor-rhdf5 
#     - bioconductor-deseq2
#     - bioconductor-genomicfeatures
RUN mkdir -p /conda-envs/8d2fb7bba52392c02883cdfdd167504c
COPY workflow/envs/deseq2_tximport.yaml /conda-envs/8d2fb7bba52392c02883cdfdd167504c/environment.yaml

# Conda environment:
#   source: workflow/envs/diamond.yaml
#   prefix: /conda-envs/4dfa9f2a63759e03fdd1f98ad42cff79
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - diamond >=2.0.11
#     - samtools >=1.10
#     - coreutils >=8.31
RUN mkdir -p /conda-envs/4dfa9f2a63759e03fdd1f98ad42cff79
COPY workflow/envs/diamond.yaml /conda-envs/4dfa9f2a63759e03fdd1f98ad42cff79/environment.yaml

# Conda environment:
#   source: workflow/envs/expansion.yaml
#   prefix: /conda-envs/f077ee9ed9ca292e10af2b6ba2a97a35
#   channels:
#     - conda-forge
#     - bioconda
#     - r
#   dependencies:
#     - coreutils >=8.31
#   #  - r-base >=3.6.0
#     - r-base >=4.1.3
#   #  - r-readr ==1.3.1
#     - r-readr >=1.3.1
#     - r-stringr
#     - r-dplyr >=1.0
#     - r-plyr
#     - r-tidyr
#     - r-doparallel
#     - r-foreach
#   #  - r-tibble
#     - r-tibble >=3.1.7
#     - r-reshape2
#     - r-stringi
#     - ucsc-fasomerecords ==377
#     - muscle ==3.8.1551
#     - trimal >=1.4
#     - fasttree >=2.1.10
#     - bash >=5.1.16
RUN mkdir -p /conda-envs/f077ee9ed9ca292e10af2b6ba2a97a35
COPY workflow/envs/expansion.yaml /conda-envs/f077ee9ed9ca292e10af2b6ba2a97a35/environment.yaml

# Conda environment:
#   source: workflow/envs/fastqc.yaml
#   prefix: /conda-envs/08d4368302a4bdf7eda6b536495efe7d
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - fastqc ==0.11.9
RUN mkdir -p /conda-envs/08d4368302a4bdf7eda6b536495efe7d
COPY workflow/envs/fastqc.yaml /conda-envs/08d4368302a4bdf7eda6b536495efe7d/environment.yaml

# Conda environment:
#   source: workflow/envs/final_tea.yaml
#   prefix: /conda-envs/a5a1fbfe0b7a52815d65a25b7eeb527b
#   channels:
#     - conda-forge
#     - r
#     - bioconda
#   dependencies:
#     - r-base >=4.1.0
#   #  - r-readr ==1.3.1
#     - r-readr
#     - r-stringr
#     - r-stringi
#     - r-ape ==5.5
#     - r-dplyr ==1.0.10
#     - r-plyr
#     - r-tidyr
#     - r-tibble
#     - r-stringi
#     - bioconductor-genomeinfodbdata==1.2.6
#     - bioconductor-deseq2
#     - bioconductor-ggtree 
#     - bioconductor-biostrings
#     - bioconductor-survcomp
#     - r-seqinr
#     - r-upsetr
#     - r-cowplot
#     - r-ggplotify
#     - r-gtools
#   #==3.0.1
#   #  - bioconductor-biostrings 
#   #==2.60.0
#   #  - r-seqinr 
#   #==3.4_5
#   #  - r-upsetr 
#   #==1.4.0
#   #  - r-cowplot 
#   #==1.1.1
#   #  - r-ggplotify 
#   #==0.0.7
RUN mkdir -p /conda-envs/a5a1fbfe0b7a52815d65a25b7eeb527b
COPY workflow/envs/final_tea.yaml /conda-envs/a5a1fbfe0b7a52815d65a25b7eeb527b/environment.yaml

# Conda environment:
#   source: workflow/envs/func_annotation.yaml
#   prefix: /conda-envs/1f9823bb4ea54784adba0c0773084167
#   channels:
#     - conda-forge
#     - r
#   dependencies:
#     - r-base
#     - r-dplyr >=1.0.0
#     - r-readr
#   #  - r-base >=4.1.0
#   #  - r-dplyr >=1.0.0
#   #  - r-readr >=1.4.0
#     - r-stringr >=1.3.1
#   #  - r::r-stringi
#   #  - readline==8.1.0
#   #  - r-tibble
RUN mkdir -p /conda-envs/1f9823bb4ea54784adba0c0773084167
COPY workflow/envs/func_annotation.yaml /conda-envs/1f9823bb4ea54784adba0c0773084167/environment.yaml

# Conda environment:
#   source: workflow/envs/genetonic.yaml
#   prefix: /conda-envs/bdff40fb0e3851052f5c011c708696e6
#   channels:
#     - conda-forge
#     - bioconda
#     - r
#   dependencies:
#   #  - r-base>=4.1.0
#   #  - bioconductor-genomeinfodbdata==1.2.6
#   #  - bioconductor-go.db==3.13.0
#     - bioconductor-genetonic==2.2
#     - bioconductor-topgo
#     - bioconductor-deseq2
#     - r-dplyr
#     - r-tidyr
#     - r-stringr
#     - r-stringi
RUN mkdir -p /conda-envs/bdff40fb0e3851052f5c011c708696e6
COPY workflow/envs/genetonic.yaml /conda-envs/bdff40fb0e3851052f5c011c708696e6/environment.yaml

# Conda environment:
#   source: workflow/envs/gffread.yaml
#   prefix: /conda-envs/513cb7fbdf7a03430541b6b159b4883c
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - gffread =0.12.1
RUN mkdir -p /conda-envs/513cb7fbdf7a03430541b6b159b4883c
COPY workflow/envs/gffread.yaml /conda-envs/513cb7fbdf7a03430541b6b159b4883c/environment.yaml

# Conda environment:
#   source: workflow/envs/hypothesis_species_tree.yaml
#   prefix: /conda-envs/3c31eb44b46b6379c1567d42e31fdb5a
#   channels:
#     - conda-forge
#     - bioconda
#     - r
#   dependencies:
#     - r-base >=4.1.3
#     - r-readr >=1.3.1
#     - r-plyr
#     - r-stringr
#     - r-dplyr >=1.0
#     - r-tibble
#     - r-tidyr
#     - r-ape >=5.5
#     - r-stringi
RUN mkdir -p /conda-envs/3c31eb44b46b6379c1567d42e31fdb5a
COPY workflow/envs/hypothesis_species_tree.yaml /conda-envs/3c31eb44b46b6379c1567d42e31fdb5a/environment.yaml

# Conda environment:
#   source: workflow/envs/kallisto.yaml
#   prefix: /conda-envs/2687ce120ad54e66f96628e8052b5229
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - kallisto ==0.45.0
RUN mkdir -p /conda-envs/2687ce120ad54e66f96628e8052b5229
COPY workflow/envs/kallisto.yaml /conda-envs/2687ce120ad54e66f96628e8052b5229/environment.yaml

# Conda environment:
#   source: workflow/envs/longest_isoforms.yaml
#   prefix: /conda-envs/7514380e99d2eb02c197b4e252198117
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - python ==3.8.6
RUN mkdir -p /conda-envs/7514380e99d2eb02c197b4e252198117
COPY workflow/envs/longest_isoforms.yaml /conda-envs/7514380e99d2eb02c197b4e252198117/environment.yaml

# Conda environment:
#   source: workflow/envs/multiqc.yaml
#   prefix: /conda-envs/80d75427ea1f8e6cede9c5f260e8fc1e
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - multiqc ==1.13
RUN mkdir -p /conda-envs/80d75427ea1f8e6cede9c5f260e8fc1e
COPY workflow/envs/multiqc.yaml /conda-envs/80d75427ea1f8e6cede9c5f260e8fc1e/environment.yaml

# Conda environment:
#   source: workflow/envs/orthofinder.yaml
#   prefix: /conda-envs/f3a5d3d9a9b0e83518f7e8ea5c907726
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - python ==3.10
#     - orthofinder =2.5.4
#     - muscle =5.1
#     - mafft =7.490
#     - fasttree =2.1.11
RUN mkdir -p /conda-envs/f3a5d3d9a9b0e83518f7e8ea5c907726
COPY workflow/envs/orthofinder.yaml /conda-envs/f3a5d3d9a9b0e83518f7e8ea5c907726/environment.yaml

# Conda environment:
#   source: workflow/envs/samtools.yaml
#   prefix: /conda-envs/d88dfc5ec33fa1063ba19945425fc600
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - samtools >=1.10
RUN mkdir -p /conda-envs/d88dfc5ec33fa1063ba19945425fc600
COPY workflow/envs/samtools.yaml /conda-envs/d88dfc5ec33fa1063ba19945425fc600/environment.yaml

# Conda environment:
#   source: workflow/envs/star.yaml
#   prefix: /conda-envs/c7b04a2f0842ab902ff601e084929aef
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - star ==2.7.8a
RUN mkdir -p /conda-envs/c7b04a2f0842ab902ff601e084929aef
COPY workflow/envs/star.yaml /conda-envs/c7b04a2f0842ab902ff601e084929aef/environment.yaml

# Conda environment:
#   source: workflow/envs/subread.yaml
#   prefix: /conda-envs/1bc593050c78940738bc64ed099dc3be
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - subread ==2.0
RUN mkdir -p /conda-envs/1bc593050c78940738bc64ed099dc3be
COPY workflow/envs/subread.yaml /conda-envs/1bc593050c78940738bc64ed099dc3be/environment.yaml

# Conda environment:
#   source: workflow/envs/trimmomatic.yaml
#   prefix: /conda-envs/74871c16a29b0fef6bdb5c883e400542
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - trimmomatic ==0.36
RUN mkdir -p /conda-envs/74871c16a29b0fef6bdb5c883e400542
COPY workflow/envs/trimmomatic.yaml /conda-envs/74871c16a29b0fef6bdb5c883e400542/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/3ad9877d8ed81c644892fb2075504468 --file /conda-envs/3ad9877d8ed81c644892fb2075504468/environment.yaml && \
    mamba env create --prefix /conda-envs/1301bba4ef77c8a38bbd2bcc61bacdcf --file /conda-envs/1301bba4ef77c8a38bbd2bcc61bacdcf/environment.yaml && \
    mamba env create --prefix /conda-envs/8d2fb7bba52392c02883cdfdd167504c --file /conda-envs/8d2fb7bba52392c02883cdfdd167504c/environment.yaml && \
    mamba env create --prefix /conda-envs/4dfa9f2a63759e03fdd1f98ad42cff79 --file /conda-envs/4dfa9f2a63759e03fdd1f98ad42cff79/environment.yaml && \
    mamba env create --prefix /conda-envs/f077ee9ed9ca292e10af2b6ba2a97a35 --file /conda-envs/f077ee9ed9ca292e10af2b6ba2a97a35/environment.yaml && \
    mamba env create --prefix /conda-envs/08d4368302a4bdf7eda6b536495efe7d --file /conda-envs/08d4368302a4bdf7eda6b536495efe7d/environment.yaml && \
    mamba env create --prefix /conda-envs/a5a1fbfe0b7a52815d65a25b7eeb527b --file /conda-envs/a5a1fbfe0b7a52815d65a25b7eeb527b/environment.yaml && \
    mamba env create --prefix /conda-envs/1f9823bb4ea54784adba0c0773084167 --file /conda-envs/1f9823bb4ea54784adba0c0773084167/environment.yaml && \
    mamba env create --prefix /conda-envs/bdff40fb0e3851052f5c011c708696e6 --file /conda-envs/bdff40fb0e3851052f5c011c708696e6/environment.yaml && \
    mamba env create --prefix /conda-envs/513cb7fbdf7a03430541b6b159b4883c --file /conda-envs/513cb7fbdf7a03430541b6b159b4883c/environment.yaml && \
    mamba env create --prefix /conda-envs/3c31eb44b46b6379c1567d42e31fdb5a --file /conda-envs/3c31eb44b46b6379c1567d42e31fdb5a/environment.yaml && \
    mamba env create --prefix /conda-envs/2687ce120ad54e66f96628e8052b5229 --file /conda-envs/2687ce120ad54e66f96628e8052b5229/environment.yaml && \
    mamba env create --prefix /conda-envs/7514380e99d2eb02c197b4e252198117 --file /conda-envs/7514380e99d2eb02c197b4e252198117/environment.yaml && \
    mamba env create --prefix /conda-envs/80d75427ea1f8e6cede9c5f260e8fc1e --file /conda-envs/80d75427ea1f8e6cede9c5f260e8fc1e/environment.yaml && \
    mamba env create --prefix /conda-envs/f3a5d3d9a9b0e83518f7e8ea5c907726 --file /conda-envs/f3a5d3d9a9b0e83518f7e8ea5c907726/environment.yaml && \
    mamba env create --prefix /conda-envs/d88dfc5ec33fa1063ba19945425fc600 --file /conda-envs/d88dfc5ec33fa1063ba19945425fc600/environment.yaml && \
    mamba env create --prefix /conda-envs/c7b04a2f0842ab902ff601e084929aef --file /conda-envs/c7b04a2f0842ab902ff601e084929aef/environment.yaml && \
    mamba env create --prefix /conda-envs/1bc593050c78940738bc64ed099dc3be --file /conda-envs/1bc593050c78940738bc64ed099dc3be/environment.yaml && \
    mamba env create --prefix /conda-envs/74871c16a29b0fef6bdb5c883e400542 --file /conda-envs/74871c16a29b0fef6bdb5c883e400542/environment.yaml && \
    mamba clean --all -y

#Step 3: install stringi for R environments via R
#currently only way to install libraries correctly...

#expansion.yaml
RUN conda init && . ~/.bashrc && conda activate /conda-envs/f077ee9ed9ca292e10af2b6ba2a97a35 && \
R -e "install.packages('stringi', repos = 'http://cran.us.r-project.org'); if (!library(stringi, logical.return=T)) quit(status=10)" && \
conda deactivate

#final_tea.yaml
RUN conda init && . ~/.bashrc && conda activate /conda-envs/a5a1fbfe0b7a52815d65a25b7eeb527b && \
R -e "install.packages('stringi', repos = 'http://cran.us.r-project.org'); if (!library(stringi, logical.return=T)) quit(status=10)" && \
conda deactivate

#func_annotation.yaml
RUN conda init && . ~/.bashrc && conda activate /conda-envs/1f9823bb4ea54784adba0c0773084167 && \
R -e "install.packages('stringi', repos = 'http://cran.us.r-project.org'); if (!library(stringi, logical.return=T)) quit(status=10)" && \
conda deactivate

##genetonic.yaml
RUN conda init && . ~/.bashrc && conda activate /conda-envs/bdff40fb0e3851052f5c011c708696e6 && \
R -e "install.packages('stringi', repos = 'http://cran.us.r-project.org'); if (!library(stringi, logical.return=T)) quit(status=10)" && \
conda deactivate

##hypothesis_species_tree.yaml
RUN conda init && . ~/.bashrc && conda activate /conda-envs/3c31eb44b46b6379c1567d42e31fdb5a && \
R -e "install.packages('stringi', repos = 'http://cran.us.r-project.org'); if (!library(stringi, logical.return=T)) quit(status=10)" && \
conda deactivate

#Step 4: cleanup
RUN mamba clean --all -y
