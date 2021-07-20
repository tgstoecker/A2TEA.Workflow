#!/bin/bash

# option 1 (you have a modern machine / slash sudo rights)
# avoid interferences by mini/conda
# use - conda deactivate - to avoid using conda's installs of gnu libs

# Option 2 - THIS IS THE RECOMMENDED WAY!
# use conda/mamba and install the compilers metapackage
# can be tricky - cblas.h problems with system vs miniconda install...
#mamba install -c anaconda gcc_linux-64
#mamba install -c conda-forge compilers
#mamba install -c anaconda openblas-devel

git clone https://github.com/tgstoecker/CAFE5
cd CAFE5/
autoconf
./configure
make
