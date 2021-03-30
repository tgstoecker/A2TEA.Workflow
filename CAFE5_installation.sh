#!/bin/bash

# option 1 (you have a modern machine / slash sudo rights)
# avoid interferences by mini/conda
# use - conda deactivate - to avoid using conda's installs of gnu libs
# then comment out the first line of option 2 and run the script

# Option 2 - THIS IS THE RECOMMENDED WAY!
# use conda/mamba and install the compilers metapackage
conda install -c conda-forge compilers
git clone https://github.com/tgstoecker/CAFE5
cd CAFE5/
autoconf
./configure
make
