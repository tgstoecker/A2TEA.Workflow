#!/bin/bash

#download the test_data archive
wget \
  --no-check-certificate \
  --content-disposition \
  "https://uni-bonn.sciebo.de/s/vBGB09QBTAGAlMT/download"

#unpack
tar -xvf A2TEA.workflow_test_data.tar.gz

#change into the test directory
cd test_data

#move the read files into rawreads/
mv \
  B73_con_1.fq.gz B73_con_2.fq.gz \
  B73_con_3.fq.gz B73_drought_1.fq.gz \
  B73_drought_2.fq.gz B73_drought_3.fq.gz \
  Chilbo_con_1.1.fq.gz Chilbo_con_1.2.fq.gz \
  Chilbo_con_2.1.fq.gz Chilbo_con_2.2.fq.gz \
  Chilbo_drought_1.1.fq.gz Chilbo_drought_1.2.fq.gz \
  Chilbo_drought_2.1.fq.gz Chilbo_drought_2.2.fq.gz \
  Scarlett_control_1.1.fq.gz Scarlett_control_1.2.fq.gz \
  Scarlett_control_2.1.fq.gz Scarlett_control_2.2.fq.gz \
  Scarlett_control_3.1.fq.gz Scarlett_control_3.2.fq.gz \
  Scarlett_drought_1.1.fq.gz Scarlett_drought_1.2.fq.gz \
  Scarlett_drought_2.1.fq.gz Scarlett_drought_2.2.fq.gz \
  Scarlett_drought_3.1.fq.gz Scarlett_drought_3.2.fq.gz \
  ../rawreads/

#unzip the fasta & annotation files
gunzip *Oryza* *Hordeum* *Zea*

#move fasta & annotation files into FS directory
mv *Oryza* *Hordeum* *Zea* ../FS/
