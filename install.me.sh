#!/bin/bash

# ***NOTE***
# Docker must already be installed!!

wkdir=$(pwd)
# Install nextflow
wget -qO- https://get.nextflow.io | bash
chmod +x nextflow
mv nextflow /usr/local/bin/.

# Pull pre-existing docker images from biocontainers & quay.io/biocontainers
docker pull quay.io/biocontainers/hisat2:2.2.1--h87f3376_4
docker pull biocontainers/fastqc:v0.11.9_cv7
docker pull quay.io/biocontainers/trim-galore:0.6.7--hdfd78af_0
docker pull biocontainers/samtools:v1.9-4-deb_cv1
docker pull quay.io/biocontainers/deeptools:3.5.1--py_0
docker pull biocontainers/htseq:v0.11.2-1-deb-py3_cv1

# Build custom docker images for hisat2+samtools, getphredencoding+TPM, & deseq2
docker build -t hisat2samtools ./custom_processes/hisat2samtools_dockerinfo
docker build -t python3only ./custom_processes/getphredencoding_tpm_dockerinfo
docker build -t deseq2 ./custom_processes/deseq2_dockerinfo
