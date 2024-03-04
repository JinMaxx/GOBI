#!/bin/bash

pip install --upgrade pip

pip install doit
pip install requests
pip install biopython

mkdir ./fasta
mkdir ./ncbi_data
mkdir ./ncbi_tools


cd ./ncbi_tools || exit 1
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat'
chmod +x datasets dataformat
cd ..
