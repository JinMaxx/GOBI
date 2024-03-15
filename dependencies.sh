#!/bin/bash

pip install --upgrade pip

pip install doit
pip install requests
pip install biopython

mkdir ./fasta
mkdir ./genomes
mkdir ./ncbi_tools
mkdir ./blast_db
mkdir ./proteins
mkdir ./proteins_filtered
mkdir ./proteins_aligned


# after here only linux:

cd ./ncbi_tools || exit 1
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets'
curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/dataformat'
chmod +x datasets dataformat
cd ..

# AliView
wget https://www.ormbunkar.se/aliview/downloads/linux/linux-version-1.28/aliview.tgz
tar -czvf aliview.tgz ./aliview/

# Mafft (Arch Repository)
# pacman -S mafft


exit 0