#!/usr/bin/python3

import subprocess

import scrape_genes
import scrape_genomes

if __name__ == "__main__":

    # installs dependencies
    subprocess.run([
        "bash",
        "./dependencies.sh",
    ])

    scrape_genes.main()  # downloads genes from all insects based on a list of gene names
    scrape_genomes.main()  # downloads all genomes of insects (only references)
