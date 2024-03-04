#!/usr/bin/python3

import subprocess
import tempfile
import json

# getting all genomes
taxonomy_id = 7227  # 50557  # 7227

datasets = "./ncbi_tools/datasets"
ncbi_data_dir = "./ncbi_data"
output_file_name = f"{ncbi_data_dir}/genome_summery.json"
download_log_file = f"{ncbi_data_dir}/ncbi_download.log"

# datasets summary genome taxon 50557


def download_summery():

    tmp = tempfile.NamedTemporaryFile()

    with open(tmp.name, 'w') as tmp_file:
        subprocess.run([
            datasets,
            "summary",
            "genome",
            "taxon", str(taxonomy_id),
            "--reference"  # only downloading the reference genome
        ], stdout=tmp_file)

    with open(tmp.name) as tmp_file:
        pretty_json = json.dumps(json.load(tmp_file), indent=2)
        with open(output_file_name, "w") as output_file:
            output_file.write(pretty_json)


def download_genomes():

    dirname = f"{ncbi_data_dir}/{taxonomy_id}_genomes"
    zip_filename = f"{dirname}.zip"

    with open(download_log_file, 'a') as log_file:

        subprocess.run([
            datasets,
            "download",
            "genome",
            "taxon", str(taxonomy_id),
            "--reference",  # only downloading the reference genome
            "--dehydrated",
            "--filename", zip_filename,
            "--include", "genome,gff3"
        ], stdout=log_file)

        subprocess.run([
            "unzip",
            zip_filename,
            "-d", dirname
        ], stdout=log_file)

        subprocess.run([
            datasets,
            "rehydrate",
            "--directory", f"{dirname}/"
        ], stdout=log_file)


def main():
    # download_summery()
    download_genomes()


if __name__ == '__main__':
    main()
