#!/usr/bin/python3

import subprocess
import tempfile
import json

# getting all genomes
taxonomy_ids = [7227]
# 50557 True Insects (~ 1.5 TB)
# 7227 Drosophilia M. (for testing)

datasets = "./ncbi_tools/datasets"
ncbi_data_dir = "./ncbi_data"
# ncbi_data_dir = "/mnt/shared_large/ncbi_data"  # Here I have more memory, linking /other/path/ncbi_data to ./ncbi_data
output_file_name = f"{ncbi_data_dir}/genome_summery.json"
download_log_file = f"{ncbi_data_dir}/ncbi_download.log"
include = "genome,gff3"


def download_summery():

    for taxonomy_id in taxonomy_ids:

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

    for taxonomy_id in taxonomy_ids:

        dirname = f"{ncbi_data_dir}/{taxonomy_id}_genomes"
        zip_filename = f"{dirname}.zip"

        with open(download_log_file, 'a') as log_file:

            print("downloading genomes ", taxonomy_id)
            subprocess.run([
                datasets,
                "download",
                "genome",
                "taxon", str(taxonomy_id),
                "--reference",  # only downloading the reference genome
                "--dehydrated",
                "--filename", zip_filename,
                "--include", include
            ])

            print("unzipping ", zip_filename)
            subprocess.run([
                "unzip",
                zip_filename,
                "-d", dirname
            ])

            print("rehydrating ", dirname)
            subprocess.run([
                datasets,
                "rehydrate",
                "--directory", f"{dirname}/"
            ])


def main():
    # download_summery()
    download_genomes()


if __name__ == '__main__':
    main()
