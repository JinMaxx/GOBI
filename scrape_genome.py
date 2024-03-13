#!/usr/bin/python3

import json
import tempfile
import subprocess

datasets = "./ncbi_tools/datasets"
ncbi_data_dir = "./ncbi_data"
output_file_name = f"{ncbi_data_dir}/genome_summery.json"
download_log_file = f"{ncbi_data_dir}/ncbi_download.log"
# include = "genome,gff3"


def display_summery(taxonomy_id):

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


def download_genome(taxonomy_id: int) -> str:

    dirname = f"{ncbi_data_dir}/genome/{taxonomy_id}"
    zip_filename = f"{dirname}.zip"

    print("downloading genomes ", taxonomy_id)
    subprocess.run([
        datasets,
        "download",
        "genome",
        "taxon", str(taxonomy_id),
        "--reference",  # only downloading the reference genome
        "--dehydrated",
        "--filename", zip_filename,
        # "--include", include
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

    return dirname


def main():
    # getting all genomes
    taxonomy_ids = [7227]
    # 50557 True Insects (~ 1.5 TB)
    # 7227 Drosophilia M. (for testing)

    for taxonomy_id in taxonomy_ids:
        display_summery(taxonomy_id)
        download_genome(taxonomy_id)


if __name__ == '__main__':
    main()
