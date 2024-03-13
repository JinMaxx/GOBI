#!/usr/bin/python3

import subprocess
import tempfile
import json

datasets = "./ncbi_tools/datasets"
ncbi_data_dir = "./ncbi_data"
output_file_name = f"{ncbi_data_dir}/genome_summery.json"
download_log_file = f"{ncbi_data_dir}/ncbi_download.log"


def display_summery(accession_id: str):

    tmp = tempfile.NamedTemporaryFile()

    with open(tmp.name, 'w') as tmp_file:
        subprocess.run([
            datasets,
            "summary",
            "gene",
            "accession", str(accession_id),
        ], stdout=tmp_file)

    with open(tmp.name) as tmp_file:
        pretty_json = json.dumps(json.load(tmp_file), indent=2)
        with open(output_file_name, "w") as output_file:
            output_file.write(pretty_json)


def download_gene(taxonomy_id: int, accession_id: str) -> str:

    dirname = f"{ncbi_data_dir}/gene/{taxonomy_id}/{accession_id}"
    zip_filename = f"{dirname}.zip"

    print("downloading gene ", accession_id)
    subprocess.run([
        datasets,
        "download",
        "gene",
        "accession", accession_id,
        "--filename", zip_filename,
    ])

    print("unzipping ", zip_filename)
    subprocess.run([
        "unzip",
        zip_filename,
        "-d", dirname
    ])

    return dirname


def main():
    # getting all genomes
    genes_to_taxonomy_id_dict = {
        7227: ["NM_079617.3"]
    }
    # NM_079617.3 Timeout in Drosophilia M. (for testing)

    for taxonomy_id, accession_ids in genes_to_taxonomy_id_dict.items():
        for accession_id in accession_ids:
            display_summery(accession_id)
            download_gene(taxonomy_id, accession_id)


if __name__ == '__main__':
    main()
