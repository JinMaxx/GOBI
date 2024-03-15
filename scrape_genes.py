#!/usr/bin/python3

import subprocess
import tempfile
import json
import os

_genes_dir = "./genes"

# no need to change that
_datasets = "./ncbi_tools/datasets"


def display_summery(accession_id: str,
                    genes_dir: str = _genes_dir):

    tmp = tempfile.NamedTemporaryFile()
    with open(tmp.name, 'w') as tmp_file:
        subprocess.run([
            _datasets,
            "summary",
            "gene",
            "accession", str(accession_id),
        ], stdout=tmp_file)

    summery_file = f"{genes_dir}/gene_summery.json"
    with open(tmp.name) as tmp_file:
        pretty_json = json.dumps(json.load(tmp_file), indent=2)
        with open(summery_file, "w") as output_file:
            output_file.write(pretty_json)


def download_gene(taxonomy_id: int,
                  accession_id: str,
                  genes_dir: str = _genes_dir) -> str:

    dirname = f"{genes_dir}/{taxonomy_id}/{accession_id}"
    os.makedirs(os.path.dirname(dirname), exist_ok=True)
    zip_filename = f"{dirname}.zip"

    print("downloading gene ", accession_id)
    subprocess.run([
        _datasets,
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


if __name__ == '__main__':
    genes_to_taxonomy_id_dict = {
        7227: ["NM_079617.3"]
    }
    # NM_079617.3 Timeout in Drosophilia M. (for testing)

    for taxonomy_id, accession_ids in genes_to_taxonomy_id_dict.items():
        for accession_id in accession_ids:
            display_summery(accession_id)
            print(download_gene(taxonomy_id, accession_id))
