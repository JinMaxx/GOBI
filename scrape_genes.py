#!/usr/bin/python3

import subprocess
import tempfile
import json
import os

_genes_dir = "./genes"

# no need to change that
_datasets = "./ncbi_tools/datasets"


def display_summery(gene_id: int = None,
                    accession_id: str = None,
                    genes_dir: str = _genes_dir):

    if gene_id is None and accession_id is None:
        raise ValueError("gene_id or accession_id must be specified.")

    _id, _keyword = (str(gene_id), "id") if gene_id is not None else (accession_id, "accession")

    tmp = tempfile.NamedTemporaryFile()
    with open(tmp.name, 'w') as tmp_file:
        subprocess.run([
            _datasets,
            "summary",
            "gene",
            _keyword, _id,
        ], stdout=tmp_file)

    summery_file = f"{genes_dir}/gene_summery.json"
    with open(tmp.name) as tmp_file:
        pretty_json = json.dumps(json.load(tmp_file), indent=2)
        with open(summery_file, "w") as output_file:
            output_file.write(pretty_json)


def download_gene(taxonomy_id: int,
                  gene_id: int = None,
                  accession_id: str = None,
                  genes_dir: str = _genes_dir) -> str:

    if gene_id is None and accession_id is None:
        raise ValueError("gene_id or accession_id must be specified.")

    _id, _keyword = (str(gene_id), "id") if gene_id is not None else (accession_id, "accession")

    dirname = f"{genes_dir}/{taxonomy_id}/{_id}"
    os.makedirs(os.path.dirname(dirname), exist_ok=True)
    zip_filename = f"{dirname}.zip"

    print("downloading gene ", _id)
    subprocess.run([
        _datasets,
        "download",
        "gene",
        _keyword, _id,
        "--filename", zip_filename
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
            display_summery(accession_id=accession_id)
            print(download_gene(taxonomy_id, accession_id=accession_id))
