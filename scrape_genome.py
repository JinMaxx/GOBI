#!/usr/bin/python3

import json
import os
import tempfile
import subprocess


_include = ["genome", "gff3"]

# no need to change those
_genomes_dir = "./genomes"
_datasets = "./ncbi_tools/datasets"


def display_summery(taxonomy_id: int,
                    genomes_dir: str = _genomes_dir):

    tmp = tempfile.NamedTemporaryFile()
    with open(tmp.name, 'w') as tmp_file:
        subprocess.run([
            _datasets,
            "summary",
            "genome",
            "taxon", str(taxonomy_id),
            "--reference"  # only downloading the reference genome
        ], stdout=tmp_file)

    summery_file = f"{genomes_dir}/genome_summery.json"
    with open(tmp.name) as tmp_file:
        pretty_json = json.dumps(json.load(tmp_file), indent=2)
        with open(summery_file, "w") as output_file:
            output_file.write(pretty_json)


def download_genome(taxonomy_id: int,
                    include: list[str] = _include,
                    genomes_dir: str = _genomes_dir) -> str:

    dirname = f"{genomes_dir}/{taxonomy_id}"
    zip_filename = f"{dirname}.zip"

    if os.path.isdir(dirname) and len(os.listdir(dirname)) > 0:
        return dirname  # ignoring already downloaded genomes

    print(f"downloading genomes {taxonomy_id} to {dirname}")
    subprocess.run([
        _datasets,
        "download",
        "genome",
        "taxon", str(taxonomy_id),
        "--reference",  # only downloading the reference genome
        "--dehydrated",
        "--filename", zip_filename,
        "--include", ",".join(include)
    ])

    print("unzipping ", zip_filename)
    subprocess.run([
        "unzip",
        zip_filename,
        "-d", dirname
    ])

    print("rehydrating ", dirname)
    subprocess.run([
        _datasets,
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
        download_genome(taxonomy_id, include=_include)


if __name__ == '__main__':
    main()
