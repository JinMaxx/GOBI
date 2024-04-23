#!/usr/bin/python3

import os
import sys
import json
import tempfile
import subprocess

from fnmatch import fnmatch
from io import StringIO
from typing import TypedDict

import Bio.SeqIO.FastaIO as FastaIO


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
                    include: list[str] = _include,      # genome_ids downloaded
                    genomes_dir: str = _genomes_dir,
                    assembly_level: str = None) -> list[str]:

    dirname = f"{genomes_dir}/{taxonomy_id}"
    zip_filename = f"{dirname}.zip"

    if os.path.isdir(dirname) and len(os.listdir(dirname)) > 0:
        # ignoring already downloaded genomes. just need to get all genome_ids

        dirname = f"{dirname}/ncbi_dataset/data"

        genome_ids: list[str] = list()
        for genome_id in os.listdir(dirname):
            if os.path.isdir(f"{dirname}/{genome_id}"): genome_ids.append(genome_id)

        return genome_ids

    print(f"downloading genomes {taxonomy_id} to {dirname}")

    params = [
        _datasets,
        "download",
        "genome",
        "taxon", str(taxonomy_id),
        "--reference",  # only downloading the reference genome
        "--dehydrated",
        "--filename", zip_filename,
        "--include", ",".join(include)
    ]

    if assembly_level is not None:
        params.append("--assembly-level")
        params.append(assembly_level)

    subprocess.run(params)

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

    return parse_file_catalog(file=f"{dirname}/ncbi_dataset/data/dataset_catalog.json")


def get_genome_seq(taxonomy_id: int,
                   genome_id: str,
                   seq_id: str) -> (str, str):

    fasta_file, gff_file = get_genome_files(taxonomy_id=taxonomy_id, genome_id=genome_id)

    with open(fasta_file, 'r') as file:
        for headline, genome_sequence in FastaIO.SimpleFastaParser(file):
            if seq_id in headline:
                # __genome_id_seq_id_seqence_dict[genome_id][seq_id] = (headline, genome_sequence)
                return headline, genome_sequence

    raise EOFError(f"Genome Sequence for {seq_id} not found in {fasta_file}")


def get_genome_files(taxonomy_id: int,
                     genome_id: str,
                     directory: str = _genomes_dir) -> (str, str):
    #                                -> fasta_file, (gff_file or None)

    genome_directory = f"{directory}/{taxonomy_id}/ncbi_dataset/data/{genome_id}"
    if not os.path.isdir(genome_directory): raise ValueError(f"{genome_directory} is not a directory.")
    else:

        fasta_file: str = None  # ignore
        gff_file: str = None  # ignore
        for file in os.listdir(genome_directory):  # iterating trough fasta file(s) in .../genome_directory/*
            file = f"{genome_directory}/{file}"
            if os.path.isfile(file) and (fnmatch(file, "*.fna") or fnmatch(file, "*.fasta")):
                fasta_file = file
            if os.path.isfile(file) and fnmatch(file, "*.gff"):
                gff_file = file

        if fasta_file is None: raise ValueError(f"fasta file not found in {genome_directory}")
        if gff_file is None: sys.stderr.write(f"gff file not found in {genome_directory}\n")

        return fasta_file, gff_file


class File(TypedDict):
    filePath: str
    fileType: str


class Assembly(TypedDict):
    accession: str
    files: File


class FileCatalog(TypedDict):
    apiVersion: str
    assemblies: list[Assembly]

    # {
    #   "apiVersion": "V2",
    #   "assemblies": [
    #     {
    #       "files": [
    #         {
    #           "filePath": "assembly_data_report.jsonl",
    #           "fileType": "DATA_REPORT"
    #         }
    #       ]
    #     },
    #     {
    #       "accession": "GCF_000001215.4",
    #       "files": [
    #         {
    #           "filePath": "GCF_000001215.4/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna",
    #           "fileType": "GENOMIC_NUCLEOTIDE_FASTA"
    #         },
    #         {
    #           "filePath": "GCF_000001215.4/genomic.gff",
    #           "fileType": "GFF3"
    #         }
    #       ]
    #     }
    #   ]
    # }


def parse_json_to_file_catalog(json_string: str) -> FileCatalog:
    json_string = StringIO(json_string)
    file_catalog: FileCatalog = json.load(json_string)
    return file_catalog


def parse_file_catalog(file: str) -> list[str]:

    genome_ids: list[str] = list()

    file_catalog: FileCatalog
    with open(file, 'r') as file_handler:
        file_catalog = parse_json_to_file_catalog(file_handler.read())
        print("file_catalog: ", file_catalog)

    for assembly in file_catalog['assemblies']:
        genome_id = assembly.get('accession')
        if genome_id is not None:
            genome_ids.append(genome_id)

    return genome_ids


if __name__ == '__main__':
    # getting all genomes
    taxonomy_ids = [7227]
    # 50557 True Insects (~ 1.5 TB)
    # 7227 Drosophilia M. (for testing)

    for taxonomy_id in taxonomy_ids:
        display_summery(taxonomy_id)
        print(f"downloaded genome_ids: {download_genome(taxonomy_id, include=_include)}")
        get_genome_seq(taxonomy_id=taxonomy_id, genome_id="GCF_000001215.4", seq_id="NT_037436.4")

