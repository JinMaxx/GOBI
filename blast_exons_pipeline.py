#!/usr/bin/python3
import tempfile

from Bio import Align


import scrape_genome
from sequence_record_parser import SequenceRecordsBundle, SequenceRecord

import os
import subprocess
from pathlib import Path
from fnmatch import fnmatch


# better to use taxonomy id that are specific for species
taxonomy_to_gene_ids_dict = {
    7227: [41615]
}
# 7227 Drosophilia M. (for testing)
# Accession Number:  [two-letter alphabetical prefix][six digits][.][version number]
# using gene_id instead as they are easier for handling. NCBI prefers those.


blastdb_directory = "./blast_db"
makeblastdb = "makeblastdb"
blastn = "blastn"


def download_genome(taxonomy_id) -> str:
    scrape_genome.display_summery(taxonomy_id)  # Download species genomes
    return scrape_genome.download_genome(taxonomy_id, include=["genome", "gff3"])


def build_blast_db(data_directory: str):
    for genome_directory in os.listdir(data_directory):  # taxonomy_id could relate to many species aka datasets
        # ignoring also files like:
        # .../data_directory/assembly_data_report.jsonl
        # .../data_directory/dataset_catalog.json
        genome_directory = f"{data_directory}/{genome_directory}"
        # print(genome_directory)
        if os.path.isdir(genome_directory):
            for fasta_file in os.listdir(genome_directory):  # iterating trough fasta file(s) in .../genome_directory/*
                fasta_file = f"{genome_directory}/{fasta_file}"
                # print(fasta_file)
                if os.path.isfile(fasta_file) and fnmatch(fasta_file, "*.fna") or fnmatch(fasta_file, "*.fasta"):
                    print(f"building blastdb from {fasta_file}")
                    genome_id = Path(fasta_file).stem
                    output_filename = f"{blastdb_directory}/{genome_id}"
                    subprocess.run([
                        makeblastdb,
                        "-in", fasta_file,
                        "-input_type", "fasta",
                        "-taxid", str(taxonomy_id),
                        "-title", genome_id,
                        "-dbtype", "nucl",
                        "-metadata_output_prefix", blastdb_directory,
                        "-out", output_filename,
                    ])
                    print(f"created blast db: {output_filename}")
                    yield genome_id


def get_exons_from_genome(data_directory: str, gene_ids: list[int]):

    # [(genome_id, gene_id, [GeneRecordsBundle])]
    per_genome_gene_records_bundle: list[(str, int, SequenceRecordsBundle)] = list()

    for genome_directory in os.scandir(data_directory):  # taxonomy_id could relate to many species aka datasets
        if genome_directory.is_dir():
            gff_file, genome_fasta_file = None, None
            genome_directory = f"{data_directory}/{genome_directory.name}"

            for file in os.listdir(genome_directory):
                if fnmatch(file, "*.gff"):
                    gff_file: str = f"{genome_directory}/{file}"
                    # print(f"gff_file:          {gff_file}")
                if fnmatch(file, "*.fna") or fnmatch(file, "*.fasta"):
                    genome_fasta_file: str = f"{genome_directory}/{file}"
                    # print(f"genome_fasta_file: {genome_fasta_file}")

            if gff_file is None or genome_fasta_file is None:
                raise ValueError(f"gff_file or fasta_file not found in {genome_directory}")

            for gene_id in gene_ids:
                per_genome_gene_records_bundle.append((
                    os.path.basename(genome_directory),
                    gene_id,
                    SequenceRecordsBundle.of_genes(
                        gff_file=gff_file,
                        genome_fasta_file=genome_fasta_file,
                        filter_types=["gene", "exon"],
                        gene_id=gene_id)
                ))

    return per_genome_gene_records_bundle


def global_alignment(reference_sequence: str, query_sequence: str):
    # TODO
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1

    for alignment in aligner.align(reference_sequence, query_sequence):
        print("Score = %.1f:" % alignment.score)
        print(alignment)


def local_alignment():
    # TODO
    pass


def find_genomic_region():
    # this function should find the approximate location of a gene with blast in another genome with BLAST
    # then a slice of that location is taken and searched with exons.
    # TODO
    pass


def blast_per_gene_records_bundle(per_genome_gene_records_bundle: list[(str, int, list[SequenceRecordsBundle])],
                                  db_list: list[str]):
    # Going to rewrite this function at it is not how I expected how it aligns exons.
    genome_id: str
    gene_id: int
    gene_records_bundle_list: list[SequenceRecordsBundle]
    for genome_id, gene_id, gene_records_bundle_list in per_genome_gene_records_bundle:

        gene_records_bundle: SequenceRecordsBundle
        for gene_records_bundle in gene_records_bundle_list:

            gene_record: SequenceRecord
            for gene_record in gene_records_bundle.sequence_records:

                # instead of using stdin using tmp files in fasta format
                query_input_file = tempfile.NamedTemporaryFile()
                with open(query_input_file.name, 'w') as query_input_file_handler:
                    query_input_file_handler.write(
                        f">gene_id: {gene_id}, type: {gene_record.gff_record.type}, "
                        f"reference_seqid: {gene_record.gff_record.seqid}, "
                        f"start: {gene_record.gff_record.start}, end: {gene_record.gff_record.end}\n"
                        f"{gene_record.sequence}"
                    )

                with open(query_input_file.name, 'r') as query_input_file_handler:
                    print(query_input_file_handler.read())

                for db in db_list:

                    output_file = tempfile.NamedTemporaryFile()
                    subprocess.run([
                        blastn,
                        '-task', 'megablast',
                        '-query', query_input_file.name,
                        '-strand', 'plus' if gene_record.gff_record.strand == '+' else 'minus',
                        '-db', db,
                        '-out', output_file.name,
                        '-outfmt', '0',  # pairwise 0, BLAST_XML 5
                        # '-window_size', '0'
                        # '-word_size', str(len(gene_record.sequence))
                        # '-taxids ', ','.join(taxids),
                        # '-perc_identity', perc_identity,
                    ])

                    with open(output_file.name, 'r') as output_file_handler:
                        print(output_file_handler.read())

                    # TODO: using BLAST to find genes similar in other species.
                    # using global alignment Needleman-Wunsch to search the exons in those regions


if __name__ == '__main__':

    for taxonomy_id, gene_id in taxonomy_to_gene_ids_dict.items():
        genomes_directory = download_genome(taxonomy_id)

        # example:
        # |   genomes_directory   |    data dir    |     genome      |    fasta_file (or .gff)
        # ./ncbi_data/genomes/7227/ncbi_dataset/data/GCF_000001215.4/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna
        data_directory = f"{genomes_directory}/ncbi_dataset/data"

        db_list = list()
        for db in build_blast_db(data_directory):
            db_list.append(f"{blastdb_directory}/{db}")

        per_genome_gene_records_bundle = get_exons_from_genome(data_directory, gene_id)
        #print(per_genome_gene_records_bundle)

        blast_per_gene_records_bundle(per_genome_gene_records_bundle, db_list)


# TODO:
#  Write Result in a special format.
#   create table like: taxonomy, start, end, exon, gene, reference_sequence, aligned_sequence