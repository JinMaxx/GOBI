#!/usr/bin/python3
import os
import tempfile
from collections.abc import Generator
from fnmatch import fnmatch

import pandas
from collections import defaultdict

import scrape_genome
from sequence_record_parser import SequenceRecordsBundle


def read_from_dir_to_sequence_records(directory: str = "./ncbi_search") -> Generator[(str, str, SequenceRecordsBundle)]:

    for file in os.listdir(directory):
        genome_directory = f"{directory}/{file}"
        for sequence_records in read_to_sequence_records(genome_directory):
            yield sequence_records


#                                            (tax_id, genome_id, _gene_id, SequenceRecordsBundle)
def read_to_sequence_records(file: str) -> Generator[(int, str, str, SequenceRecordsBundle)]:

    dataframe = pandas.read_csv(file, sep='\t')
    dataframe.drop_duplicates(inplace=True)
    print(f"read table {file}:\n{dataframe.columns}\n{dataframe.head(5)}")
    dataframe = dataframe[['tax_id', 'geneID']]

    # build dict for getting genes
    taxonomy_to_gene_ids_dict: dict[int, set[int]] = defaultdict(set)
    for index, row in dataframe.iterrows():
        taxonomy_to_gene_ids_dict[int(row['tax_id'])].add(int(row['geneID']))
        print(int(row['tax_id']), row['geneID'])

    for taxonomy_id, gene_ids in taxonomy_to_gene_ids_dict.items():
        print(f"processing taxonomy_id:[gene_ids]| {taxonomy_id}:{gene_ids}")

        with tempfile.TemporaryDirectory() as tmp_dir:
            # # # # # # # downloading genomes # # # # # # #
            scrape_genome.display_summery(taxonomy_id)  # Download species genomes (ignores if file already exist)
            data_directory = scrape_genome.download_genome(taxonomy_id, include=['genome', 'gff3'], genomes_dir=tmp_dir)
            #                   ./genomes/<taxonomy_id>
            data_directory = f"{data_directory}/ncbi_dataset/data"
            print(f"genome data for {taxonomy_id} in {data_directory}")

            # Could have multiple genomes downloaded (if taxonomy id is not specific)
            for _genome_id in os.listdir(data_directory):
                genome_directory = f"{data_directory}/{_genome_id}"

                # ignoring also files like:
                # .../<data_directory>/assembly_data_report.jsonl
                # .../<data_directory>/dataset_catalog.json
                if os.path.isdir(genome_directory):

                    fasta_file: str = None  # ignore
                    gff_file: str = None  # ignore
                    for file in os.listdir(genome_directory):  # iterating trough fasta file(s) in .../genome_directory/*
                        file = f"{genome_directory}/{file}"
                        # print(file)
                        if os.path.isfile(file) and (fnmatch(file, "*.fna") or fnmatch(file, "*.fasta")):
                            fasta_file = file
                        if os.path.isfile(file) and fnmatch(file, "*.gff"):
                            gff_file = file

                    if fasta_file is None: raise ValueError(f"fasta file not found in {genome_directory}")
                    else: print(f"for {_genome_id} found {fasta_file}")

                    # # # # # # # parsing genome_sequence_records # # # # # # #
                    if gff_file is None: raise ValueError(f"gff file not found in {genome_directory}")
                    else:
                        print(f"parsing {gff_file} with gene_ids: {gene_ids}")  # optional?
                        for _gene_id in gene_ids:
                            print(f"parsing for {_gene_id}")
                            yield (
                                taxonomy_id,
                                os.path.basename(_genome_id),
                                _gene_id,
                                SequenceRecordsBundle.of_genes(
                                    taxonomy_id=taxonomy_id,
                                    genome_id=os.path.basename(_genome_id),
                                    gff_file=gff_file,
                                    genome_fasta_file=fasta_file,
                                    filter_types=["gene", "exon"],
                                    gene_id=_gene_id)
                            )

# TODO: create parser for those tsv files
#  parser: panda?
#  parser -> genome_gene_id_dict
#  -> Sequence_record_bundle
#  genome -> filepath?
#  temporary genome_files