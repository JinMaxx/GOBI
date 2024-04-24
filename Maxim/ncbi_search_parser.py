#!/usr/bin/python3
import os
import tempfile
from collections.abc import Generator

import pandas
from collections import defaultdict

import scrape_genome
from sequence_record_parser import SequenceRecordsBundle

debug = True


def read_all_to_sequence_records_bundle(directory: str) -> Generator[SequenceRecordsBundle]:
    for file in os.listdir(directory):
        for sequence_records_bundles in __read_to_sequence_records_bundles(f"{directory}/{file}"):
            for sequence_record_bundle in sequence_records_bundles:
                yield sequence_record_bundle


#    (tax_id, genome_id, _gene_id, SequenceRecordsBundle)
# -> Generator[(int, str, str, SequenceRecordsBundle)]:
def __read_to_sequence_records_bundles(file: str) -> Generator[list[SequenceRecordsBundle]]:

    dataframe = pandas.read_csv(file, sep='\t')
    dataframe.drop_duplicates(inplace=True)
    if debug: print(f"read table {file}:\n{dataframe.columns}\n{dataframe.head(5)}")
    dataframe = dataframe[['tax_id', 'GeneID']]

    # build dict for getting genes
    taxonomy_to_gene_ids_dict: dict[int, set[int]] = defaultdict(set)
    for index, row in dataframe.iterrows():
        taxonomy_to_gene_ids_dict[int(row['tax_id'])].add(int(row['GeneID']))
        if debug: print(int(row['tax_id']), row['GeneID'])

    for taxonomy_id, gene_ids in taxonomy_to_gene_ids_dict.items():
        if debug: print(f"processing taxonomy_id:[gene_ids]| {taxonomy_id}:{gene_ids}")

        with tempfile.TemporaryDirectory() as tmp_dir:
            # # # # # # # downloading genomes (to temporary directory) # # # # # # #
            scrape_genome.display_summery(taxonomy_id=taxonomy_id)
            genome_ids: list[str] = scrape_genome.download_genome(taxonomy_id=taxonomy_id, include=['genome', 'gff3'], genomes_dir=tmp_dir)
            # Could have multiple genomes downloaded (if taxonomy id is not specific)
            for genome_id in genome_ids:
                # genome_directory = f"{tmp_dir}/{taxonomy_id}/ncbi_dataset/data/{genome_id}"
                fasta_file, gff_file = scrape_genome.get_genome_files(taxonomy_id=taxonomy_id, genome_id=genome_id, directory=tmp_dir)

                if gff_file is None: continue
                # raise ValueError(f"gff file not found for {taxonomy_id}, {genome_id} in {tmp_dir}")

                # print(f"parsing {gff_file} with gene_ids: {gene_ids}")
                for gene_id in gene_ids:
                    # print(f"parsing for {gene_id}")
                    yield SequenceRecordsBundle.of_genes(
                            taxonomy_id=taxonomy_id,
                            genome_id=os.path.basename(genome_id),
                            gff_file=gff_file,
                            genome_fasta_file=fasta_file,
                            filter_types=["gene", "exon"],
                            gene_id=gene_id)


def read_all_to_sequence_records_bundle_neighbours(directory: str, region_width: int) -> Generator[SequenceRecordsBundle]:
    for file in os.listdir(directory):
        for sequence_records_bundles in __read_to_sequence_records_bundles_neighbours(f"{directory}/{file}", region_width):
            for sequence_record_bundle in sequence_records_bundles:
                yield sequence_record_bundle


# i was lazy so i copied the code from above
def __read_to_sequence_records_bundles_neighbours(file: str, region_width: int) -> Generator[list[SequenceRecordsBundle]]:

    dataframe = pandas.read_csv(file, sep='\t')
    dataframe.drop_duplicates(inplace=True)
    if debug: print(f"read table {file}:\n{dataframe.columns}\n{dataframe.head(5)}")
    dataframe = dataframe[['tax_id', 'GeneID']]

    # build dict for getting genes
    taxonomy_to_gene_ids_dict: dict[int, set[int]] = defaultdict(set)
    for index, row in dataframe.iterrows():
        taxonomy_to_gene_ids_dict[int(row['tax_id'])].add(int(row['GeneID']))
        if debug: print(int(row['tax_id']), row['GeneID'])

    for taxonomy_id, gene_ids in taxonomy_to_gene_ids_dict.items():
        if debug: print(f"processing taxonomy_id:[gene_ids]| {taxonomy_id}:{gene_ids}")

        with tempfile.TemporaryDirectory() as tmp_dir:
            # # # # # # # downloading genomes (to temporary directory) # # # # # # #
            scrape_genome.display_summery(taxonomy_id=taxonomy_id)
            genome_ids: list[str] = scrape_genome.download_genome(taxonomy_id=taxonomy_id, include=['genome', 'gff3'], genomes_dir=tmp_dir)
            # Could have multiple genomes downloaded (if taxonomy id is not specific)
            for genome_id in genome_ids:
                # genome_directory = f"{tmp_dir}/{taxonomy_id}/ncbi_dataset/data/{genome_id}"
                fasta_file, gff_file = scrape_genome.get_genome_files(taxonomy_id=taxonomy_id, genome_id=genome_id, directory=tmp_dir)

                if gff_file is None: continue
                # raise ValueError(f"gff file not found for {taxonomy_id}, {genome_id} in {tmp_dir}")

                # print(f"parsing {gff_file} with gene_ids: {gene_ids}")
                for gene_id in gene_ids:
                    # print(f"parsing for {gene_id}")
                    yield SequenceRecordsBundle.of_gene_neighbours(
                            taxonomy_id=taxonomy_id,
                            genome_id=os.path.basename(genome_id),
                            gff_file=gff_file,
                            genome_fasta_file=fasta_file,
                            filter_types=["gene"],
                            gene_id=gene_id,
                            region_width=region_width)


if __name__ == '__main__':
    for __sequence_record_bundle in read_all_to_sequence_records_bundle("./ncbi_search_test"):
        print(__sequence_record_bundle)
        exit(0)

# TODO: create parser for those tsv files
#  parser: panda?
#  parser -> genome_gene_id_dict
#  -> Sequence_record_bundle
#  genome -> filepath?
#  temporary genome_files
