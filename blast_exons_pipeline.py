#!/usr/bin/python3

#  blastn [-h] [-help] [-import_search_strategy filename]
#     [-export_search_strategy filename] [-task task_name] [-db database_name]
#     [-dbsize num_letters] [-gilist filename] [-seqidlist filename]
#     [-negative_gilist filename] [-negative_seqidlist filename]
#     [-taxids taxids] [-negative_taxids taxids] [-taxidlist filename]
#     [-negative_taxidlist filename] [-no_taxid_expansion]
#     [-entrez_query entrez_query] [-db_soft_mask filtering_algorithm]
#     [-db_hard_mask filtering_algorithm] [-subject subject_input_file]
#     [-subject_loc range] [-query input_file] [-out output_file]
#     [-evalue evalue] [-word_size int_value] [-gapopen open_penalty]
#     [-gapextend extend_penalty] [-perc_identity float_value]
#     [-qcov_hsp_perc float_value] [-max_hsps int_value]
#     [-xdrop_ungap float_value] [-xdrop_gap float_value]
#     [-xdrop_gap_final float_value] [-searchsp int_value] [-penalty penalty]
#     [-reward reward] [-no_greedy] [-min_raw_gapped_score int_value]
#     [-template_type type] [-template_length int_value] [-dust DUST_options]
#     [-filtering_db filtering_database]
#     [-window_masker_taxid window_masker_taxid]
#     [-window_masker_db window_masker_db] [-soft_masking soft_masking]
#     [-ungapped] [-culling_limit int_value] [-best_hit_overhang float_value]
#     [-best_hit_score_edge float_value] [-subject_besthit]
#     [-window_size int_value] [-off_diagonal_range int_value]
#     [-use_index boolean] [-index_name string] [-lcase_masking]
#     [-query_loc range] [-strand strand] [-parse_deflines] [-outfmt format]
#     [-show_gis] [-num_descriptions int_value] [-num_alignments int_value]
#     [-line_length line_length] [-html] [-sorthits sort_hits]
#     [-sorthsps sort_hsps] [-max_target_seqs num_sequences]
#     [-num_threads int_value] [-mt_mode int_value] [-remote] [-version]

import scrape_genes
import scrape_genome

import os
import subprocess
from pathlib import Path


# better to use taxonomy id that are specific for species
genes_to_taxonomy_id_dict = {
    7227: [41615]
}
# 50557 True Insects (~ 1.5 TB)
# 7227 Drosophilia M. (for testing)
# Accession Number:  [two-letter alphabetical prefix][six digits][.][version number]
# using gene_id instead as they are easier for handling. NCBI prefers those.


blastdb_directory_prefix = "./blastdb"  # os.path.abspath("./blastdb")


def download_genome(taxonomy_id) -> str:
    scrape_genome.display_summery(taxonomy_id)  # Download species genomes
    return scrape_genome.download_genome(taxonomy_id, include=["genome", "gff3"])


def build_blast_db(genomes_directory: str):

    data_directory = f"{genomes_directory}/ncbi_dataset/data"

    # example:
    # |   genomes_directory   |   after join    |     data      |    fasta_file
    # ./ncbi_data/genomes/7227/ncbi_dataset/data/GCF_000001215.4/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna
    for data in os.scandir(data_directory):  # taxonomy_id could relate to many species aka datasets
        if data.is_dir():
            for fasta_file in os.scandir(data):  # iterating trough fasta file(s) in .../data/*
                if fasta_file.is_file() and fasta_file.name.endswith(".fna"):
                    print(f"building blastdb from {fasta_file}")
                    output_filename = f"{blastdb_directory_prefix}/{taxonomy_id}/{Path(fasta_file.path).stem}.db"
                    subprocess.run([
                        "makeblastdb",
                        "-input_file", fasta_file,
                        "-input_type", "fasta"
                                       "-taxid", taxonomy_id,
                        "-title", taxonomy_id,
                        "-dbtype" "nt"
                        "-metadata_output_prefix", blastdb_directory_prefix,
                        "-out", output_filename
                    ])
                    yield output_filename


def get_genes(taxonomy_id: int, gene_ids: list):

    # Download all genes
    for gene_id in gene_ids:
        scrape_genes.display_summery(gene_id=gene_id)
        gene_directory = scrape_genes.download_gene(taxonomy_id, gene_id=gene_id)


        # >NM_079617.3 timeout [organism=Drosophila melanogaster] [GeneID=41615]
        # parse from fasta

        # here I have to download

        # TODO
        #   parse exons
        #    for this I will create a tuple (sequence, GFF_Record)
        #    with SeqIO i could parse the fasta file of the genome
        #    with start/end I could extract the sequence
        #    also I could use gff3 to get the location of the gene (type: "gene")


if __name__ == '__main__':

    #for taxonomy_id in genes_to_taxonomy_id_dict.keys():
    #    gene_directory = download_genome(taxonomy_id)
    #    build_blast_db(gene_directory)

    for taxonomy_id, accession_ids in genes_to_taxonomy_id_dict.items():
        get_genes(taxonomy_id, accession_ids)

    # TODO:
    #   blast exons against each db except its own
    #   if {blastdb_directory_prefix}/{taxonomy_id}.db == file: continue
    #   create table like: start, end, taxonomy, exon, gene
    #   also safe alignments


#    genomes_directory = os.fsencode(f"{genomes_directory}/ncbi_dataset/data")

# subprocess.run([
#     blastn,
#     '-query', query_file,
#     # '-query_loc', query_location,
#     '-strand', strand,
#     '-db', db,
#     # '-out', tmp.name,
#     '-outfmt', outfmt,
#     # '-taxids ', ','.join(taxids),
#     # '-perc_identity', perc_identity,
# ])

# Download Genes of interest. Preferably as Exons.

# Build Blast DB of those species.
# TODO make in blast.

# Blast Genes against each the db.

# Write Result in a special format.
# TODO this can be done in python. Parsing and such.
# Need start and end of those genes.
