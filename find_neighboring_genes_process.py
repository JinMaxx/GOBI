#!/usr/bin/python3

# API Reference https://ncbi.github.io/blast-cloud/dev/api.html

# Payload
# QUERY: NM_001260153.1
# db: nucleotide
# QUERY_FROM:
# QUERY_TO:
# QUERYFILE: (binary)
# GENETIC_CODE: 1
# JOB_TITLE: test
# ADV_VIEW: true
# SUBJECTS:
# stype: nucleotide
# SUBJECTS_FROM:
# SUBJECTS_TO:
# SUBJECTFILE: (binary)
# DBTYPE: gc
# DATABASE: refseq_representative_genomes
# EQ_MENU: Drosophila melangaster (taxid:7227)
# NUM_ORG: 1
# EQ_TEXT:
# BLAST_PROGRAMS: megaBlast
# PHI_PATTERN:
# MAX_NUM_SEQ: 100
# SHORT_QUERY_ADJUST: on
# EXPECT: 0.05
# WORD_SIZE: 28
# HSP_RANGE_MAX: 0
# MATRIX_NAME: PAM30
# MATCH_SCORES: 1,-2
# GAPCOSTS: 0 0
# COMPOSITION_BASED_STATISTICS: 0
# FILTER: L
# REPEATS: repeat_9606
# FILTER: m
# TEMPLATE_LENGTH: 0
# TEMPLATE_TYPE: 0
# PSSM: (binary)
# I_THRESH:
# DI_THRESH:
# PSI_PSEUDOCOUNT:
# SHOW_OVERVIEW: true
# SHOW_LINKOUT: true
# GET_SEQUENCE: true
# FORMAT_OBJECT: Alignment
# FORMAT_TYPE: HTML
# ALIGNMENT_VIEW: Pairwise
# MASK_CHAR: 2
# MASK_COLOR: 1
# DESCRIPTIONS: 100
# ALIGNMENTS: 100
# LINE_LENGTH: 60
# NEW_VIEW: false
# NCBI_GI: false
# SHOW_CDS_FEATURE: false
# NUM_OVERVIEW: 100
# FORMAT_EQ_TEXT:
# FORMAT_ORGANISM:
# EXPECT_LOW:
# EXPECT_HIGH:
# PERC_IDENT_LOW:
# PERC_IDENT_HIGH:
# QUERY_INDEX:
# FORMAT_NUM_ORG: 1
# CONFIG_DESCR: ClustMemNbr,ClustComn,Ds,Sc,Ms,Ts,Cov,Eval,Idnt,AccLen,Acc
# CLIENT: web
# SERVICE: plain
# CMD: request
# PAGE: MegaBlast
# PROGRAM: blastn
# MEGABLAST: on
# RUN_PSIBLAST:
# WWW_BLAST_TYPE:
# TWO_HITS:
# DEFAULT_PROG: megaBlast
# DB_DISPLAY_NAME: refseq_representative_genomes
# ORG_DBS: orgDbsOnly_primer1_dbvers5
# SHOW_ORGANISMS: on
# DBTAXID:
# SAVED_PSSM:
# SELECTED_PROG_TYPE: megaBlast
# SAVED_SEARCH: true
# BLAST_SPEC:
# MIXED_DATABASE:
# QUERY_BELIEVE_DEFLINE:
# DB_DIR_PREFIX:
# CHECKSUM:
# USER_DATABASE:
# USER_WORD_SIZE:
# USER_MATCH_SCORES:
# USER_FORMAT_DEFAULTS:
# NO_COMMON:
# NUM_DIFFS: 1
# NUM_OPTS_DIFFS: 0
# UNIQ_DEFAULTS_NAME: A_SearchDefaults_1rjgQi_26qi_ducLcACGHCi5_GTMVb_8Rqcp
# PAGE_TYPE: BlastSearch
# USER_DEFAULT_PROG_TYPE: megaBlast
# USER_DEFAULT_MATCH_SCORES: 0



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


genes_to_taxonomy_id_dict = {
    7227: ["NM_079617.3"]
}
# 50557 True Insects (~ 1.5 TB)
# 7227 Drosophilia M. (for testing)
# Accession Number:  [two-letter alphabetical prefix][six digits][.][version number]

blastdb_directory_prefix = "./blastdb"  # os.path.abspath("./blastdb")

def build_blast_db(taxonomy_id):

    # Download species genomes
    scrape_genome.display_summery(taxonomy_id)
    genomes_directory = scrape_genome.download_genome(taxonomy_id)

    genomes_directory = os.path.join(genomes_directory, "ncbi_dataset", "data")
    print(genomes_directory)

    # example:
    # |   genomes_directory  |   os.path.join  |     data      |    fasta_file
    # ./ncbi_data/genome/7227/ncbi_dataset/data/GCF_000001215.4/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna
    for data in os.scandir(genomes_directory):  # taxonomy_id could relate to many species aka datasets
        if data.is_dir():
            for fasta_file in os.scandir(genomes_directory):  # iterating trough fasta file(s) in .../data/*
                if fasta_file.is_file() and fasta_file.name.endswith(".fna"):

                    print(f"building blastdb from {fasta_file}")
                    # build db
                    subprocess.run([
                        "makeblastdb",
                        "-input_file", fasta_file,
                        "-input_type", "fasta"
                        "-taxid", taxonomy_id,
                        "-title", taxonomy_id,
                        "-dbtype" "nt"
                        "-metadata_output_prefix", blastdb_directory_prefix,
                        "-out", f"{blastdb_directory_prefix}/{taxonomy_id}.db"
                    ])


def get_genes(taxonomy_id: int, accession_ids: list):

    # Download all genes
    for accession_id in accession_ids:
        scrape_genes.display_summery(accession_id)
        gene_directory = scrape_genes.download_gene(taxonomy_id, accession_id)

        # TODO
        #   parse exons
        #   blast exons against each db except its own
        #   if {blastdb_directory_prefix}/{taxonomy_id}.db == file: continue
        #   create table like: start, end, taxonomy, exon, gene
        #   also safe alignments


if __name__ == '__main__':
    for taxonomy_id in genes_to_taxonomy_id_dict.keys():
        build_blast_db(taxonomy_id)

    for taxonomy_id, accession_ids in genes_to_taxonomy_id_dict.items():
        get_genes(taxonomy_id, accession_ids)




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
