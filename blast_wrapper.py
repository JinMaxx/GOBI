#!/usr/bin/python3

import os
import json
import tempfile
import subprocess
from io import StringIO
from pathlib import Path
from typing import TypedDict

__blastdb_directory = "./blast_db"
__makeblastdb = "makeblastdb"
__blastn = "blastn"

# <blastdb_directory>/genome_id -> genome_id
__genome_db_dict: dict[str, str] = dict()
# This might look weird but there is one case where I want to get the genome_id.
# Most of the time I use it as a list of tuples.


# adding to genome_db_dict databases that where already build (used as a check later)
for db_file in os.listdir(__blastdb_directory):
    __genome_id = Path(db_file).stem
    __genome_db_dict[f"{__blastdb_directory}/{__genome_id}"] = __genome_id
for __genome_id in __genome_db_dict.values():
    print(f"found blast database: {__genome_id}")


class Params(TypedDict):
    expect: int
    sc_match: int
    sc_mismatch: int
    gap_open: int
    gap_extend: int
    filter: str


class QueryMasking(TypedDict):
    _from: int  # have to replace from as it is a python keyword
    to: int


class Description(TypedDict):
    id: str
    accession: str
    title: str
    taxid: str


class HSPS(TypedDict):  # most important type to iterate
    num: int
    bit_score: float
    score: int
    evalue: float
    identity: int
    query_from: int
    query_to: int
    query_strand: str
    hit_from: int
    hit_to: int
    hit_strand: str
    align_len: int
    gaps: int
    qseq: str
    hseq: str
    midline: str


class Hit(TypedDict):
    num: int
    description: list[Description]
    len: int
    hsps: list[HSPS]


class Search(TypedDict):
    query_id: str
    query_title: str
    query_len: int
    query_masking: list[QueryMasking]
    hits: list[Hit]


class Results(TypedDict):
    search: Search


class Report(TypedDict):
    program: str
    version: str
    reference: str
    search_target: dict[str, str]
    params: Params
    results: Results


class SimpleBlastReport:

    def __init__(self,
                 query_sequence: str,
                 hit_sequence: str,
                 mid_line: str,
                 gaps: int,
                 align_len: int,
                 identity: int,
                 query_from: int,
                 query_to: int,
                 query_strand: str,
                 hit_from: int,
                 hit_to: int,
                 hit_strand: str,
                 bit_score: float,
                 score: int,
                 evalue: float,
                 seqid: str,
                 genome_id: str):
        self.query_sequence = query_sequence
        self.hit_sequence = hit_sequence
        self.mid_line = mid_line
        self.gaps = gaps
        self.align_len = align_len
        self.identity = identity
        self.query_from = query_from
        self.query_to = query_to
        self.query_strand = query_strand
        self.hit_from = hit_from
        self.hit_to = hit_to
        self.hit_strand = hit_strand
        self.bit_score = bit_score
        self.score = score
        self.evalue = evalue
        self.seqid = seqid
        self.genome_id = genome_id

    @staticmethod
    def __transform_report_to_simple(report: Report, genome_id: str) -> list['SimpleBlastReport']:
        print(report)
        output_list: list[SimpleBlastReport] = list()
        hit: Hit
        for hit in report['results']['search']['hits']:
            hit_title = hit['description'][0]['title']
            for hsps in hit['hsps']:
                # print("hsps: \n", hsps)
                output_list.append(SimpleBlastReport.__transform_hsps_to_simple(hsps, hit_title, genome_id))
        return output_list

    @staticmethod
    def __transform_hsps_to_simple(hsps: HSPS, hit_title: str, genome_id: str) -> 'SimpleBlastReport':
        return SimpleBlastReport(
            query_sequence=hsps['qseq'],
            hit_sequence=hsps['hseq'],
            mid_line=hsps['midline'],
            gaps=hsps['gaps'],
            align_len=hsps['align_len'],
            identity=hsps['identity'],
            query_from=int(hsps['query_from']),
            query_to=int(hsps['query_to']),
            query_strand='+' if hsps['query_strand'].lower() == "plus" else '-',
            hit_from=int(hsps['hit_from']),
            hit_to=int(hsps['hit_to']),
            hit_strand='+' if hsps['hit_strand'].lower() == "plus" else '-',
            bit_score=hsps['bit_score'],
            score=hsps['score'],
            evalue=hsps['evalue'],
            seqid=hit_title.split(' ', maxsplit=1)[0],
            genome_id=genome_id
        )

    @staticmethod
    def parse_json(json_string: str, genome_id: str) -> list['SimpleBlastReport']:

        json_string = StringIO(json_string.replace('"from"', '"_from"'))

        output_simple_list: list[SimpleBlastReport] = list()

        data = json.load(json_string)
        output: dict[str, dict[str: list[Report]]] = data
        # {
        #   "BlastOutput2": [
        #       {
        #           "report": {
        #               ...

        for blast_output_list in output.values():
            for report_dict in blast_output_list:
                for report in report_dict.values():
                    output_simple_list.extend(SimpleBlastReport.__transform_report_to_simple(report, genome_id))

        return output_simple_list

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return (
            f"h{self.hit_strand}: {self.hit_sequence}\n"
            f"    {self.mid_line}\n"
            f"q{self.query_strand}: {self.query_sequence}\n"
            f"genome_id: {self.genome_id}, seqid: {self.seqid}\n"
            f"location:\n"
            f" {self.hit_from:10.0f} - {self.hit_to:10.0f}\n"
            f" {self.query_from:10.0f} - {self.query_to:10.0f}\n"
            f"identity: {self.identity}/{self.align_len}, gaps: {self.gaps}\n"
            f"bit_score: {self.bit_score}, score: {self.score}, e-value: {self.evalue}"
        )


def create_blast_db(genome_id: str, fasta_file: str, taxonomy_id: int = None, other_folder: str = None):

    directory: str

    if other_folder is not None:
        directory = other_folder
        db_set = set()
        for db_file in os.listdir(other_folder):
            __genome_id = Path(db_file).stem
            db_set.add(f"{other_folder}/{__genome_id}")
    else:
        directory = __blastdb_directory
        db_set = set(__genome_db_dict.keys())  # use default db

    output_filename = f"{directory}/{genome_id}"

    if output_filename in db_set:
        print(f"db already exists {output_filename}. skipping")  # checking if db was already created
    else:
        print(f"building blast db: {genome_id}")

        params: list[str] = [
            __makeblastdb,
            "-in", fasta_file,
            "-input_type", "fasta",
            "-title", genome_id,
            "-dbtype", "nucl",
            "-metadata_output_prefix", __blastdb_directory,
            "-out", output_filename,
        ]

        if taxonomy_id is not None:
            params.append("-taxid")
            params.append(str(taxonomy_id))

        subprocess.run(params)
        print(f"created blast db: {output_filename}")


def blast(query_fasta_file: str, strand: str = 'both', other_blast_db_directory: str = None) -> list[SimpleBlastReport]:

    if strand == '+':
        strand = 'plus'
    if strand == '-':
        strand = 'minus'

    db_dict: dict[str, str]

    if other_blast_db_directory is not None:
        db_dict = dict()
        for db_file in os.listdir(other_blast_db_directory):
            __genome_id = Path(db_file).stem
            db_dict[f"{other_blast_db_directory}/{__genome_id}"] = __genome_id
    else:
        db_dict = __genome_db_dict  # use default db

    for db in db_dict.keys():
        output_file = tempfile.NamedTemporaryFile()
        subprocess.run([
            __blastn,
            '-task', 'megablast',
            '-query', query_fasta_file,
            '-strand', strand,
            '-db', db,
            '-out', output_file.name,
            '-outfmt', '15',  # pairwise 0, BLAST_XML 5, Seqalign (JSON) 12, Single-file BLAST JSON 15

        ])

        with open(output_file.name, 'r') as output_file_handler:
            # print(output_file_handler.read())
            return SimpleBlastReport.parse_json(output_file_handler.read(), db_dict[db])


if __name__ == '__main__':

    # for testing
    for simple_blast_output in blast("./blast_test.fasta", "+"):
        print(simple_blast_output)

