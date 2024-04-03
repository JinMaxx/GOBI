#!/usr/bin/python3

import os
import re
import sys

import scrape_genome

from fnmatch import fnmatch
from collections import defaultdict
from tempfile import NamedTemporaryFile, TemporaryDirectory
from ncbi_search_parser import read_from_dir_to_sequence_records
from blast_wrapper import SimpleBlastReport, blast, create_blast_db
from sequence_record_parser import SequenceRecordsBundle, SequenceRecord

import Bio.SeqIO.FastaIO as FastaIO

# for the resulting .csv
_data_format = (
    # This is a format containing all parameters

    # "hit_from,hit_to,query_from,query_to,hit_gaps,query_gaps,align_from,align_to,"
    # "query_sequence,hit_sequence,mid_line,align_len,identity,query_strand,hit_strand,"
    # "bit_score,evalue,hit_seqid,hit_genome_id,query_genome_id,hit_taxonomy_id,query_taxonomy_id,full_query_sequence,"
    # "gff_seqid,gff_source,gff_type,gff_start,gff_end,gff_score,gff_strand,gff_phase,gff_attributes,"
    # gff_attribute:GeneID\n"

    "hit_from,hit_to,query_from,query_to,hit_gaps,query_gaps,align_from,align_to,"
    "align_len,identity,query_strand,hit_strand,"
    "bit_score,evalue,hit_seqid,hit_genome_id,query_genome_id,hit_taxonomy_id,query_taxonomy_id,"
    "gff_seqid,gff_source,gff_type,gff_start,gff_end,gff_score,gff_strand,gff_phase,"
    "gff_attribute:GeneID\n"

    # maybe find a way to also include query_taxid, hit_taxid
)

# better to use taxonomy id that are specific for species

#  shaggy

# species to look against:
taxonomies: list[int] = [
    # datasets:
    218467,  # Bark scorpion (Centruroides sculpturatus)
    114398,  # Common house spider (Parasteatoda tepidariorum)
    7091,    # Domestic silkworm (Bombyx mori)
    7460,    # Honey bee (Apis mellifera)
    13037,   # Monarch butterfly (Danaus plexippus)
    7227,    # Fruit fly (Drosophila melanogaster)
    13037,   # Monarch butterfly (Danaus plexippus)
    7111,    # Cabbage looper (Trichoplusia ni)
    13686,   # Red fire ant (Solenopsis invicta)
    6945,    # Black-legged tick (Ixodes scapularis)
    7070,    # Red flour beetle (Tribolium castaneum)
    110193,  # (Nicrophorus vespilloides)

    # interesting insectos to compare to:
    7119,    # Chinese oak silkmoth (Antheraea pernyi)
    189913,  # butterfly family Pieridae[white, yellow, sulphur] (Leptidea sinapis)

    # other insects:
    7165,    # African malaria mosquito (Anopheles gambiae)
    79782,   # Bed bug (Cimex lectularius)
]


__region_width: int = 300  # maybe higher 10.000

blastdb_directory = "./blast_db"
genomes_directory = "./genomes"
makeblastdb = "makeblastdb"
blastn = "blastn"
alignments_bundle_output_file = "./alignments_bundle_output.csv"

# Some global variables for easier handling.

# genome_id -> seqid -> (headline, sequence)    saving already read sequence for given seqid per genome_id
__genome_id_seqid_seqence_dict: dict[str, dict[str, (str, str)]] = defaultdict(dict)

# [(genome_id, gene_id, [GeneRecordsBundle])]
__per_genome_sequence_records_bundle: list[(str, int, SequenceRecordsBundle)] = list()

# building blast_db
for taxonomy_id in taxonomies:
    print(f"processing taxonomy_id: {taxonomy_id}")

    # # # # # # # downloading genomes # # # # # # #
    scrape_genome.display_summery(taxonomy_id)  # Download species genomes (ignores if file already exist)
    data_directory = scrape_genome.download_genome(taxonomy_id, include=['genome'])
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
            for file in os.listdir(genome_directory):  # iterating trough fasta file(s) in .../genome_directory/*
                file = f"{genome_directory}/{file}"
                # print(file)
                if os.path.isfile(file) and (fnmatch(file, "*.fna") or fnmatch(file, "*.fasta")): fasta_file = file

            if fasta_file is None:
                raise ValueError(f"fasta file not found in {genome_directory}")
            else:  # # # # # # building blast db # # # # # #
                print(f"for {_genome_id} found {fasta_file}")
                create_blast_db(_genome_id, fasta_file, taxonomy_id)


class Alignment:
    #         0      - region_width   align_from   hit_from      hit_to       align_to    + region_width      n
    # refer:  |-~~~--^----------------^------------|-------------|------------^-----------^---------------~~~-|
    # blast alignment:                             |||||||||||||||
    # query:                          |------------|-------------|------------|
    #                                 0            query_from    query_to     sequence_len
    #   can be negative!   align_from = hit_from - query_from        align_to = hit_to + (sequence_len - query_to)
    #   region_from = align_from - region_width >= 0                           region_to = align_to + region_width <= n

    __gff_attributes_pattern = re.compile(r"gff_attribute:\w+")

    def __init__(self,
                 hit_from: int, hit_to: int,
                 query_from: int, query_to: int,
                 hit_gaps: int, query_gaps: int,
                 align_from: int, align_to: int,
                 simple_blast_report: SimpleBlastReport,
                 query_sequence_record: SequenceRecord):
        self.hit_from = hit_from
        self.hit_to = hit_to
        self.query_from = query_from
        self.query_to = query_to
        self.hit_gaps = hit_gaps
        self.query_gaps = query_gaps
        self.align_from = align_from
        self.align_to = align_to
        self.simple_blast_report = simple_blast_report
        self.sequence_record = query_sequence_record

    # I need to parse the fasta file here with the genomic_id
    @staticmethod
    def from_simple_blast_report(query_sequence_record: SequenceRecord,
                                 simple_blast_report: SimpleBlastReport) -> 'Alignment':

        hit_from = simple_blast_report.hit_from
        hit_to = simple_blast_report.hit_to
        query_from = simple_blast_report.query_from
        query_to = simple_blast_report.query_to

        # location flipped based on strand if - not + because we want to have the location on the reference genome
        if hit_from > hit_to: hit_from, hit_to = hit_to, hit_from
        if query_from > query_to: query_from, query_to = hit_to, hit_from

        # they begin counting from 1 (but x_to is correct and only after flipping)
        hit_from -= 1
        query_from -= 1

        # counting gaps on both strands
        hit_gaps = simple_blast_report.hit_sequence.count('-')
        query_gaps = simple_blast_report.query_sequence.count('-')

        # calculating alignment position on genome
        align_from = hit_from - query_from
        align_to = hit_to + (len(query_sequence_record.sequence) - query_to) + hit_gaps - query_gaps

        # print(f"align_from: {align_from}")
        # print(f"align_to: {align_to}")

        # headline, sequence = get_genome_seq(simple_blast_report.genome_id, simple_blast_report.seqid)
        # print(f"referencing {simple_blast_report.genome_id} {simple_blast_report.seqid}[{hit_from}:{hit_to}]: {sequence[hit_from:hit_to]}")

        # TODO: some chromosomes are split by multiples sequences with seqids.
        #  need to think about whaht to do when using regions as they could span over multiple sequences...

        return Alignment(hit_from=hit_from, hit_to=hit_to,
                         query_from=query_from, query_to=query_to,
                         hit_gaps=hit_gaps, query_gaps=query_gaps,
                         align_from=align_from,
                         align_to=align_to,
                         simple_blast_report=simple_blast_report,
                         query_sequence_record=query_sequence_record)

    def __get_region(self, region_width: int, ref_genome_seq_len: int) -> (int, int):

        # ref_genome_seq = self.__get_reference_genome_seq()

        # applying region_width
        region_from = self.align_from - region_width
        region_to = self.align_to + region_width

        # checking for boundaries
        if region_from <= 0: region_from = 0
        if region_to >= ref_genome_seq_len - 1: region_to = ref_genome_seq_len - 1  # TODO check if last index is correct

        return region_from, region_to

    def as_region_fasta_entry(self, region_width: int) -> (str, str):
        headline: str
        sequence: str
        headline, sequence = get_genome_seq(genome_id=self.simple_blast_report.genome_id,
                                            seqid=self.simple_blast_report.seqid)
        region_from, region_to = self.__get_region(region_width, len(sequence))
        # print(f"headline: {headline}")
        # print(f"seqid: {self.simple_blast_report.seqid}")
        # if self.simple_blast_report.seqid in headline:
        if headline.startswith(self.simple_blast_report.seqid):
            second_split = headline.split(' ', maxsplit=1)[1]
            headline = f"{self.simple_blast_report.seqid}|{region_from}:{region_to}"
            if second_split is not None:
                headline = f"{headline} {second_split}"
            # headline.replace(self.simple_blast_report.seqid,
            #                  f"{self.simple_blast_report.seqid}|{region_from}:{region_to}")
            # print(f"parsed: {headline}")
        else:
            raise ValueError(f"headline does not start with {self.simple_blast_report.seqid}\n{headline}")
        sequence = sequence[region_from:region_to]
        return headline, sequence

    def as_table_record(self, data_format: str) -> str:

        data_format = data_format.lower()  # yeah... there might be a more elegant solution to this

        # alignment
        if "hit_from" in data_format: data_format = data_format.replace("hit_from", str(self.hit_from))
        if "hit_to" in data_format: data_format = data_format.replace("hit_to", str(self.hit_to))
        if "query_from" in data_format: data_format = data_format.replace("query_from", str(self.query_from))
        if "query_to" in data_format: data_format = data_format.replace("query_to", str(self.query_to))
        if "hit_gaps" in data_format: data_format = data_format.replace("hit_gaps", str(self.hit_gaps))
        if "query_gaps" in data_format: data_format = data_format.replace("query_gaps", str(self.query_gaps))
        if "align_from" in data_format: data_format = data_format.replace("align_from", str(self.align_from))
        if "align_to" in data_format: data_format = data_format.replace("align_to", str(self.align_to))

        # simple blast report
        if "query_sequence" in data_format:
            data_format = data_format.replace("query_sequence", str(self.simple_blast_report.query_sequence))
        if "hit_sequence" in data_format:
            data_format = data_format.replace("hit_sequence", str(self.simple_blast_report.hit_sequence))
        if "mid_line" in data_format:
            data_format = data_format.replace("mid_line", str(self.simple_blast_report.mid_line))
        if "align_len" in data_format:
            data_format = data_format.replace("align_len", str(self.simple_blast_report.align_len))
        if "identity" in data_format:
            data_format = data_format.replace("identity", str(self.simple_blast_report.identity))
        if "query_strand" in data_format:
            data_format = data_format.replace("query_strand", str(self.simple_blast_report.query_strand))
        if "hit_strand" in data_format:
            data_format = data_format.replace("hit_strand", str(self.simple_blast_report.hit_strand))
        if "bit_score" in data_format:
            data_format = data_format.replace("bit_score", str(self.simple_blast_report.bit_score))
        if "evalue" in data_format:
            data_format = data_format.replace("evalue", str(self.simple_blast_report.evalue))
        if "hit_seqid" in data_format:
            data_format = data_format.replace("hit_seqid", str(self.simple_blast_report.seqid))
        if "hit_genome_id" in data_format:
            data_format = data_format.replace("hit_genome_id", str(self.simple_blast_report.genome_id))
        if "hit_taxonomy_id" in data_format:
            data_format = data_format.replace("hit_taxonomy_id", str(self.simple_blast_report.taxonomy_id))

        # sequence record
        if "query_genome_id" in data_format:
            data_format = data_format.replace("query_genome_id", str(self.sequence_record.get_genome_id()))
        if "query_taxonomy_id" in data_format:
            data_format = data_format.replace("query_taxonomy_id", str(self.sequence_record.get_taxonomy_id()))
        if "full_query_sequence" in data_format:
            data_format = data_format.replace("full_query_sequence", str(self.sequence_record.sequence))

        # gff record
        if "gff_seqid" in data_format:
            data_format = data_format.replace("seqid", str(self.sequence_record.gff_record.seqid))
        if "gff_source" in data_format:
            data_format = data_format.replace("gff_source", str(self.sequence_record.gff_record.source))
        if "gff_type" in data_format:
            data_format = data_format.replace("gff_type", str(self.sequence_record.gff_record.type))
        if "gff_start" in data_format:
            data_format = data_format.replace("gff_start", str(self.sequence_record.gff_record.start))
        if "gff_end" in data_format:
            data_format = data_format.replace("gff_end", str(self.sequence_record.gff_record.end))
        if "gff_score" in data_format:
            data_format = data_format.replace("gff_score", str(self.sequence_record.gff_record.score))
        if "gff_strand" in data_format:
            data_format = data_format.replace("gff_strand", str(self.sequence_record.gff_record.strand))
        if "gff_phase" in data_format:
            data_format = data_format.replace("gff_phase", str(self.sequence_record.gff_record.phase))

        for match in Alignment.__gff_attributes_pattern.finditer(data_format):
            attribute = match.group(0).split(':', maxsplit=1)[1]
            attribute_value = self.sequence_record.gff_record.get_attribute(attribute)
            if attribute_value is None: sys.stderr.write(f"No value found for attribute {attribute}")
            data_format = data_format[:match.start(0)] + str(attribute_value) + data_format[match.end(0):]

        return data_format

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        simple_blast_report_str = str(self.simple_blast_report).replace('\n', '\n\t')
        gff_record_str = str(self.sequence_record.gff_record).replace('\n', '\n\t')
        return (
            f"alignment: {self.align_from}-{self.align_from}\n"
            f"query: {self.query_from}-{self.query_to}, gaps: {self.query_gaps}"
            f"refer: {self.hit_from}-{self.hit_to}, gaps: {self.hit_gaps}"
            "simple_blast_report:\n"
            f"\t{simple_blast_report_str}\n"
            "sequence_record.gff_record:\n"
            f"\t{gff_record_str}"
        )


class AlignmentsBundle:  # remove overlaps and concatenate regions?

    # genome_id -> seqid -> (genome_sequence, list[gene_region])
    def __init__(self):
        self.genome_alignments_dict: dict[str, dict[str, list[Alignment]]] = defaultdict(dict)

    @staticmethod
    def from_simple_blast_reports(sequence_record: SequenceRecord,
                                  simple_blast_reports: list[SimpleBlastReport]) -> 'AlignmentsBundle':

        alignment_bundle: AlignmentsBundle = AlignmentsBundle()

        for simple_blast_report in simple_blast_reports:

            # because i cant defaultdict(defaultdict(list))
            if alignment_bundle.genome_alignments_dict[simple_blast_report.genome_id].get(
                    simple_blast_report.seqid) is None:
                alignment_bundle.genome_alignments_dict[simple_blast_report.genome_id][
                    simple_blast_report.seqid] = list()

            (alignment_bundle.genome_alignments_dict[simple_blast_report.genome_id][simple_blast_report.seqid].
             append(Alignment.from_simple_blast_report(query_sequence_record=sequence_record,
                                                       simple_blast_report=simple_blast_report)))

        return alignment_bundle

    def regions_as_fasta(self, region_width: int) -> list[NamedTemporaryFile]:
        # one file per genome
        genomes_regions_tmp_fasta_files: list[NamedTemporaryFile] = list()

        genome_id: str
        alignments_dict: dict[str, list[Alignment]]
        for genome_id, alignments_dict in self.genome_alignments_dict.items():

            tmp_fasta_file = NamedTemporaryFile()
            genomes_regions_tmp_fasta_files.append(tmp_fasta_file)

            alignments: list[Alignment]
            for alignments in alignments_dict.values():

                alignment: Alignment
                for alignment in alignments:
                    headline, sequence = alignment.as_region_fasta_entry(region_width)
                    to_tmp_fasta(headline, sequence, tmp_fasta_file)

        return genomes_regions_tmp_fasta_files

    def append(self, alignments_bundle: 'AlignmentsBundle'):
        alignments_dict: dict[str, list[Alignment]]
        for genome_id, alignments_dict in alignments_bundle.genome_alignments_dict.items():

            alignments: list[Alignment]
            for seqid, alignments in alignments_dict.items():

                if self.genome_alignments_dict[genome_id].get(seqid) is None:
                    self.genome_alignments_dict[genome_id][seqid] = list()

                self.genome_alignments_dict[genome_id][seqid].extend(alignments)

    def __repr__(self):
        return self.__str__()

    def __str__(self):

        output_string: str = str()

        alignments_dict: dict[str, list[Alignment]]
        for genome_id, alignments_dict in self.genome_alignments_dict.items():

            output_string = (
                f"{output_string}\n"
                f"{genome_id}:\n"
            )

            alignments: list[Alignment]
            for seqid, alignments in alignments_dict.items():

                output_string = (
                    f"{output_string}\n"
                    f"\n\t{seqid}:\n\n"
                )

                alignment: Alignment
                for alignment in alignments:
                    alignment_str = str(alignment).replace('\n', '\n\t\t')
                    output_string = (
                        f"{output_string}"
                        f"\t\t{alignment_str}\n\n"
                    )

        return output_string

    def write_table_to_file(self, file: str, format: str):

        with open(file, 'a') as file_handler:

            alignments_dict: dict[str, list[Alignment]]
            for _, alignments_dict in self.genome_alignments_dict.items():

                alignments: list[Alignment]
                for _, alignments in alignments_dict.items():

                    alignment: Alignment
                    for alignment in alignments:
                        file_handler.write(alignment.as_table_record(data_format=format))


def get_genome_seq(genome_id: str,
                   seqid: str) -> (str, str):
    if __genome_id_seqid_seqence_dict[genome_id].get(seqid) is not None:
        return __genome_id_seqid_seqence_dict[genome_id][seqid]
    else:
        genome_directory = f"{data_directory}/{genome_id}"
        if os.path.isdir(genome_directory):
            fasta_file: str = None  # ignore
            for file in os.listdir(genome_directory):  # iterating trough fasta file(s) in .../genome_directory/*
                file = f"{genome_directory}/{file}"
                if os.path.isfile(file) and (fnmatch(file, "*.fna") or fnmatch(file, "*.fasta")):
                    fasta_file = file

            if fasta_file is None:
                raise ValueError(f"fasta file not found in {genome_directory}")
            else:
                with open(fasta_file, 'r') as file:
                    for headline, genome_sequence in FastaIO.SimpleFastaParser(file):
                        if seqid in headline:
                            __genome_id_seqid_seqence_dict[genome_id][seqid] = (headline, genome_sequence)
                            return headline, genome_sequence

                raise EOFError(f"Genome Sequence for {seqid} not found in {fasta_file}")


def to_tmp_fasta(headline: str,  # without leading >
                 sequence: str,
                 tmp_fasta_file: NamedTemporaryFile = None) -> NamedTemporaryFile:
    if tmp_fasta_file is None: tmp_fasta_file = NamedTemporaryFile()  # new file, otherwise append
    with open(tmp_fasta_file.name, 'a') as tmp_fasta_file_handler:
        tmp_fasta_file_handler.write(
            f">{headline}\n"
            f"{sequence}\n"
        )
    return tmp_fasta_file


def __sequence_record_to_tmp_fasta(sequence_record: SequenceRecord):
    return to_tmp_fasta(headline=(
        f"gene_id: {_gene_id}, type: {sequence_record.gff_record.type}, "
        f"reference_seqid: {sequence_record.gff_record.seqid}, "
        f"start: {sequence_record.gff_record.start}, end: {sequence_record.gff_record.end}"),
        sequence=sequence_record.sequence)


def gene_alignment(sequence_record: SequenceRecord) -> list[SimpleBlastReport]:
    query_input_file = __sequence_record_to_tmp_fasta(sequence_record)
    return blast(query_fasta_file=query_input_file.name)
    #           strand=sequence_record.gff_record.strand


# I might have to create multiple temporary db_tables
def exon_alignment_to_regions(genome_id: str, region_width: int,
                              sequence_record: SequenceRecord,
                              alignments_bundle: AlignmentsBundle) -> list[SimpleBlastReport]:
    # using tmp files in fasta format instead of using stdin
    query_input_file = __sequence_record_to_tmp_fasta(sequence_record)

    tmp_blast_db_dir: TemporaryDirectory = TemporaryDirectory()

    for genome_region_fasta_file in alignments_bundle.regions_as_fasta(region_width=region_width):
        create_blast_db(genome_id=genome_id,
                        fasta_file=genome_region_fasta_file.name,
                        other_folder=tmp_blast_db_dir.name)

    return blast(query_fasta_file=query_input_file.name, other_blast_db_directory=tmp_blast_db_dir.name)
    #            strand=sequence_record.gff_record.strand


if __name__ == '__main__':

    with open(alignments_bundle_output_file, 'w') as _output_file_handler:
        _output_file_handler.write(_data_format)  # overwriting with header

    # _alignments_bundle_output: AlignmentsBundle = AlignmentsBundle()  # this can be my memory killer...

    _genome_id: str
    _gene_id: int  # when there are many different reference genomes
    _sequence_records_bundle_list: list[SequenceRecordsBundle]
    for _taxonomy_id, _genome_id, _gene_id, _sequence_records_bundle_list \
            in read_from_dir_to_sequence_records(directory="./ncbi_search_conservative"):

        _sequence_records_bundle: SequenceRecordsBundle
        for _sequence_records_bundle in _sequence_records_bundle_list:
            # for 1 gene there is 1 sequence_record_bundle with many sequences for exons and one sequence for the gene

            # Getting the AlignmentsBundle just for type == "gene"
            _alignments_bundle: AlignmentsBundle = None  # ignore
            _sequence_record: SequenceRecord
            for _sequence_record in _sequence_records_bundle:
                if _sequence_record.gff_record.type == "gene":
                    _alignments_bundle = (AlignmentsBundle  # finding the region where the gene exists
                                          .from_simple_blast_reports(sequence_record=_sequence_record,
                                                                     simple_blast_reports=gene_alignment(_sequence_record)))
                    break                                                               # ^ alignment is done here
            if _alignments_bundle is None: raise ValueError("No gene annotation found in sequence_records")
            # _alignments_bundle_output.append(_alignments_bundle)

            for _sequence_record in _sequence_records_bundle.sequence_records:

                if _sequence_record.gff_record.type == "exon":
                    _simple_blast_reports: list[SimpleBlastReport] = (
                        exon_alignment_to_regions(genome_id=_genome_id,
                                                  region_width=__region_width,
                                                  sequence_record=_sequence_record,
                                                  alignments_bundle=_alignments_bundle))

                    _simple_blast_report: SimpleBlastReport
                    for _simple_blast_report in _simple_blast_reports:
                        # transforming regional blast hits to global

                        # print("before")
                        # print(f"_seqid: {_simple_blast_report.seqid}")
                        # print(f"hit_from: {_simple_blast_report.hit_from}")
                        # print(f"hit_to: {_simple_blast_report.hit_to}")

                        _seqid = _simple_blast_report.seqid
                        # parsing
                        _seqid = _seqid.split(sep=" ", maxsplit=1)[0]  # just in case
                        _seqid, _range = _seqid.split(sep="|", maxsplit=1)
                        _regional_hit_from, _regional_hit_to = _range.split(sep=":", maxsplit=1)
                        # updating values
                        _simple_blast_report.seqid = _seqid
                        _simple_blast_report.hit_from += int(_regional_hit_from)
                        _simple_blast_report.hit_to += int(_regional_hit_from)

                        # print("after")
                        # print(f"_seqid: {_seqid}")
                        # print(f"hit_from: {_simple_blast_report.hit_from}")
                        # print(f"hit_to: {_simple_blast_report.hit_to}")
                        # _headline, _sequence = get_genome_seq(_genome_id, _seqid)
                        # print(f"{_simple_blast_report.hit_strand}")
                        # print(f"{_simple_blast_report.hit_sequence.replace('-', '')}")
                        # print(f"{_sequence[_simple_blast_report.hit_from-1:_simple_blast_report.hit_to]}")
                        # print()

                    # with open(f"{alignments_bundle_output_file}.txt", 'w') as _output_file_handler:
                    #     _output_file_handler.write(str(_out_alignment_bundle))

                    (AlignmentsBundle.from_simple_blast_reports(sequence_record=_sequence_record,
                                                                simple_blast_reports=_simple_blast_reports)
                     .write_table_to_file(alignments_bundle_output_file, _data_format))

# TODO use this to find gene region.
#  visualize in R (alpha or colour gradient of exon arrow should be score)

# TODO: create code to parse exons and genes to a handy table.

# TODO:


# for taxonomy_id, gene_ids in taxonomy_to_gene_ids_dict.items():
#     print(f"processing taxonomy_id:[gene_ids]| {taxonomy_id}:{gene_ids}")
#
#     # # # # # # # downloading genomes # # # # # # #
#     scrape_genome.display_summery(taxonomy_id)  # Download species genomes (ignores if file already exist)
#     data_directory = scrape_genome.download_genome(taxonomy_id, include=['genome', 'gff3'])
#     #                   ./genomes/<taxonomy_id>
#     data_directory = f"{data_directory}/ncbi_dataset/data"
#     print(f"genome data for {taxonomy_id} in {data_directory}")
#
#     # Could have multiple genomes downloaded (if taxonomy id is not specific)
#     for _genome_id in os.listdir(data_directory):
#         genome_directory = f"{data_directory}/{_genome_id}"
#
#         # ignoring also files like:
#         # .../<data_directory>/assembly_data_report.jsonl
#         # .../<data_directory>/dataset_catalog.json
#         if os.path.isdir(genome_directory):
#
#             fasta_file: str = None  # ignore
#             gff_file: str = None  # ignore
#             for file in os.listdir(genome_directory):  # iterating trough fasta file(s) in .../genome_directory/*
#                 file = f"{genome_directory}/{file}"
#                 # print(file)
#                 if os.path.isfile(file) and (fnmatch(file, "*.fna") or fnmatch(file, "*.fasta")):
#                     fasta_file = file
#                 if os.path.isfile(file) and fnmatch(file, "*.gff"):
#                     gff_file = file
#
#             if fasta_file is None:
#                 raise ValueError(f"fasta file not found in {genome_directory}")
#             else:
#                 print(f"for {_genome_id} found {fasta_file}")
#
#                 # # # # # # # building blast db # # # # # # #
#                 create_blast_db(_genome_id, fasta_file, taxonomy_id)
#
#             # # # # # # # parsing genome_sequence_records # # # # # # #
#             if gff_file is None:
#                 sys.stderr.write(f"gff file not found in {genome_directory}")
#                 sys.stderr.write("skipping parsing of gff file")
#             else:
#                 print(f"parsing {gff_file} with gene_ids: {gene_ids}")  # optional?
#                 for _gene_id in gene_ids:
#                     print(f"parsing for {_gene_id}")
#                     genome_sequence_records = (
#                         os.path.basename(_genome_id),
#                         _gene_id,
#                         SequenceRecordsBundle.of_genes(
#                             genome_id=os.path.basename(_genome_id),
#                             gff_file=gff_file,
#                             genome_fasta_file=fasta_file,
#                             filter_types=["gene", "exon"],
#                             gene_id=_gene_id)
#                     )
#                     __per_genome_sequence_records_bundle.append(genome_sequence_records)


# taxonomy_to_gene_ids_dict: dict[int,list[int]] = {
#     7227: [  # Drosophila melanogaster (fruit fly)
#         41615  # timout  https://www.ncbi.nlm.nih.gov/gene/41615
#     ],
#
#     7460: [  # Apis mellifera (honey bee)   or (Apis mellifera Linnaeus, 1758)?
#         550811  # timout  https://www.ncbi.nlm.nih.gov/gene/550811
#
#     ],
#
#     7070: [  # Tribolium castaneum (red flour beetle)
#         659838  # timeout  https://www.ncbi.nlm.nih.gov/gene/659838
#     ],
#
#     454923: [  # Diachasma alloeum (wasps, ants & bees)
#         107037605  # timeout  https://www.ncbi.nlm.nih.gov/gene/107037605
#     ],
#
#     7425: [  # Nasonia vitripennis (jewel wasp)
#         100121375  # timeout homolog https://www.ncbi.nlm.nih.gov/gene/100121375
#     ],
#
# }
# using gene_id as they are easier for handling. NCBI prefers those.
