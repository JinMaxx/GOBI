#!/usr/bin/python3

from collections.abc import Generator, Iterable
from typing import Self

import scrape_genome

from collections import defaultdict
from blast_wrapper import SimpleBlastReport
from sequence_record_parser import SequenceRecord


genomes_directory = "./genomes"


class Alignment:
    #         0      - region_width   align_from   hit_from      hit_to       align_to    + region_width      n
    # refer:  |-~~~--^----------------^------------|-------------|------------^-----------^---------------~~~-|
    # blast alignment:                             |||||||||||||||
    # query:                          |------------|-------------|------------|
    #                                 0            query_from    query_to     sequence_len
    #   can be negative!   align_from = hit_from - query_from        align_to = hit_to + (sequence_len - query_to)
    #   region_from = align_from - region_width >= 0                           region_to = align_to + region_width <= n

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

        return Alignment(hit_from=hit_from, hit_to=hit_to,
                         query_from=query_from, query_to=query_to,
                         hit_gaps=hit_gaps, query_gaps=query_gaps,
                         align_from=align_from,
                         align_to=align_to,
                         simple_blast_report=simple_blast_report,
                         query_sequence_record=query_sequence_record)

    def __get_hit_genome_headline_sequence(self) -> (str, str):
        return scrape_genome.get_genome_seq(taxonomy_id=self.simple_blast_report.taxonomy_id,
                                            genome_id=self.simple_blast_report.genome_id,
                                            seq_id=self.simple_blast_report.seqid)

    def get_span(self, region_width: int, ref_genome_seq_len: int = None) -> (int, int):

        if ref_genome_seq_len is None:
            _, ref_genome_seq = self.__get_hit_genome_headline_sequence()
            ref_genome_seq_len = len(ref_genome_seq)

        # applying region_width
        region_from = self.align_from - region_width
        region_to = self.align_to + region_width

        # checking for boundaries
        if region_from <= 0: region_from = 0
        if region_to >= ref_genome_seq_len - 1: region_to = ref_genome_seq_len - 1

        return region_from, region_to

    def __get_region_sequence(self, region_width: int) -> (int, int):
        _, sequence = self.__get_hit_genome_headline_sequence()
        region_from, region_to = self.get_span(region_width, len(sequence))
        return sequence[region_from:region_to]

    def as_region_fasta_entry(self, region_width: int) -> (str, str):
        headline: str
        sequence: str
        headline, sequence = self.__get_hit_genome_headline_sequence()
        region_from, region_to = self.get_span(region_width, len(sequence))
        if headline.startswith(self.simple_blast_report.seqid):
            second_split = headline.split(' ', maxsplit=1)[1]
            headline = f"{self.simple_blast_report.seqid}|{region_from}:{region_to}"
            if second_split is not None: headline = f"{headline} {second_split}"
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
        if "region_sequence" in data_format:
            data_format = data_format.replace("region_sequence", str(self.__get_region_sequence(1_000)))

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
            data_format = data_format.replace("query_genome_id", "sr_genome_id")
        if "query_taxonomy_id" in data_format:
            data_format = data_format.replace("query_taxonomy_id", "sr_taxonomy_id")
        if "full_query_sequence" in data_format:
            data_format = data_format.replace("full_query_sequence", "sr_sequence")
        if "query_headline" in data_format:
            data_format = data_format.replace("query_headline", "sr_headline")

        data_format = self.sequence_record.as_table_record(data_format)

        return data_format

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        simple_blast_report_str = str(self.simple_blast_report).replace('\n', '\n\t')
        gff_record_str = str(self.sequence_record.gff_record).replace('\n', '\n\t')
        return (
            f"alignment: [{self.align_from:10.0f}:{self.align_from:10.0f}]\n"
            f"query: [{self.query_from:10.0f}:{self.query_to:10.0f}], gaps: {self.query_gaps}\n"
            f"  hit: [{self.hit_from:10.0f}:{self.hit_to:10.0f}], gaps: {self.hit_gaps}\n"
            "simple_blast_report:\n"
            f"\t{simple_blast_report_str}\n"
            "sequence_record.gff_record:\n"
            f"\t{gff_record_str}"
        )


class AlignmentsBundle:  # remove overlaps and concatenate regions?

    def __init__(self):            # hit_genome_id -> hit_seqid -> Alignment
        self.genome_alignments_dict: dict[str, dict[str, list[Alignment]]] = defaultdict(lambda: defaultdict(list))
        # I somehow forgot to include taxonomy_id, so it is taxonomy_id -> genome_id -> seq_id -> alignment
        # but rewriting the code is too late at this point...

    @staticmethod
    def from_simple_blast_reports(sequence_record: SequenceRecord,
                                  simple_blast_reports: Iterable[SimpleBlastReport]) -> 'AlignmentsBundle':

        alignment_bundle: AlignmentsBundle = AlignmentsBundle()

        for simple_blast_report in simple_blast_reports:

            (alignment_bundle.genome_alignments_dict[simple_blast_report.genome_id][simple_blast_report.seqid].
             append(Alignment.from_simple_blast_report(query_sequence_record=sequence_record,
                                                       simple_blast_report=simple_blast_report)))

        return alignment_bundle

    #                                                  # [(taxonomy_id, genome_id, (headline, sequence))]
    def regions_as_fasta_entries(self, region_width: int) -> Generator[(str, str, list[(str, str)])]:
        taxonomy_genome_fasta_entries_dict: dict[int, dict[str, list[(str, str)]]] = defaultdict(lambda: defaultdict(list))

        # because I forgot taxonomy as mentioned earlier I will lose some performance by building this dict here.
        alignment: Alignment
        for alignment in self.alignments():
            taxonomy_id: int = alignment.simple_blast_report.taxonomy_id
            genome_id: str = alignment.simple_blast_report.genome_id
            taxonomy_genome_fasta_entries_dict[taxonomy_id][genome_id].append(alignment.as_region_fasta_entry(region_width))

        taxonomy_id: int
        genome_fasta_entries_dict: dict[str, list[(str, str)]]
        for taxonomy_id, genome_fasta_entries_dict in taxonomy_genome_fasta_entries_dict.items():
            genome_id: str
            fasta_entries: list[(str, str)]
            for genome_id, fasta_entries in genome_fasta_entries_dict.items():
                yield taxonomy_id, genome_id, fasta_entries

    def alignments(self) -> Generator[Alignment]:
        alignments_dict: dict[str, list[Alignment]]
        for alignments_dict in self.genome_alignments_dict.values():
            alignments: list[Alignment]
            for alignments in alignments_dict.values():
                alignment: Alignment
                for alignment in alignments:
                    yield alignment

    def append(self, alignments_bundle: 'AlignmentsBundle') -> Self:
        alignments_seqid_dict: dict[str, list[Alignment]]
        for genome_id, alignments_seqid_dict in alignments_bundle.genome_alignments_dict.items():
            alignments: list[Alignment]
            for seqid, alignments in alignments_seqid_dict.items():
                self.genome_alignments_dict[genome_id][seqid].extend(alignments)
        return self

    def is_empty(self):
        return len(self.genome_alignments_dict) <= 0

    # def __repr__(self):
    #     return self.__str__()
    #
    # def __str__(self):  # this can get huge
    #
    #     output_string: str = str()
    #
    #     alignments_dict: dict[str, list[Alignment]]
    #     for genome_id, alignments_dict in self.genome_alignments_dict.items():
    #
    #         output_string = (
    #             f"{output_string}\n"
    #             f"{genome_id}:\n"
    #         )
    #
    #         alignments: list[Alignment]
    #         for seqid, alignments in alignments_dict.items():
    #
    #             output_string = (
    #                 f"{output_string}\n"
    #                 f"\n\t{seqid}:\n\n"
    #             )
    #
    #             alignment: Alignment
    #             for alignment in alignments:
    #                 alignment_str = str(alignment).replace('\n', '\n\t\t')
    #                 output_string = (
    #                     f"{output_string}"
    #                     f"\t\t{alignment_str}\n\n"
    #                 )
    #
    #     return output_string

    def print(self):

        alignments_dict: dict[str, list[Alignment]]
        for genome_id, alignments_dict in self.genome_alignments_dict.items():
            print(f"{genome_id}:")

            alignments: list[Alignment]
            for seqid, alignments in alignments_dict.items():
                print(f"\n\t{seqid}:\n")

                alignment: Alignment
                for alignment in alignments:
                    alignment_str = str(alignment).replace('\n', '\n\t\t')
                    print(f"\t\t{alignment_str}\n")

    def write_table_to_file(self, file: str, format: str):
        with open(file, 'a') as file_handler:
            alignment: Alignment
            for alignment in self.alignments():
                file_handler.write(alignment.as_table_record(data_format=format))