#!/usr/bin/python3

import os
import sys
import scrape_genome

from fnmatch import fnmatch
from collections import defaultdict
from tempfile import NamedTemporaryFile, TemporaryDirectory
from blast_wrapper import SimpleBlastReport, blast, create_blast_db
from sequence_record_parser import SequenceRecordsBundle, SequenceRecord

import Bio.SeqIO.FastaIO as FastaIO


# better to use taxonomy id that are specific for species
taxonomy_to_gene_ids_dict = {
    7227: [41615]
}
# 7227 Drosophilia M. (for testing)
# using gene_id as they are easier for handling. NCBI prefers those.

__region_width: int = 300  # maybe higher 10.000

blastdb_directory = "./blast_db"
genomes_directory = "./genomes"
makeblastdb = "makeblastdb"
blastn = "blastn"
alignments_bundle_output_file = "./alignments_bundle_output_file"

# Some global variables for easier handling.

# genome_id -> seqid -> (headline, sequence)    saving already read sequence for given seqid per genome_id
genome_id_seqid_seqence_dict: dict[str, dict[str, (str, str)]] = defaultdict(dict)

# [(genome_id, gene_id, [GeneRecordsBundle])]
per_genome_sequence_records_bundle: list[(str, int, SequenceRecordsBundle)] = list()


for taxonomy_id, gene_ids in taxonomy_to_gene_ids_dict.items():
    print(f"processing taxonomy_id:[gene_ids]| {taxonomy_id}:{gene_ids}")

    # # # # # # # downloading genomes # # # # # # #
    scrape_genome.display_summery(taxonomy_id)  # Download species genomes (ignores if file already exist)
    data_directory = scrape_genome.download_genome(taxonomy_id, include=['genome', 'gff3'])
    #                   ./genomes/<taxonomy_id>
    data_directory = f"{data_directory}/ncbi_dataset/data"
    print(f"genome data for {taxonomy_id} in {data_directory}")

    # Could have multiple genomes downloaded (if taxonomy id is not specific)
    for genome_id in os.listdir(data_directory):
        genome_directory = f"{data_directory}/{genome_id}"

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

            if fasta_file is None:
                raise ValueError(f"fasta file not found in {genome_directory}")
            else:
                print(f"for {genome_id} found {fasta_file}")

                # # # # # # # building blast db # # # # # # #
                create_blast_db(genome_id, fasta_file, taxonomy_id)

            # # # # # # # parsing genome_sequence_records # # # # # # #
            if gff_file is None:
                sys.stderr.write(f"gff file not found in {genome_directory}")
                sys.stderr.write("skipping parsing of gff file")
            else:
                print(f"parsing {gff_file} with gene_ids: {gene_ids}")  # optional?
                for gene_id in gene_ids:
                    print(f"parsing for {gene_id}")
                    genome_sequence_records = (
                        os.path.basename(genome_id),
                        gene_id,
                        SequenceRecordsBundle.of_genes(
                            gff_file=gff_file,
                            genome_fasta_file=fasta_file,
                            filter_types=["gene", "exon"],
                            gene_id=gene_id)
                    )
                    per_genome_sequence_records_bundle.append(genome_sequence_records)


class Alignment:
    #         0      - region_width   align_from   hit_from      hit_to       align_to    + region_width      n
    # refer:  |-~~~--^----------------^------------|-------------|------------^-----------^---------------~~~-|
    # blast alignment:                             |||||||||||||||
    # query:                          |------------|-------------|------------|
    #                                 0            query_from    query_to     sequence_len
    #   can be negative!   align_from = hit_from - query_from        align_to = hit_to + (sequence_len - query_to)
    #   region_from = align_from - region_width >= 0                           region_to = align_to + region_width <= n

    def __init__(self,
                 align_from: int,
                 align_to: int,
                 simple_blast_report: SimpleBlastReport,
                 query_sequence_record: SequenceRecord):

        self.align_from = align_from
        self.align_to = align_to
        self.simple_blast_report = simple_blast_report
        self.sequence_record = query_sequence_record

    # I need to parse the fasta file here with the genomic_id
    @staticmethod
    def from_simple_blast_report(query_sequence_record: SequenceRecord,
                                 simple_blast_report: SimpleBlastReport) -> 'Alignment':

        print(query_sequence_record.sequence)
        print(f"report:\n{simple_blast_report}")

        hit_from = simple_blast_report.hit_from
        hit_to = simple_blast_report.hit_to
        query_from = simple_blast_report.query_from
        query_to = simple_blast_report.query_to

        # location flipped based on strand if - not +
        if hit_from > hit_to: hit_from, hit_to = hit_to, hit_from
        if query_from > query_to: query_from, query_to = hit_to, hit_from

        hit_gaps = simple_blast_report.hit_sequence.count('-')
        query_gaps = simple_blast_report.query_sequence.count('-')

        # calculating alignment position on genome
        align_from = hit_from - query_from
        print(f"align_from: {align_from}")
        align_to = hit_to + (len(query_sequence_record.sequence) - query_to) - hit_gaps
        print(f"align_to: {align_to}")

        # part1 = simple_blast_report.hit_from - align_from
        # print(f"hit_from -  align_from  = {simple_blast_report.hit_from} \t- {align_from}\t= {part1} = part1")
        # part2 = simple_blast_report.hit_to - simple_blast_report.hit_from
        # print(f"hit_to   -  hit_from    = {simple_blast_report.hit_to} \t- {simple_blast_report.hit_from}\t= {part2} = part2")
        # part3 = align_to - simple_blast_report.hit_to
        # print(f"align_to -  hit_to      = {align_to} \t- {simple_blast_report.hit_to}\t= {part3} = part3")
        # print(f"part1 + part2 + part3   = {part1} + {part2} + {part3} = {part1 + part2 + part3} == {len(query_sequence_record.sequence)} = query_sequence_len")

        # (align_to - align_from) - query_sequence_len
        # 76229 - 75225 = 1004         880 -        371 509 23
        # 76084 - 75225 = 859          430 -          0 430 15
        # 75755 - 75225 = 530          880 -        608 272 23
        # 75237 - 75225 = 12           708 -       1021 -313 15
        # 75237 - 75225 = 12           705 -        988 -283 15
        # 75501 - 75225 = 276          138 -          0 138 2
        # 75497 - 75225 = 272          136 -          0 136 2

        print(f"query_sequence_len = {len(query_sequence_record.sequence)}")
        print(f"align_to - align_from = {align_to} - {align_from} = {align_to - align_from}")
        print(f"align_len: {simple_blast_report.align_len}")
        print(f"gaps: {simple_blast_report.gaps}")

        # TODO: Mind the gaps!!
        #  reference is smaller then hit_sequence!
        #  input sequence is smaller then query_sequence!

        # TODO!! some chromosomes are split by multiples sequences with seqids.
        #  need to think about that too when using regions as they could span over multiple sequences

        return Alignment(align_from=align_from,
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
        if region_to >= ref_genome_seq_len-1: region_to = ref_genome_seq_len-1  # TODO check if last index is correct

        return region_from, region_to

    def as_region_fasta_entry(self, region_width: int) -> (str, str):
        headline: str
        sequence: str
        headline, sequence = get_genome_seq(genome_id=self.simple_blast_report.genome_id, seqid=self.simple_blast_report.seqid)
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

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        simple_blast_report_str = str(self.simple_blast_report).replace('\n', '\n\t')
        gff_record_str = str(self.sequence_record.gff_record).replace('\n', '\n\t')
        return (
            f"alignment: {self.align_from}-{self.align_from}\n"
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
            if alignment_bundle.genome_alignments_dict[simple_blast_report.genome_id].get(simple_blast_report.seqid) is None:
                alignment_bundle.genome_alignments_dict[simple_blast_report.genome_id][simple_blast_report.seqid] = list()

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
        for genome_id, alignments_dict in alignments_bundle.genome_alignments_dict.items():

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


def get_genome_seq(genome_id: str,
                   seqid: str) -> (str, str):
    
    if genome_id_seqid_seqence_dict[genome_id].get(seqid) is not None:
        return genome_id_seqid_seqence_dict[genome_id][seqid]
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
                            genome_id_seqid_seqence_dict[genome_id][seqid] = (headline, genome_sequence)
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
        f"gene_id: {gene_id}, type: {sequence_record.gff_record.type}, "
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

    alignments_bundle_output: AlignmentsBundle = AlignmentsBundle()

    genome_id: str
    gene_id: int
    sequence_records_bundle_list: list[SequenceRecordsBundle]
    for genome_id, gene_id, sequence_records_bundle_list in per_genome_sequence_records_bundle:
        # per_genome_gene_records_bundle: when there are many different reference genomes

        sequence_records_bundle: SequenceRecordsBundle
        for sequence_records_bundle in sequence_records_bundle_list:
            # for 1 gene there is 1 sequence_record_bundle with many sequences for exons and one sequence for the gene

            alignments_bundle: AlignmentsBundle = None  # ignore

            sequence_record: SequenceRecord
            for sequence_record in sequence_records_bundle.sequence_records:
                # either gene or exon

                if sequence_record.gff_record.type == "gene":
                    # finding the region where the gene exists
                    alignments_bundle = (AlignmentsBundle
                                         .from_simple_blast_reports(sequence_record=sequence_record,
                                                                    simple_blast_reports=gene_alignment(sequence_record)))
                    break
            if alignments_bundle is None: raise ValueError("No gene annotation found in sequence_records")
            alignments_bundle_output.append(alignments_bundle)

            # TESTING
            # for aligma_dict in alignments_bundle.genome_alignments_dict.values():
            #     for aligmas in aligma_dict.values():
            #         for aligma in aligmas:
            #             # hehe ligma balls
            #             seq1 = get_genome_seq(aligma.simple_blast_report.genome_id, aligma.simple_blast_report.seqid)[1][aligma.align_from:aligma.align_to]
            #             seq2 = f"{'-' * (aligma.simple_blast_report.hit_from - aligma.align_from)} {aligma.simple_blast_report.query_sequence}"
            #             print(f"seq1: {seq1}")
            #             print(f"seq2: {seq2}")
            #             # align_from   hit_from      hit_to       align_to
            # exit(0)

            for sequence_record in sequence_records_bundle.sequence_records:

                if sequence_record.gff_record.type == "exon":
                    simple_blast_reports: list[SimpleBlastReport] = (
                        exon_alignment_to_regions(genome_id=genome_id, 
                                                  region_width=__region_width, 
                                                  sequence_record=sequence_record, 
                                                  alignments_bundle=alignments_bundle))
                    simple_blast_report: SimpleBlastReport
                    for simple_blast_report in simple_blast_reports:
                        # transforming regional blast hits to global 
                        seqid = simple_blast_report.seqid
                        # print(f"parsing: {seqid}")
                        # parsing
                        seqid = seqid.split(sep=" ", maxsplit=1)[0]  # just in case
                        seqid, range = seqid.split(sep="|", maxsplit=1)
                        regional_hit_from, regional_hit_to = range.split(sep=":", maxsplit=1)
                        # updating values
                        simple_blast_report.seqid = seqid
                        simple_blast_report.hit_from -= int(regional_hit_from)
                        simple_blast_report.hit_to -= int(regional_hit_to)
                    alignments_bundle_output.append(AlignmentsBundle
                                                    .from_simple_blast_reports(sequence_record=sequence_record,
                                                                               simple_blast_reports=simple_blast_reports))

    with open(alignments_bundle_output_file, 'w') as output_file_handler:
        output_file_handler.write(str(alignments_bundle_output))


# TODO:
#  * make temporary directory for fasta_files
#  * make temporary fasta_files: genome_id.fasta
#   entries:
#    >seqid|region_from-region_to
#    <genomic_region_sequence>
#  * make temporary directory for blast_db
#  * create db
#  * blast against that db
#  * db_output should be parsed ->


# TODO use this to find gene region.
#  create table
#  visualize in R (alpha or colour gradient of exon arrow should be score)
#  create table like: taxonomy, start, end, exon, gene, reference_sequence, aligned_sequence
