#!/usr/bin/python3
import os
from pathlib import Path

from alignment import AlignmentsBundle

import scrape_genome

from tempfile import NamedTemporaryFile, TemporaryDirectory
from ncbi_search_parser import read_all_to_sequence_records_bundle
from blast_wrapper import SimpleBlastReport, blastn, tblastx, create_blast_db
from sequence_record_parser import SequenceRecordsBundle, SequenceRecord


# for the resulting .csv
_data_format = (
    # This is a format containing all parameters

    # "hit_from,hit_to,query_from,query_to,hit_gaps,query_gaps,align_from,align_to,"
    # "query_sequence,hit_sequence,mid_line,align_len,identity,query_strand,hit_strand,"
    # "bit_score,evalue,hit_seqid,hit_genome_id,query_genome_id,hit_taxonomy_id,query_taxonomy_id,full_query_sequence,"
    # "gff_seqid,gff_source,gff_type,gff_start,gff_end,gff_score,gff_strand,gff_phase,gff_attributes,"
    # "gff_attribute:parent,gff_gene_id\n"

    "hit_from,hit_to,query_from,query_to,hit_gaps,query_gaps,align_from,align_to,"
    "align_len,identity,query_strand,hit_strand,"
    "bit_score,evalue,hit_seqid,hit_genome_id,query_genome_id,hit_taxonomy_id,query_taxonomy_id,"
    "gff_seqid,gff_source,gff_type,gff_start,gff_end,gff_score,gff_strand,gff_phase,gff_gene_id,gff_attribute:ID"
    # "query_sequence,hit_sequence,mid_line"
    "\n"

    # maybe find a way to also include query_taxid, hit_taxid
)

# species to look against:
taxonomies: list[int] = [  # better to use taxonomy id that are specific for species
    # datasets:
    218467,  # Bark scorpion (Centruroides sculpturatus)
    114398,  # Common house spider (Parasteatoda tepidariorum)
    7091,    # Domestic silkworm (Bombyx mori)
    7460,    # Honey bee (Apis mellifera)
    13037,   # Monarch butterfly (Danaus plexippus)
    7227,    # Fruit fly (Drosophila melanogaster)
    7111,    # Cabbage looper (Trichoplusia ni)
    13686,   # Red fire ant (Solenopsis invicta)
    6945,    # Black-legged tick (Ixodes scapularis)
    7070,    # Red flour beetle (Tribolium castaneum)
    110193,  # (Nicrophorus vespilloides)

    # interesting insects to compare to:
    7119,    # Chinese oak silkmoth (Antheraea pernyi)
    189913,  # butterfly family Pieridae[white, yellow, sulphur] (Leptidea sinapis)

    # other insects:
    7165,    # African malaria mosquito (Anopheles gambiae)
    79782,   # Bed bug (Cimex lectularius)
]

# Some global variables for easier handling.
__region_width: int = 10_000  # maybe higher 10_000
blastdb_directory = "./blast_db"
genomes_directory = "./genomes"
makeblastdb = "makeblastdb"
blastn = "blastn"
alignments_bundle_output_file = "./alignments_bundle_output.csv"

debug = True


# genome_id -> seqid -> (headline, sequence)    saving already read sequence for given seqid per genome_id
# __genome_id_seqid_seqence_dict: dict[str, dict[str, (str, str)]] = defaultdict(dict)
# THIS BUFFER COULD BLOW UP THE MEMORY

# [(genome_id, gene_id, [GeneRecordsBundle])]
# __per_genome_sequence_records_bundle: list[(str, int, SequenceRecordsBundle)] = list()


# building blast_db
if debug: print("building blast db")
for _taxonomy_id in taxonomies:
    if debug: print(f"processing taxonomy_id: {_taxonomy_id}")

    # # # # # # # downloading genomes # # # # # # #
    scrape_genome.display_summery(_taxonomy_id)  # Download species genomes (ignores if file already exist)
    genome_ids: list[str] = scrape_genome.download_genome(_taxonomy_id, include=['genome', 'gff3'])
    if debug: print(f"taxonomy|[genome_ids]: {_taxonomy_id}|{genome_ids}")

    for _genome_id in genome_ids:
        fasta_file, _ = scrape_genome.get_genome_files(_taxonomy_id, _genome_id)
        # # # # # # building blast db # # # # # #
        create_blast_db(_genome_id, fasta_file, _taxonomy_id)


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


def __sequence_record_to_tmp_fasta(sequence_record: SequenceRecord) -> NamedTemporaryFile:
    return to_tmp_fasta(sequence=sequence_record.sequence, headline=(
        f"gene_id: {sequence_record.gff_record.get_attribute('GeneID')}, type: {sequence_record.gff_record.type}, "
        f"reference_seqid: {sequence_record.gff_record.seqid}, "
        f"start: {sequence_record.gff_record.start}, end: {sequence_record.gff_record.end}"))


def gene_alignment(sequence_records_bundle: SequenceRecordsBundle) -> AlignmentsBundle:
    # Getting the Alignment just for type "gene"

    _alignments_bundle_gene: AlignmentsBundle = AlignmentsBundle()  # ignore

    _sequence_record: SequenceRecord
    for _sequence_record in sequence_records_bundle.sequence_records:

        if _sequence_record.gff_record.type == "gene":

            if debug: print("gene alignment of taxonomy_id|genome_id|gene_id: "
                            f"{sequence_records_bundle.taxonomy_id}|"
                            f"{sequence_records_bundle.genome_id}|"
                            f"{_sequence_record.gff_record.get_attribute_dict_value('dbxref', 'GeneID')}")

            _simple_blast_reports: list[SimpleBlastReport]
            with __sequence_record_to_tmp_fasta(_sequence_record) as tmp_fasta_file:
                _simple_blast_reports = list(blastn(query_fasta_file=tmp_fasta_file.name,
                                                    negative_genome_ids=[_sequence_record.get_genome_id()]))

            if debug: print(f"number of blast results (genes) {len(_simple_blast_reports)}\n"
                            f"genomes hit: {set([_s_b_report.genome_id for _s_b_report in _simple_blast_reports])}\n")

            _alignments_bundle_gene.append(AlignmentsBundle  # finding the region where the gene exists
                                           .from_simple_blast_reports(sequence_record=_sequence_record,
                                                                      simple_blast_reports=_simple_blast_reports))
            # break  # ^ alignment is done here

    return _alignments_bundle_gene

    # if _alignments_bundle_gene.is_empty(): raise ValueError("No gene annotation found in sequence_records")
    # else: return _alignments_bundle_gene  # This seems to be erroneous if no alignments are found


def exon_alignment_to_regions(region_width: int,
                              sequence_records_bundle: SequenceRecordsBundle,  # exon which is going to be used as blast query
                              alignments_bundle: AlignmentsBundle  # multiple alignments of 1 gene to many genomes
                              ) -> AlignmentsBundle:  # it is the same from before but with combined with other
    # using tmp files in fasta format instead of using stdin

    _tmp_blast_db_dir: TemporaryDirectory = TemporaryDirectory()

    _taxonomy_id: int
    _genome_id: str
    _fasta_entries: (str, str)
    for _taxonomy_id, _genome_id, _fasta_entries in alignments_bundle.regions_as_fasta_entries(region_width):
        with NamedTemporaryFile() as _tmp_fasta_file:
            _headline: str
            _sequence: str
            for _headline, _sequence in _fasta_entries:
                to_tmp_fasta(headline=_headline, sequence=_sequence, tmp_fasta_file=_tmp_fasta_file)
            # TODO: check if the file is written as should be

            # # # # # Building temporary blast dbs # # # # #
            if debug: print(f"creating temporary blast db taxonomy_id|genome_id: {_taxonomy_id}|{_genome_id}")
            create_blast_db(taxonomy_id=_taxonomy_id,
                            genome_id=_genome_id,
                            fasta_file=_tmp_fasta_file.name,
                            other_folder=_tmp_blast_db_dir.name)

    if debug: print(f"created following temporary blast databases in {_tmp_blast_db_dir.name}: "
                    f"{set([Path(__db_file).stem for __db_file in os.listdir(_tmp_blast_db_dir.name)])}")

    # TODO somewhere here it looks for an gff file for genome. maybe by loading the fasta file.

    # _alignment: Alignment
    # for _alignment in alignments_bundle.alignments():
    #
    #     _taxonomy_id: int = _alignment.simple_blast_report.taxonomy_id
    #     _genome_id: str = _alignment.simple_blast_report.genome_id
    #     _headline: str
    #     _sequence: str
    #     _headline, _sequence = _alignment.as_region_fasta_entry(region_width)
    #
    #     # # # # # Building temporary blast dbs # # # # #
    #     # if debug: print(f"creating temporary blast db taxonomy_id|genome_id|hit_seq_id[span]: "  # to slow
    #     #                 f"{_taxonomy_id}|{_genome_id}|{_hit_seq_id}[{_alignment.get_span(region_width)}]")
    #     if debug: print(f"creating temporary blast db taxonomy_id|genome_id|(seq_id|from:to): "
    #                     f"{_taxonomy_id}|{_genome_id}|({_headline.split(' ', maxsplit=1)[0]})")
    #     with to_tmp_fasta(headline=_headline, sequence=_sequence) as _tmp_fasta_file:
    #         create_blast_db(taxonomy_id=_taxonomy_id,
    #                         genome_id=_genome_id,
    #                         fasta_file=_tmp_fasta_file.name,
    #                         other_folder=_tmp_blast_db_dir.name)
    #     # TODO I have to bundle them into a single per genome_id fasta file...

    _sequence_record: SequenceRecord
    for _sequence_record in sequence_records_bundle.sequence_records:
        if _sequence_record.gff_record.type == "exon":
            if debug: print("exon alignment of taxonomy_id|genome_id|gene_id|gff_id: "
                            f"{sequence_records_bundle.taxonomy_id}|"
                            f"{sequence_records_bundle.genome_id}|"
                            f"{_sequence_record.gff_record.get_attribute_dict_value('dbxref', 'GeneID')}|"
                            f"{_sequence_record.gff_record.get_attribute('ID')}")
            with __sequence_record_to_tmp_fasta(_sequence_record) as tmp_fasta_file:

                _simple_blast_reports: list[SimpleBlastReport] = list()
                _simple_blast_report: SimpleBlastReport
                for _simple_blast_report in tblastx(query_fasta_file=tmp_fasta_file.name,
                                                   other_blast_db_directory=_tmp_blast_db_dir.name):
                    # transforming regional blast hits to global
                    _seq_id = _simple_blast_report.seqid
                    # parsing
                    _seq_id = _seq_id.split(sep=" ", maxsplit=1)[0]  # just in case
                    _seq_id, _range = _seq_id.split(sep="|", maxsplit=1)
                    _regional_hit_from, _regional_hit_to = _range.split(sep=":", maxsplit=1)
                    # updating values
                    _simple_blast_report.seqid = _seq_id
                    _simple_blast_report.hit_from += int(_regional_hit_from)
                    _simple_blast_report.hit_to += int(_regional_hit_from)
                    if debug: print(f"blast_result exon taxonomy_id|genome_id|seq_id[span]: "
                                    f"{_simple_blast_report.taxonomy_id}|{_simple_blast_report.taxonomy_id}|{_seq_id}"
                                    f"[({_simple_blast_report.hit_from}:{_simple_blast_report.hit_to})]")
                    # to list
                    _simple_blast_reports.append(_simple_blast_report)

                if debug: print(f"number of blast results (exons): {len(_simple_blast_reports)}\n")
                alignments_bundle.append(
                    AlignmentsBundle.from_simple_blast_reports(sequence_record=_sequence_record,
                                                               simple_blast_reports=_simple_blast_reports))

    _tmp_blast_db_dir.cleanup()
    return alignments_bundle


if __name__ == '__main__':

    with open(alignments_bundle_output_file, 'w') as _output_file_handler:
        _output_file_handler.write(_data_format)  # overwriting with header

    _sequence_records_bundle: SequenceRecordsBundle
    for _sequence_records_bundle in read_all_to_sequence_records_bundle(directory="./ncbi_search_conservative"):

        # if _sequence_records_bundle.taxonomy_id in taxonomies: continue  # skip aligning with itself.

        # parsed from ./ncbi_search_*/*.tsvto_tmp_fasta
        # generates sequence_record_bundles per gene_id containing gene and exons sequences
        if debug: print(f"sequence_records_bundle:\n{_sequence_records_bundle}")

        _gene_alignments_bundle = gene_alignment(_sequence_records_bundle)
        # if debug: print(f"gene_alignments_bundle: \n{_gene_alignments_bundle}")

        _alignments_bundle = exon_alignment_to_regions(region_width=__region_width,
                                                       sequence_records_bundle=_sequence_records_bundle,
                                                       alignments_bundle=_gene_alignments_bundle)
        if debug: _alignments_bundle.print()

        _alignments_bundle.write_table_to_file(alignments_bundle_output_file, _data_format)

# TODO use this to find gene region.
#  visualize in R (alpha or colour gradient of exon arrow should be score)
#  create code to parse exons and genes to a handy table.
