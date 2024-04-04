#!/usr/bin/python3
from collections.abc import Generator

from alignment import AlignmentsBundle

import scrape_genome

from tempfile import NamedTemporaryFile, TemporaryDirectory
from ncbi_search_parser import read_all_to_sequence_records_bundle
from blast_wrapper import SimpleBlastReport, blast, create_blast_db
from sequence_record_parser import SequenceRecordsBundle, SequenceRecord


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
    "gff_attribute:GeneID,query_sequence,hit_sequence,mid_line\n"

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

# Some global variables for easier handling.
__region_width: int = 1_000  # maybe higher 10_000
blastdb_directory = "./blast_db"
genomes_directory = "./genomes"
makeblastdb = "makeblastdb"
blastn = "blastn"
alignments_bundle_output_file = "./alignments_bundle_output.csv"


# genome_id -> seqid -> (headline, sequence)    saving already read sequence for given seqid per genome_id
# __genome_id_seqid_seqence_dict: dict[str, dict[str, (str, str)]] = defaultdict(dict)
# THIS BUFFER COULD BLOW UP THE MEMORY

# [(genome_id, gene_id, [GeneRecordsBundle])]
# __per_genome_sequence_records_bundle: list[(str, int, SequenceRecordsBundle)] = list()


# building blast_db
for taxonomy_id in taxonomies:
    print(f"processing taxonomy_id: {taxonomy_id}")

    # # # # # # # downloading genomes # # # # # # #
    scrape_genome.display_summery(taxonomy_id)  # Download species genomes (ignores if file already exist)
    genome_ids: list[str] = scrape_genome.download_genome(taxonomy_id, include=['genome'])

    for _genome_id in genome_ids:
        fasta_file, _ = scrape_genome.get_genome_files(taxonomy_id, _genome_id)
        # # # # # # building blast db # # # # # #
        create_blast_db(_genome_id, fasta_file, taxonomy_id)


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
    return to_tmp_fasta(sequence=sequence_record.sequence, headline=(
        f"gene_id: {sequence_record.gff_record.get_attribute('GeneID')}, type: {sequence_record.gff_record.type}, "
        f"reference_seqid: {sequence_record.gff_record.seqid}, "
        f"start: {sequence_record.gff_record.start}, end: {sequence_record.gff_record.end}"))


def gene_alignment(sequence_record: SequenceRecord) -> Generator[SimpleBlastReport]:
    query_input_file = __sequence_record_to_tmp_fasta(sequence_record)
    yield from blast(query_fasta_file=query_input_file.name)


# I might have to create multiple temporary db_tables
def exon_alignment_to_regions(region_width: int,
                              sequence_record: SequenceRecord,  # exon which is going to be used as blast query
                              alignments_bundle: AlignmentsBundle  # multiple alignments of 1 gene to many genomes
                              ) -> list[SimpleBlastReport]:
    # using tmp files in fasta format instead of using stdin

    tmp_blast_db_dir: TemporaryDirectory = TemporaryDirectory()

    genome_id: str
    fasta_entries: list[(str, str)]
    for genome_id, fasta_entries in alignments_bundle.regions_as_fasta_entries(region_width=region_width):
        for headline, sequence in fasta_entries:
            create_blast_db(genome_id=genome_id,
                            fasta_file=to_tmp_fasta(headline=headline, sequence=sequence).name,
                            other_folder=tmp_blast_db_dir.name)

    yield from blast(query_fasta_file=__sequence_record_to_tmp_fasta(sequence_record).name,
                     other_blast_db_directory=tmp_blast_db_dir.name)


if __name__ == '__main__':

    with open(alignments_bundle_output_file, 'w') as _output_file_handler:
        _output_file_handler.write(_data_format)  # overwriting with header

    _sequence_records_bundle: SequenceRecordsBundle
    for _sequence_records_bundle in read_all_to_sequence_records_bundle(directory="./ncbi_search_conservative"):
        # parsed from ./ncbi_search_*/*.tsv
        # generates sequence_record_bundles per gene_id containing gene and exons sequences

        # Getting the AlignmentsBundle just for type == "gene"
        _alignments_bundle: AlignmentsBundle = None  # ignore
        _sequence_record: SequenceRecord
        for _sequence_record in _sequence_records_bundle.sequence_records:
            if _sequence_record.gff_record.type == "gene":
                _alignments_bundle = (AlignmentsBundle  # finding the region where the gene exists
                                      .from_simple_blast_reports(sequence_record=_sequence_record,
                                                                 simple_blast_reports=gene_alignment(_sequence_record)))
                break                                                               # ^ alignment is done here
        if _alignments_bundle is None: raise ValueError("No gene annotation found in sequence_records")
        else:
            for _sequence_record in _sequence_records_bundle.sequence_records:

                if _sequence_record.gff_record.type == "exon":
                    _simple_blast_reports: list[SimpleBlastReport] = list(
                        exon_alignment_to_regions(region_width=__region_width,
                                                  sequence_record=_sequence_record,
                                                  alignments_bundle=_alignments_bundle))

                    _simple_blast_report: SimpleBlastReport
                    for _simple_blast_report in _simple_blast_reports:  # TypeError: 'NoneType' object is not iterable
                        # transforming regional blast hits to global
                        _seqid = _simple_blast_report.seqid
                        # parsing
                        _seqid = _seqid.split(sep=" ", maxsplit=1)[0]  # just in case
                        _seqid, _range = _seqid.split(sep="|", maxsplit=1)
                        _regional_hit_from, _regional_hit_to = _range.split(sep=":", maxsplit=1)
                        # updating values
                        _simple_blast_report.seqid = _seqid
                        _simple_blast_report.hit_from += int(_regional_hit_from)
                        _simple_blast_report.hit_to += int(_regional_hit_from)

                    _alignments_bundle.append(
                        AlignmentsBundle.from_simple_blast_reports(sequence_record=_sequence_record,
                                                                   simple_blast_reports=_simple_blast_reports))

            _alignments_bundle.write_table_to_file(alignments_bundle_output_file, _data_format)

# TODO use this to find gene region.
#  visualize in R (alpha or colour gradient of exon arrow should be score)
#  create code to parse exons and genes to a handy table.
