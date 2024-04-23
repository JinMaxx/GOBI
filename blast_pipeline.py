#!/usr/bin/python3

from alignment import AlignmentsBundle

import scrape_genome

from tempfile import NamedTemporaryFile
from ncbi_search_parser import read_all_to_sequence_records_bundle, read_all_to_sequence_records_bundle_neighbours
from blast_wrapper import SimpleBlastReport, tblastx, create_blast_db
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
    "query_sequence,hit_sequence,mid_line"
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
__region_width: int = 100_000  # maybe higher 10_000
blastdb_directory = "./blast_db"
genomes_directory = "./genomes"
makeblastdb = "makeblastdb"
blastn = "blastn"
alignments_bundle_output_file = "./aligned_exons.csv"
ncbi_search_directory = "./ncbi_search_continued"  # "./ncbi_search_conservative"

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


def alignment(sequence_records_bundle: SequenceRecordsBundle) -> AlignmentsBundle:

    _alignments_bundle: AlignmentsBundle = AlignmentsBundle()

    _sequence_record: SequenceRecord
    for _sequence_record in sequence_records_bundle.sequence_records:

        if _sequence_record.gff_record.type == "gene":  # or just put exon here or whatever

            if debug: print("alignment of taxonomy_id|genome_id|gene_id: "
                            f"{sequence_records_bundle.taxonomy_id}|"
                            f"{sequence_records_bundle.genome_id}|"
                            f"{_sequence_record.gff_record.get_attribute_dict_value('dbxref', 'GeneID')}")

            _simple_blast_reports: list[SimpleBlastReport]
            with __sequence_record_to_tmp_fasta(_sequence_record) as tmp_fasta_file:
                _simple_blast_reports = list(tblastx(query_fasta_file=tmp_fasta_file.name,
                                                     negative_genome_ids=[_sequence_record.get_genome_id()]))

            if debug: print(f"number of blast results {len(_simple_blast_reports)}\n"
                            f"genomes hit: {set([_s_b_report.genome_id for _s_b_report in _simple_blast_reports])}\n")

            _alignments_bundle.append(AlignmentsBundle  # finding the region where the gene exists
                                      .from_simple_blast_reports(sequence_record=_sequence_record,
                                                                 simple_blast_reports=_simple_blast_reports))
            # break  # ^ alignment is done here

    return _alignments_bundle

    # if _alignments_bundle_gene.is_empty(): raise ValueError("No gene annotation found in sequence_records")
    # else: return _alignments_bundle_gene  # This seems to be erroneous if no alignments are found


if __name__ == '__main__':

    with open(alignments_bundle_output_file, 'w') as _output_file_handler:
        _output_file_handler.write(_data_format)  # overwriting with header

    _sequence_records_bundle: SequenceRecordsBundle
    # for _sequence_records_bundle in read_all_to_sequence_records_bundle(directory=ncbi_search_directory): TODO CHANGE
    for _sequence_records_bundle in read_all_to_sequence_records_bundle_neighbours(directory=ncbi_search_directory, region_width=__region_width):

        # if _sequence_records_bundle.taxonomy_id in taxonomies: continue  # skip aligning with itself.

        # parsed from ./ncbi_search_*/*.tsv to_tmp_fasta
        if debug: print(f"sequence_records_bundle:\n{_sequence_records_bundle}")

        _alignments_bundle = alignment(_sequence_records_bundle)
        if debug: _alignments_bundle.print()

        _alignments_bundle.write_table_to_file(alignments_bundle_output_file, _data_format)

# TODO use this to find gene region.
#  visualize in R (alpha or colour gradient of exon arrow should be score)
#  create code to parse exons and genes to a handy table.
