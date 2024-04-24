#!/usr/bin/python3
import scrape_genome
from sequence_record_parser import SequenceRecordsBundle

# get exons and genes
# print them as csv like in blast exons
# no alignments

# gene_name -> [ (taxid, gene_id, assembly_level) ]
gene_dict: dict[str, list[(int, int, str)]] = {
    "timeout": [
        (218467, 111642155, None),           # Centruroides sculpturatus
        (114398, 107438930, None),           # Parasteatoda tepidariorum
        (7091,   101740051, "chromosome"),   # Bombyx mori
        (7460,   550811, "chromosome"),      # Apis mellifera
        (7227,   41615, "chromosome"),       # Drosophila melanogaster
        (13037,  116766764, "chromosome"),   # Danaus plexippus
        (7111,   113500526, "chromosome"),   # Trichoplusia ni
        (13686,  105205832, "chromosome"),   # Solenopsis invicta
        (6945,   8038579, None),             # Ixodes scapularis
        (110193, 108567601, None),           # Nicrophorus vespilloides 1
        (110193, 108558671, None),           # Nicrophorus vespilloides 2
        (7070,   659838, "chromosome"),      # Tribolium castaneum
        (7165,   5668064, "chromosome"),     # Anopheles gambiae
        (7425,   100121375, "chromosome"),   # Nasonia vitripennis
        (32391,  123275409, "chromosome"),   # Cotesia glomerata
        (69319,  103574560, "chromosome"),   # Microplitis demolitor
    ],

    "timeless": [
        (218467, None, None),                # Centruroides sculpturatus
        (114398, 107440384, None),           # Parasteatoda tepidariorum
        (7091,   733065, "chromosome"),      # Bombyx mori
        (7460,   None, "chromosome"),        # Apis mellifera
        (7227,   33571, "chromosome"),       # Drosophila melanogaster
        (13037,  116765334, "chromosome"),   # Danaus plexippus
        (7111,   113501470, "chromosome"),   # Trichoplusia ni
        (13686,  None, "chromosome"),        # Solenopsis invicta
        (6945,   120838164, None),           # Ixodes scapularis
        (110193, 108557741, None),           # Nicrophorus vespilloides
        (7070,   657878, "chromosome"),      # Tribolium castaneum
        (7165,   5668288, "chromosome"),     # Anopheles gambiae
        (7425,   None, "chromosome"),        # Nasonia vitripennis
        (32391,  None, "chromosome"),        # Cotesia glomerata
        (69319,  None, "chromosome"),        # Microplitis demolitor
    ],
}

genomes_dir: str = "./genomes"
output_file: str = "./get_exons_and_neighbours.csv"  # "./get_exons_and_genes.csv"
region_width: int = 100_000

data_format: str = (
    "gene_name,"
    "sr_taxonomy_id,sr_genome_id,sr_sequence,"
    # "sr_genome_sequence,"
    "gff_seqid,gff_source,gff_type,gff_start,gff_end,gff_score,gff_strand,gff_phase,gff_gene_id,"
    "gff_attribute:ID,gff_attribute:Name,gff_attribute:description,"
    "sr_headline"
    "\n"
)

debug = True

if __name__ == '__main__':

    with open(output_file, 'w') as output_file_handler:
        output_file_handler.write(data_format)  # overwriting with header

        gene_name: str
        tax_gene_id_list: list[int, int, str]

        for gene_name, tax_gene_id_list in gene_dict.items():

            if debug: print(f"processing {gene_name}")

            tax_id: int
            gene_id: int
            assembly_level: str
            for tax_id, gene_id, assembly_level in tax_gene_id_list:

                if debug: print(f"taxonomy_id: {tax_id}, gene_id: {gene_id}")

                if gene_id is None: continue

                for genome_id in \
                        scrape_genome.download_genome(taxonomy_id=tax_id, include=["genome", "gff3"],
                                                      genomes_dir=genomes_dir, assembly_level=assembly_level):

                    if debug: print(f"genome_id: {genome_id}")

                    fasta_file: str
                    gff_file: str
                    fasta_file, gff_file = scrape_genome.get_genome_files(taxonomy_id=tax_id, genome_id=genome_id,
                                                                          directory=genomes_dir)

                    if fasta_file is None or gff_file is None:
                        raise ValueError(f"fasta {fasta_file} or gff {gff_file} is None")
                    elif debug: print(f"parsing fasta {fasta_file} and gff {gff_file}")

                    sequence_records_bundle_list: list[SequenceRecordsBundle] = SequenceRecordsBundle.of_genes(
                        taxonomy_id=tax_id, genome_id=genome_id, gene_id=gene_id,
                        genome_fasta_file=fasta_file, gff_file=gff_file,
                        filter_types=["gene"]  # ["gene", "exon"]
                    )
                    for sequence_records_bundle in sequence_records_bundle_list:
                        for sequence_record in sequence_records_bundle.sequence_records:
                            line = sequence_record.as_table_record(data_format).replace("gene_name", gene_name)
                            if debug: print(f"\t{line}")
                            output_file_handler.write(line)

                    sequence_records_bundle_list = SequenceRecordsBundle.of_gene_neighbours(
                        taxonomy_id=tax_id, genome_id=genome_id, gene_id=gene_id,
                        genome_fasta_file=fasta_file, gff_file=gff_file,
                        filter_types=["gene"], region_width=region_width
                    )
                    for sequence_records_bundle in sequence_records_bundle_list:
                        for sequence_record in sequence_records_bundle.sequence_records:
                            line = (sequence_record
                                    .as_table_record(data_format)
                                    .replace("gene_name", f"neighbour_{gene_name}"))
                            if debug: print(f"\t{line}")
                            output_file_handler.write(line)
