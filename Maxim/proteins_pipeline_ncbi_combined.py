#!/usr/bin/python3
import os.path

from Bio import SeqIO
from Bio.SeqIO import FastaIO

import scrape_genes
import multiple_sequence_alignment

# this script combines several scripts as a pipeline for easier handling
# downloads proteins from ncbi, align them with mafft, display with them with AliView

# target species
gene_ids_dict: dict[str, list[int]] = {
    "timeout_and_timeless": [
        111642155,  # Centruroides sculpturatus
        107438930,  # Parasteatoda tepidariorum
        101740051,  # Bombyx mori
        550811,     # Apis mellifera
        41615,      # Drosophila melanogaster
        116766764,  # Danaus plexippus
        113500526,  # Trichoplusia ni
        105205832,  # Solenopsis invicta
        8038579,    # Ixodes scapularis
        108567601,  # Nicrophorus vespilloides 1
        108558671,  # Nicrophorus vespilloides 2
        659838,     # Tribolium castaneum
        5668064,    # Anopheles gambiae
        100121375,  # Nasonia vitripennis
        123275409,  # Cotesia glomerata
        103574560,  # Microplitis demolitor
        None,  # Centruroides sculpturatus
        107440384,  # Parasteatoda tepidariorum
        733065,  # Bombyx mori
        None,  # Apis mellifera
        33571,  # Drosophila melanogaster
        116765334,  # Danaus plexippus
        113501470,  # Trichoplusia ni
        None,  # Solenopsis invicta
        120838164,  # Ixodes scapularis
        108557741,  # Nicrophorus vespilloides
        657878,  # Tribolium castaneum
        None,  # Ixodes scapularis
        5668288,  # Anopheles gambiae
        None,  # Nasonia vitripennis
        None,  # Cotesia glomerata
        None,  # Microplitis demolitor
    ]
}


# no need to change those
genes_dir = "./genes"
aligned_proteins_dir = "./proteins_aligned_ncbi"


# just execute the script, easy peasy :P
if __name__ == '__main__':

    subdir: str
    gene_ids_list: list[int]
    for subdir, gene_ids_list in gene_ids_dict.items():

        download_directory: str = f"{genes_dir}/{subdir}"
        combined_fasta: str = f"{genes_dir}/{subdir}.fasta"
        protein_fasta_list: list[str] = list()

        gene_id: int
        for gene_id in gene_ids_list:

            # downloads all proteins from uniprot to ./proteins
            if gene_id is not None:
                print(f"scraping {gene_id} to {download_directory}")
                prefix_directory = scrape_genes.download_gene(gene_id=gene_id, directory=download_directory)
                protein_fasta = f"{prefix_directory}/ncbi_dataset/data/protein.faa"
                if os.path.isfile(protein_fasta): protein_fasta_list.append(protein_fasta)
                else: raise ValueError(f"{protein_fasta} is not a file")

        # from pathlib import Path
        # files = list(Path(".").rglob("*.[tT][xX][tT]"))

        # reading all fasta sequences to a single file.
        print("combining all fasta files")
        with open(combined_fasta, 'w') as combined_fasta_handle:
            writer = FastaIO.FastaWriter(combined_fasta_handle)

            for protein_fasta in protein_fasta_list:
                with open(protein_fasta, 'r') as protein_fasta_handle:
                    for record in SeqIO.parse(protein_fasta_handle, "fasta"):
                        writer.write_record(record)

    # generate alignments (requires mafft to be installed, see dependencies,sh)
    print(f"aligning proteins from {genes_dir} to {aligned_proteins_dir}")
    multiple_sequence_alignment.generate_alignments(input_directory=genes_dir,
                                                    output_directory=aligned_proteins_dir)

    # displays alignments (requires AliView to be installed, see dependencies.sh)
    print(f"displaying aligned proteins in {aligned_proteins_dir}")     # will open a lot of windows!
    multiple_sequence_alignment.display_alignments(alignments_directory=aligned_proteins_dir)
