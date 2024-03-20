#!/usr/bin/python3

import scrape_proteins
import filter_proteins
import multiple_sequence_alignment

# this script combines several scripts as a pipeline for easier handling
# downloads proteins from uniprot, remove sequences by keywords, align them with mafft, display with them with AliView

# target species
taxonomy_ids = [50557]  # can also be empty to include all species
#              true insects

# genes of interest
gene_names = ["timeless", "period", "clock", "cycle", "cwo", "vrille", "pdp1e", "shaggy", "cry1", "cry2", "timeout"]

# Gene name | Synonyms
# timeless  | tim, tim1, tim-1, timeless1
# period    | per
# clock     | clk
# cycle     |
# cwo       | clockwork orange
# vrille    | vri
# pdp1e     | Pdp1 DM1, DM32, mel_pdp1, pdp, PDP-1epsilon, Pdp1 epsilon, Pdp1epsilon, CG17888, Dmel_CG17888
# shaggy    | sgg, Dmel_CG2621, Dmsgg3, DMZ3K25Z GSK, Gsk-3, GSK-3B, Gsk-3beta, sgg-1, sgg-zw3, ZW-3
# cry1      |
# cry2      |
# timeout   | CG14381; CG14382; CG7855; CG8148; Dmel\CG7855; dtim2; tim; Tim-2; tim2; TIM2


filter_keywords = ["kaput", "fragment"]  # can also use id like "tr|A0A034V0Z2|A0A034V0Z2_BACDO"


# no need to change those
download_proteins_dir = "./proteins"
filtered_proteins_dir = "./proteins_filtered"
aligned_proteins_dir = "./proteins_aligned"


# just execute the script, easy peasy :P
if __name__ == '__main__':

    if True:  # downloads all proteins from uniprot to ./proteins
        print(f"#### scraping proteins ### \n{taxonomy_ids}, {gene_names}, {download_proteins_dir}")
        scrape_proteins.download_proteins(taxonomy_ids=taxonomy_ids,
                                          gene_names=gene_names,
                                          directory=download_proteins_dir)

    if True:  # filters proteins based on title (fasta) if includes the keywords
        print(f"#### filtering proteins ### \n{filter_keywords}, {download_proteins_dir}, {filtered_proteins_dir}")
        filter_proteins.filter_proteins(filter_keywords=filter_keywords,
                                        input_proteins_dir=download_proteins_dir,
                                        output_proteins_dir=filtered_proteins_dir)

    if True:  # generate alignments (requires mafft to be installed, see dependencies,sh)
        print(f"#### aligning filtered proteins ### \n{filtered_proteins_dir}, {aligned_proteins_dir}")
        multiple_sequence_alignment.generate_alignments(input_directory=filtered_proteins_dir,
                                                        output_directory=aligned_proteins_dir)

    if True:  # displays alignments (requires AliView to be installed, see dependencies.sh)
        print(f"#### displaying aligned proteins ### \n{aligned_proteins_dir}")     # will open a lot of windows!
        multiple_sequence_alignment.display_alignments(alignments_directory=aligned_proteins_dir)


# TODO:do this but with blasted genome trascribed translated and compare those proteins.
