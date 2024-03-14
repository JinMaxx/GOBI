#!/usr/bin/python3
import os

import requests
import time
import Bio.SeqIO.FastaIO as FastaIO

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

#   DBT: Doubletime https://www.uniprot.org/uniprotkb/O76324/entry

# target species
taxonomy_ids = [50557]
#              true insects

# Uniprot API
# for reference: https://www.uniprot.org/help/api_queries
url = "https://uniprot.org/uniprotkb/search"
params = {
    '?format': 'fasta',
    'compressed': 'false'
}
# because requests is fucking stupid and removes the ? from urls as it thinks that the params might be empty...
# ? cant be at the end of the url because it gets parsed out.
# url encoding does not work.
# using a prepared url does not work.
# So instead using dirty fix by putting it in front of the first param.

directory = "./proteins"
files_list = os.listdir(directory)


def main():

    query = ""
    for taxonomy_id in taxonomy_ids:
        query = " AND ".join(f"(taxonomy_id: {taxonomy_id})")

    for gene_name in gene_names:

        filename = f"{gene_name.lower()}.fasta"

        if filename not in files_list:
            params['query'] = f"{query} AND (gene:{gene_name})"
            with requests.get(url, params=params) as response:  # stream=True
                print(response.url)  # for testing if url is correctly encoded.
                print(response.status_code)
                response.raise_for_status()
                # print(response.text)
                filename = f"{directory}/{filename}"
                print("file: ", filename)

                with open(filename, 'w') as file:
                    file.write(response.text)

                with open(filename, 'r') as file:
                    for (headline, sequence) in FastaIO.SimpleFastaParser(file):
                        print(headline)

            time.sleep(5)  # rate limiting


if __name__ == '__main__':
    main()
