#!/usr/bin/python3
import os

import requests
import time
import Bio.SeqIO.FastaIO as FastaIO


# target species
_taxonomy_ids = [50557]
#              true insects

# genes of interest
_gene_names = ["timeless", "period", "clock", "cycle", "cwo", "vrille", "pdp1e", "shaggy", "cry1", "cry2", "timeout"]

# TODO: update synonyms
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


# No need to change the following variables

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

_download_directory = "./proteins"
_do_overwrite = True


def download_proteins(taxonomy_ids: list[int] = _taxonomy_ids,
                      gene_names: list[str] = _gene_names,
                      directory: str = _download_directory) -> None:

    files_list = os.listdir(directory)

    query = ""
    if taxonomy_ids:  # check if list is empty
        query = " AND ".join([f"(taxonomy_id: {taxonomy_id})" for taxonomy_id in taxonomy_ids])

    for gene_name in gene_names:

        filename = f"{gene_name.lower()}.fasta"

        if _do_overwrite or filename not in files_list:
            params['query'] = f"(gene:{gene_name})" if query is not "" else f"{query} AND (gene:{gene_name})"
            print(f"search query: {params['query']}")

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
    download_proteins(_taxonomy_ids, _gene_names)
