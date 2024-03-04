#!/usr/bin/python3

import requests
import time
import Bio.SeqIO.FastaIO as FastaIO

# genes of interest
gene_names = ["Period", "Timeless", "Clock"]

# Gene name | Synonyms
# timeless  | tim, tim1, tim-1, timeless1
# period    | per
# clock     | clk
# more genes to be added

#   CWO: Clockwork Orange
#   VRI: Vrille
# PDP1e: Par domain protein 1e
#   SGG: Shaggy
#   DBT: Doubletime
#  CRY1 / d-CRY:  Cryptochrome in Drosophilia
#  CRY2 / m-CRY:  Mammalian-type Chryptochrome

# targets
taxonomy_ids = [50557]
#              true insects

# Uniprot API
# for reference: https://www.uniprot.org/help/api_queries
url = "https://uniprot.org/uniprotkb/search"
params = {
    '?format': 'fasta',
    # because requests is fucking stupid and removes the ? from urls as it thinks that the query might be empty...
    'compressed': 'false'
}
directory = "./fasta"


def main():
    # get all genomes as jason

    for taxonomy_id in taxonomy_ids:
        for gene_name in gene_names:

            params['query'] = f"(taxonomy_id:{taxonomy_id}) AND (gene:{gene_name})"

            with requests.get(url, params=params) as response:  # stream=True
                print(response.url)  # for testing if url is correctly encoded.
                print(response.status_code)
                response.raise_for_status()
                # print(response.text)
                filename = f"{directory}/{gene_name}_{taxonomy_id}.fasta"
                print("file: ", filename)

                with open(filename, 'w') as file:
                    file.write(response.text)

                with open(filename, 'r') as file:
                    for (headline, sequence) in FastaIO.SimpleFastaParser(file):
                        print(headline)

            time.sleep(5)  # rate limiting


if __name__ == '__main__':
    main()
