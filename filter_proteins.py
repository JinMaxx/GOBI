#!/usr/bin/python3

import os
import Bio.SeqIO as SeqIO
import Bio.SeqIO.FastaIO as FastaIO


# will overwrite files in filtered_proteins

_filter_keywords = ["isoform", "kaput", "fragment"]  # can also use id like "tr|A0A034V0Z2|A0A034V0Z2_BACDO"


# No need to change those

_proteins_dir = "./proteins"
_proteins_filtered_dir = "./proteins_filtered"


def filter_proteins(filter_keywords: list[str] = _filter_keywords,
                    input_proteins_dir: str = _proteins_dir,
                    output_proteins_dir: str = _proteins_filtered_dir):

    for file_input in os.scandir(input_proteins_dir):

        if file_input.is_file() and file_input.name.endswith(".fasta"):
            file_output_relative_name = f"{output_proteins_dir}/{file_input.name}"
            # open(file_output_relative_name, 'a').close()  # create empty file

            with open(file_input) as input_file_handle, open(file_output_relative_name, "w") as output_file_handle:
                writer = FastaIO.FastaWriter(output_file_handle)

                for record in SeqIO.parse(input_file_handle, "fasta"):
                    if any(keyword in record.description for keyword in filter_keywords):
                        print(f"removing: {record.description}")
                        continue  # skips if fasta title contains one keyword
                    else:
                        # TODO remove signal peptides with SignalP6 here
                        #  get software with their fucking stupid registration form
                        writer.write_record(record)


if __name__ == '__main__':
    filter_proteins(_filter_keywords)
