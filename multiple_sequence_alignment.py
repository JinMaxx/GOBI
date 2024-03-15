#!/usr/bin/python3

import os
import subprocess

# doas overwrite files!

mafft = "mafft-linsi"
aliview = "./aliview/aliview"

_input_directory = "./proteins_filtered"  # or ./proteins
_output_directory = "./proteins_aligned"


def generate_alignments(input_directory: str = _input_directory,
                        output_directory: str = _output_directory):

    output_directory_files = os.listdir(output_directory)

    for input_file in os.listdir(input_directory):
        if input_file not in output_directory_files:

            output_file = f"{output_directory}/{input_file}"
            input_file = f"{input_directory}/{input_file}"

            print(mafft, input_file, ">", output_file)
            with open(output_file, 'w') as output_file:
                subprocess.run([mafft, input_file], stdout=output_file)


def display_alignments(alignments_directory: str = _output_directory):
    for output_file in os.listdir(alignments_directory):
        output_file = f"{alignments_directory}/{output_file}"
        subprocess.run(["bash", aliview, output_file])


if __name__ == '__main__':
    generate_alignments(_input_directory, _output_directory)
    display_alignments(_output_directory)
