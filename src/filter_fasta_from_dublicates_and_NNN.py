#!/usr/bin/python3
# usage: python3 {sys.argv[0]} path_to_input_fasta [path_to_filtered_output_fasta]

"""
drop seqs from initial fasta file before multiple alignnment
"""

import sys


def read_fasta_generator(filepath):
    """read fasta without deleting '\n' from line ends to write that
    in the next step
    """
    with open(filepath) as fin:
        seq = ''
        for line in fin:
            if line.startswith(">"):
                if seq != '':
                    yield header, seq
                header = line
                seq = ''
            else:
                seq += line
        yield header, seq


def filter_NNN_and_dublicates_and_write(input_filepath, output_filepath):
    known_seqs, known_uids = set(), set()
    base_seq_set = {'\n', 'A', 'C', 'G', 'T'}
    with open(output_filepath, 'w') as fout:
        for header, seq in read_fasta_generator(input_filepath):
            uid = header.split("|")[1]
            # if "N" in seq:
            if set(seq) != base_seq_set:
                continue
            if uid in known_uids or seq in known_seqs:
                continue
            fout.write(header)
            fout.write(seq)


def main():
    if len(sys.argv) == 1:
        print(
            "usage: python3 {sys.argv[0]} path_to_input_fasta [path_to_filtered_output_fasta]"
        )
        return
    inp_fasta_path = sys.argv[1]
    if len(sys.argv) > 2:
        out_fasta_path = sys.argv[2]
    else:
        out_fasta_path = inp_fasta_path.replace('.fasta', '.filtered.fasta')
    filter_NNN_and_dublicates_and_write(inp_fasta_path, out_fasta_path)


if __name__ == "__main__":
    main()
