#!/usr/bin/python3
# usage: python3 ./ss_union.py ss_file [ss_file]* > union_file

import sys
import re

assert len(sys.argv) > 1, 'arguments are nesessary'

choose_score_header = ">consensus_pairs_unique_single_line_0.9"
pattern = re.compile("SARSCoV2-RNA-(.+?)_consensus_1000\.ss")

filenames = sys.argv[1:]
union_file_content = []

for fp in filenames:
	positions = re.search(pattern, fp).groups()[0]
	with open(fp) as fin:
		line = ""
		while not line.startswith(choose_score_header):
			line = fin.readline()

		ss_seq = fin.readline()
		new_header = f">consensus_{positions}\n"

		union_file_content.append(new_header)
		union_file_content.append(ss_seq)


# union_fasta_ss = "./SARSCoV2-RNA.ss"
content = ''.join(union_file_content)
sys.stdout.write(content)
