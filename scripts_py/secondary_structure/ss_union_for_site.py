#!/usr/bin/python3
# usage: python3 ./ss_union.py ss_file [ss_file]* > union_file

import sys
import re

assert len(sys.argv) > 1, 'arguments are nesessary'

pattern = re.compile("motif_(.+?)\.db")

filenames = sys.argv[1:]
union_file_content = []

for fp in filenames:
    positions = re.search(pattern, fp).groups()[0]
    with open(fp) as fin:
        header = fin.readline()
        seq = fin.readline()
        ss_seq = fin.readline()
        
        union_file_content.append((header, ss_seq))

union_file_content = sorted(
    union_file_content, 
    key=lambda p: int(p[0].lstrip(">motif_").split('-')[0])
)
union_file_content = [x for p in union_file_content for x in p]  # linearization

content = ''.join(union_file_content)
sys.stdout.write(content)
