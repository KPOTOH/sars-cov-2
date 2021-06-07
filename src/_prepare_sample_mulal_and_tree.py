"""
create sample multiple alignment and tree from big dataset that doesn't
fit into 500gb memory

`n` is covid count in sample dataset

useless after reading docs of prank: we can select flag to use only seqs in tree
"""
import os
from random import randint

from ete3 import PhyloTree

from filter_fasta_from_dublicates_and_NNN import read_fasta_generator

n = 50_000

mulal_path = "../mulal_gisaid_2021-01-22.filtered.fasta"
tree_path =  "../mulal_gisaid_2021-01-22.filtered.fasta.tre"

assert os.path.exists(mulal_path)
assert os.path.exists(tree_path)

fasta = read_fasta_generator(mulal_path)
tree = PhyloTree(tree_path)
full_nodes = set(tree.get_leaf_names())

sample_mulal_path = f"sample_mulal_{n}.fasta"
sample_tree_path = f"sample_tree_{n}.tre-simple"

node_names = []
with open(sample_mulal_path, 'w') as fout:
    for i, (header, seq) in enumerate(fasta):
        # simple randomizer, reference including
        if i == 0 or randint(0, 1) == 0:
            name = header.strip('\n>')
            if name not in full_nodes:
                continue

            node_names.append(name)
            fout.write(header)
            fout.write(seq)

            if len(node_names) == n:
                break

tree.prune(node_names)
tree.write(outfile=sample_tree_path, format=9) # tre-simple
