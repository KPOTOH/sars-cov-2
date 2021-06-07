"""
read tree and multiple alignment from prank output and
write json file that contains all horisontal substitutions from tree

json format:
[
    (parent node,
    child node,
    substitution [(position_0_based, parent_nucl, child_nucl), ...]), ...
]
"""

import json
import os
from queue import Queue
import sys

from ete3 import PhyloTree

PRANK_OUTPUT_PATH = "../data/mulal_gisaid_2021-01-22.filtered.fasta.prank.anc"
TREE_PATH = PRANK_OUTPUT_PATH + ".dnd"
FASTA_PATH = PRANK_OUTPUT_PATH + ".fas"

MUTATION_STORAGE_PATH = '../data/overall_mutations_3.json'


def read_fasta_generator(filepath: str):
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


fasta_generator = read_fasta_generator(FASTA_PATH)
fasta_storage = dict()


def get_sequence(node_name: str) -> str:
    """read fasta file online and write seqs to dict if another one are needed"""
    if node_name in fasta_storage:
        seq = fasta_storage[node_name]
        return seq

    for new_node_name, new_seq in fasta_generator:
        new_node_name = new_node_name.strip(">\n")
        new_seq = new_seq.replace("\n", "")
        fasta_storage[new_node_name] = new_seq
        if new_node_name == node_name:
            return new_seq


def node_parent(node):
    try:
        return next(node.iter_ancestors())
    except BaseException:
        return None


def trim_two_seqs(seq1, seq2):
    """there are '------' in the start and end of seqs. Drop it!"""
    n = len(seq1)
    start_pos, stop_pos = 0, n
    for start_pos in range(n):
        if not (seq1[start_pos] == '-' or seq2[start_pos] == "-"):
            break

    for stop_pos in range(n - 1, -1, -1):
        if not (seq1[stop_pos] == '-' or seq2[stop_pos] == "-"):
            break
    stop_pos += 1
    return start_pos, stop_pos


def release_mutations_from_two_seqs(parent_seq: str, child_seq: str):
    start_pos, stop_pos = trim_two_seqs(parent_seq, child_seq)

    mutations = []
    for pos in range(start_pos, stop_pos):
        sourse_nucl = parent_seq[pos]
        mutant_nucl = child_seq[pos]
        if sourse_nucl != mutant_nucl:
            mutations.append((pos, sourse_nucl, mutant_nucl))
    return mutations


def get_mutations(node):
    parent = node_parent(node)
    if parent is None:
        return
    seq_of_parent = get_sequence(parent.name)
    seq_of_child = get_sequence(node.name)
    assert len(seq_of_child) == len(seq_of_parent), (
        "parent and child seq lenghts aren't equal"
    )
    mutations = release_mutations_from_two_seqs(
        seq_of_parent,
        seq_of_child,
    )
    return parent.name, node.name, mutations


def bfs_to_extract_mutspec(path_to_tree):
    """ BFS for extraction """
    tree = PhyloTree(path_to_tree, format=1)
    print("tree readed", file=sys.stderr)
    discovered_nodes = set()
    Q = Queue()
    Q.put(tree)
    discovered_nodes.add(tree.name)

    overall_mutations = []
    while not Q.empty():
        cur_node = Q.get()
        for child in cur_node.children:
            Q.put(child)

        if cur_node.name not in discovered_nodes:
            discovered_nodes.add(cur_node.name)

            # main process starts here
            cur_mutations_of_one = get_mutations(cur_node)
            overall_mutations.append(cur_mutations_of_one)
    print(
        f"n_discovered_nodes: {len(discovered_nodes)}, "
        f"n_overall_mutations: {len(overall_mutations)}",
        file=sys.stderr
    )

    return overall_mutations


if __name__ == '__main__':
    overall_mutations = bfs_to_extract_mutspec(TREE_PATH)

    with open(MUTATION_STORAGE_PATH, 'w') as fout:
        json.dump(overall_mutations, fout)
