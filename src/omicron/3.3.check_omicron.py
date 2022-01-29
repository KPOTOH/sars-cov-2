import os
from multiprocessing import Process, Pool
from typing import List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import tqdm

from utils import release_mutations_from_two_seqs


PATH_TO_MULAL = "../data/omicron/mulal.fasta"
PATH_TO_REF = "../data/raw/ideal_table.csv"
NSEQS = 133493
THREADS = 24

reader = SeqIO.parse(PATH_TO_MULAL, 'fasta')
ref = next(reader)


def get_nums(rec):
    seq1, seq2 = str(ref.seq), str(rec.seq)
    mutations = release_mutations_from_two_seqs(seq1, seq2)
    return mutations, rec.name


def release_mutations():
    with Pool(THREADS) as p:
        muts = p.map(get_nums, reader)
    return muts


def load_ref_annot(path_to_ref):
    cols = [
        'Pos', 'GenName', 'GenType', 'CodonNumber', 'RefNuc', 
        'AltNuc', 'RefCodon', 'AltCodon', 'NucInCodon', 'RefAa', 'AltAa'
    ]
    ref = pd.read_csv(path_to_ref, usecols=cols)

    # TODO drop columns and add nsp labels and aa positions
    # TODO from table to set of distinct signatures

    return ref


def process_one_mut(mut: List[tuple]) -> float:
    """ check one record of mutations and return prob of "omicronity" """

    for pos, refnuc, altnuc in mut:
        Pos = pos + 1

    return 0.0 

def check_mutations(muts):
    ref_annot = load_ref_annot(PATH_TO_REF)

    for mut, name in muts:
        pass  # TODO