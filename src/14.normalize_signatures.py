from collections import defaultdict
import json

import numpy as np
import pandas as pd

PATH_TO_HUMAN_COUNTS = "./data/signatures/GRCh37_codon_counts.json"
PATH_TO_SIGNATURES = "./data/external/reference_sig_GRCh37/COSMIC_v3.2_SBS_GRCh37.txt"
PATH_TO_SIGNATURES_NORMALIZED = "./data/signatures/COSMIC_v3.2_SBS_GRCh37_normalized.txt"


def load_triplet_counts(path: str):
    """ read and collapse raw trinucleotide counts """
    with open(path) as fin:
        counts = json.load(fin)

    new_counts = defaultdict(int)
    for trinuc, num in counts.items():
        standart_trinuc = trinuc.upper()
        if len(set(standart_trinuc).difference("ATGC")) == 0:
            new_counts[standart_trinuc] += num
    return new_counts


def load_reference_sig(path: str):
    df = pd.read_csv(path, sep="\t")
    return df


def main():
    triplets = load_triplet_counts(PATH_TO_HUMAN_COUNTS)
    signatures = load_reference_sig(PATH_TO_SIGNATURES)

    ref_triplet = signatures.Type.apply(lambda s: "".join([s[0], s[2], s[-1]]))
    counts_raw = ref_triplet.map(triplets)
    counts = counts_raw / counts_raw.sum()

    signatures = signatures.set_index("Type")
    signatures_normed = signatures.div(counts.values, axis=0)
    signatures_normed = signatures_normed.div(
        signatures_normed.sum(axis=1).values, axis=0)

    signatures_normed.reset_index().to_csv(PATH_TO_SIGNATURES_NORMALIZED, index=None, sep="\t")


if __name__ == "__main__":
    main()
