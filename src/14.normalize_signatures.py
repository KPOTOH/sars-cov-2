
from collections import defaultdict
import json

PATH_TO_HUMAN_COUNTS = "./data/signatures/GRCh37_codon_counts.json"
PATH_TO_SIGNATURES = "./data/"

def read_triplet_counts(path: str):
    """ read and collapse raw trinucleotide counts """
    with open(path) as fin:
        counts = json.load(fin)
    
    new_counts = defaultdict(int)
    for trinuc, num in counts.items():
        standart_trinuc = trinuc.upper()
        if len(set(standart_trinuc).difference("ATGC")) == 0:
            new_counts[standart_trinuc] += num
    return new_counts


def main():
    triplets = read_triplet_counts(PATH_TO_HUMAN_COUNTS)



if __name__ == "__main__":
    main()
