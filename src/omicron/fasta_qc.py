
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
import pandas as pd


PATH_TO_FASTA = "./data/gisaid_omicron_13-01-22.fasta"
PATH_TO_OUT = "./data/gisaid_omicron_13-01-22.filtered.fasta"

MIN_SEQ_LEN = 29001


def main():
    sequences = SeqIO.parse(PATH_TO_FASTA, "fasta")
    n_contain_N = 0
    n_less_min = 0
    n_amounts = []

    rec: SeqRecord = None
    with open(PATH_TO_OUT, "w") as fout:
        for rec in tqdm(sequences, total=262000):
            seq = rec.seq
            n_contain_N += int("N" in seq)
            n_less_min += int(len(seq) < MIN_SEQ_LEN)
            n_amounts.append(seq.count("N"))


            # SeqIO.write(rec, fout, "fasta")

    print(n_less_min, n_contain_N)
    print(pd.Series(n_amounts).value_counts().sort_index())


if __name__ == "__main__":
    main()
