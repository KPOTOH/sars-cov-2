
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
import pandas as pd


PATH_TO_FASTA = "./data/gisaid_omicron_13-01-22.fasta"
PATH_TO_OUT = "./data/gisaid_omicron_13-01-22.filtered.fasta"

MIN_SEQ_LEN = 29001
N_SHARE_CUTOFF = 0.01


def main():
    sequences = SeqIO.parse(PATH_TO_FASTA, "fasta")
    ndrops = 0
    with open(PATH_TO_OUT, "w") as fout:
        for rec in tqdm(sequences, total=262000):
            seq = str(rec.seq)
            n_count = seq.count("N")
            lenght = len(seq)
            if lenght < MIN_SEQ_LEN or n_count / lenght > N_SHARE_CUTOFF:
                ndrops += 1
                continue

            SeqIO.write(rec, fout, "fasta")

    print("Dropped {} seqs".format(ndrops))


if __name__ == "__main__":
    main()
