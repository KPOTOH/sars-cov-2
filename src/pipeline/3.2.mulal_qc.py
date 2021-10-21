from os import write
import random

import click
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from src.utils import count_two_seqs_diff

PATH_TO_MULAL_IN = "./data/mulal.fasta"
PATH_TO_MULAL_OUT = "./data/mulal.filtered.fasta"
MAX_MUT_NUM = 80
DROP_PROB = 0.6
NSEQS = 1139387 # number of records in input fasta


def get_mut_num(rec: SeqRecord, ref: SeqRecord) -> int:
    """ return number of mutations relatively to refseq """
    mut_num = count_two_seqs_diff(str(ref.seq), str(rec.seq))
    return mut_num


def mulal_filtrator(inpath: str, max_mut: int, drop_prob: float) -> SeqRecord:
    """
    generator of filtered records

    filtration:
        - drop seqs highly different from refseq
        - drop 60% seqs by random 

    params:
        inpath - path to input multiple alingnment fasta
    """
    reader = SeqIO.parse(inpath, "fasta")
    ref = next(reader)

    for rec in tqdm(reader, total=NSEQS):
        mut_num = get_mut_num(rec, ref)
        if mut_num > max_mut:
            continue
        if random.random() < drop_prob:
            continue

        yield rec


def fasta_writer(seqs, handle):
    """ write records to file

    params:
        seqs - iterable of seqs
        handle - filepath or file handle
    """
    SeqIO.write(seqs, handle, "fasta-2line")


@click.command("qc")
@click.option("--inpath", default=PATH_TO_MULAL_IN, help="path to input mulal fasta")
@click.option("--outpath", default=PATH_TO_MULAL_OUT, help="path to output mulal fasta")
@click.option("--max-mut", default=MAX_MUT_NUM, help="maximum allovable num of mutations")
@click.option("--drop-prob", default=DROP_PROB, help="probability to drop random sequence")
def main(inpath: str, outpath: str, max_mut: int, drop_prob: float):
    filtered_seqs = mulal_filtrator(inpath, max_mut, drop_prob)
    print(type(filtered_seqs))
    fasta_writer(filtered_seqs, outpath)


if __name__ == "__main__":
    main()
