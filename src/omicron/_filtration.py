
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

PATH_TO_FASTA = "./data/gisaid_omicron_13-01-22.fasta"
PATH_TO_OUT = "./data/mod_gisaid_omicron_13-01-22.fasta"


def main():
    sequences = SeqIO.parse(PATH_TO_FASTA, "fasta")
    fcollection = set()

    rec: SeqRecord = None
    with open(PATH_TO_OUT, "w") as fout:
        for rec in sequences:
            if rec.name in fcollection:
                continue

            fcollection.add(rec.name)
            SeqIO.write(rec, fout, "fasta")


if __name__ == "__main__":
    main()
