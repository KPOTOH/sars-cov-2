import imp
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

PATH_TO_GB = "./data/share/NC_045512.2.gb"

with open(PATH_TO_GB) as fin:
    rec: SeqRecord = None
    fea: SeqFeature = None
    for rec in SeqIO.parse(fin, format="genbank"):
        for fea in rec.features:
            if fea.type == "gene":
                gene = fea.qualifiers["gene"][0]
                seqrec = fea.extract(rec.seq)

                print()
        