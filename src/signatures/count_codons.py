

from collections import defaultdict
import json

PATH_TO_GENOME = "./data/human/GCF_000001405.25_GRCh37.p13_genomic.fna"
PATH_TO_JSON_OUT = "./data/human/codon_counts.json"


def one_seq_count(seq: str, counts: defaultdict):
    for i in range(len(seq) - 2):
        codon = seq[i: i + 3]
        if "N" in codon or "n" in codon:
            continue
        counts[codon] += 1


def main():
    chr_num = 0
    codon_counts = defaultdict(int)
    with open(PATH_TO_GENOME) as fasta:
        for line in fasta:
            # if chr_num == 2:
            #     break

            if line.startswith(">"):
                header = line.strip("\n>")
                print("Processing", header)
                prefix = ""
                chr_num += 1
            else:
                seq =  prefix + line.strip()
                if set(seq) == set("N"):
                    continue

                for i in range(len(seq) - 2):
                    codon = seq[i: i + 3]
                    if "N" in codon or "n" in codon:
                        continue
                    codon_counts[codon] += 1

                prefix = seq[-2:]

    print(codon_counts)
    with open(PATH_TO_JSON_OUT, "w") as fout:
        json.dump(codon_counts, fout)


if __name__ == "__main__":
    main()
