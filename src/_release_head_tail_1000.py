from queue import Queue
from Bio import SeqIO

NSEQS = 1000
PATH_TO_SORTED_HEADERS = "../data/headers_sorted_by_date.txt"
PATH_TO_MULAL = "../data/mulal_gisaid_2021-01-22.filtered.twice.fasta"
PATH_TO_OUTPUT = "../data/first_last_1000_mulal.fasta"


def release_head_tails():
    collection_of_headers = [">NC_045512.2\n"]
    Q = Queue(maxsize=NSEQS)

    with open(PATH_TO_SORTED_HEADERS) as fin:
        while len(collection_of_headers) < NSEQS + 1:
            header = fin.readline()
            if header.count('|') != 2:
                continue
            date = header.strip().split('|')[-1]
            if date.count('-') != 2:
                continue
            collection_of_headers.append(header)

        for header in fin:
            if header.count('|') != 2:
                continue
            date = header.strip().split('|')[-1]
            if date.count('-') != 2:
                continue

            if Q.full():
                Q.get()

            Q.put(header)

        while not Q.empty():
            header = Q.get()
            collection_of_headers.append(header)

    return collection_of_headers


def write_to_file(filepath, collection_of_headers):
    reader = SeqIO.parse(PATH_TO_MULAL, 'fasta')
    seqs = {header: '' for header in collection_of_headers}
    for record in reader:
        header = '>' + record.id + '\n'
        seq = str(record.seq) + '\n'

        if header in collection_of_headers:
            seqs[header] = seq

    output = [x for hs in seqs.items() for x in hs]

    with open(PATH_TO_OUTPUT, 'w') as fout:
        fout.write(''.join(output))


def main():
    collection_of_headers = release_head_tails()
    write_to_file(PATH_TO_OUTPUT, collection_of_headers)


if __name__ == '__main__':
    main()
