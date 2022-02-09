"""plot known ss distribution on genome structures"""
from collections import defaultdict
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.append('../sourse/')
from utils import read_ss_file

SEC_STRUCT_PATH_v1 = "../data/structure_data/SARSCoV2-RNA.ss"
SEC_STRUCT_PATH_from_site = "../data/structure_data/SARSCoV2-RNA_from_site_0based.ss"
ANNOTATION_PATH = "../data/genes_annotation.csv"


def indicate_pairs_alogn_genome(ss_clusters):
    genome_is_stem = defaultdict(int)
    genome_ss_pairs = defaultdict(lambda: None)

    for j in range(len(ss_clusters)):
        cur_ss = ss_clusters[j]
        start, stop = cur_ss[1:3]
        types = cur_ss[5]
        pairs_dct = cur_ss[4]

        for i in range(start, stop + 1):
            ss_idx = i - start
            genome_is_stem[i] = 0
            if types[ss_idx] == 1:
                genome_is_stem[i] = 1
                genome_ss_pairs[i] = pairs_dct[ss_idx] + start
    return genome_is_stem, genome_ss_pairs


def plot_ss_distribution(df):
    # https://matplotlib.org/stable/gallery/lines_bars_and_markers/barchart.html

    def autolabel(rects):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for rect in rects:
            height = rect.get_height()
            ax.annotate('{:.2f}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')

    # labels = (gr.get_group(0).GenName.value_counts() /
    #           df.GenName.value_counts()).index
    labels = [
        'ORF1ab', 'S', 'ORF3a', 'E', 'M', 'ORF6', 'ORF7a', 'ORF7b',
        'ORF8', 'N', 'ORF10', '5UTR', '3UTR', 'ORF1ab-S_space', 'S-ORF3a_space',
        'ORF3a-E_space', 'E-M_space', 'M-ORF6_space', 'ORF6-ORF7a_space',
        'ORF7a-ORF7b_overlay', 'ORF7b-ORF8_space', 'ORF8-N_space', 'N-ORF10_space',
    ]
    ss_structures = [
        "free or loop",
        "stem",
    ]

    gr = df.groupby("IsStem")
    known_ss_shares = [(gr.get_group(i).GenName.value_counts() /
                        df.GenName.value_counts())[labels].values for i in range(gr.ngroups)]

    x = np.arange(len(labels))  # the label locations
    width = 0.25  # the width of the bars

    fig, ax = plt.subplots(figsize=(14, 6))
    collection_of_rects = []
    for i, v in enumerate(known_ss_shares):
        rects = ax.bar(x + (-1)**i * width / 2, v,
                       width, label=ss_structures[i])
        collection_of_rects.append(rects)

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Share of nucl amount')
    ax.set_xlabel("Genome structures")
    ax.set_title(
        'Amount of nucleotides in stem or not across all genome structures')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=55)
    ax.legend()

    for r in collection_of_rects:
        autolabel(r)

    fig.tight_layout()
    plt.show()


def ss_to_df(df, ss_clusters):
    df = df.copy()
    genome_is_stem, genome_ss_pairs = indicate_pairs_alogn_genome(ss_clusters)
    df['IsStem'] = pd.Series(genome_is_stem)
    df['SsPairs'] = pd.Series(genome_ss_pairs)

    # check if ss_pairs is not nan only where ss_type == 1
    assert ~np.any(df[df.IsStem == 1].SsPairs.isna())
    assert np.any(df[df.IsStem != 1].SsPairs.isna())

    return df


def main():
    ss_clusters = read_ss_file(SEC_STRUCT_PATH_from_site)
    df = pd.read_csv(ANNOTATION_PATH)
    df = ss_to_df(df, ss_clusters)
    plot_ss_distribution(df)


if __name__ == "__main__":
    main()
