"""merge reliable ss and part of rough ss, that contains intergenes
and many free nucleotides (our control)

Посмотрел на грязные данные вторички. Там есть еще 2 межгенника. Однако 1 из них 
несет в себе стемы, сформированные через-чур далекими регионами, не верю им.
Второй нормальный, выделил кусок для данного межгенника и вокруг и руками вставил 
в надежную вторичку. 

Слегка дополненная вторичка лежит в `../data/structure_data/SARSCoV2-RNA_from_site_plus_0based.ss`

"""

from ss2genome_structures import ss_to_df, plot_ss_distribution
import sys
import pandas as pd

sys.path.append('../sourse/')
import utils
from utils import (
    read_ss_file,
    SEC_STRUCT_PATH_v1,
    SEC_STRUCT_PATH_from_site,
)

ANNOTATION_PATH = "../data/genes_annotation.csv"


def _print_intergenes_info(ss_clusters_reliable, ss_clusters_rough, df_annot):
    print('Reliable:')
    for i, fragment in enumerate(ss_clusters_reliable):
        start, stop = fragment[1:3]
        genome_structures = df_annot[
            (df_annot.Pos > start) &
            (df_annot.Pos <= stop + 1)
        ].GenName.unique()
        igenes = [x for x in genome_structures if '_' in x]
        if len(igenes) > 0:
            print(i, fragment[0], igenes)

    print('Rough:')
    for i, fragment in enumerate(ss_clusters_rough):
        start, stop = fragment[1:3]
        genome_structures = df_annot[
            (df_annot.Pos > start) &
            (df_annot.Pos <= stop + 1)
        ].GenName.unique()
        igenes = [x for x in genome_structures if '_' in x]
        if len(igenes) > 0:
            print(i, fragment[0], igenes)

    intergenes_names = [g for g in df_annot.GenName.unique() if "_" in g]
    print()
    print(df_annot.groupby('GenName').Pos.describe()[['min', 'max']].loc[intergenes_names].astype('int'))


def release_sub_ss(ss, start, end):
    """ return part of given ss-fragment

    params:
        - ss: fragment info
        - start: start of releasing region, 0-based, including
        - end: end of releasing region, 0-based, excluding
    """
    new_start = ss[1] + start
    new_end = ss[1] + end - 1
    new_name = f">consensus_{new_start + 1}-{new_end + 1}"
    pos_variants = set(range(start, end))
    new_pairs = [(x[0] - start, x[1] - start) for x in ss[3] if x[0] in pos_variants or x[1] in pos_variants]
    new_pairs_dict = utils.dict_from_pairs(new_pairs)
    new_sstypes = ss[5][start: end]
    new_ss_seq = ss[6][start: end]

    return (
        new_name, 
        new_start,
        new_end,
        new_pairs,
        new_pairs_dict,
        new_sstypes,
        new_ss_seq,
    )

df_annot = pd.read_csv(ANNOTATION_PATH)
ss_clusters_rough = read_ss_file(SEC_STRUCT_PATH_v1)
ss_clusters_reliable = read_ss_file(SEC_STRUCT_PATH_from_site, True)

# utils.print_ss_fasta_to_draw(ss_clusters_rough, 59, 61)

new_region1_ss = release_sub_ss(ss_clusters_rough[59], 385, 438)
new_region2_ss = release_sub_ss(ss_clusters_rough[60], 27, 410)

utils.print_ss_fasta_to_draw([new_region1_ss, new_region2_ss], 0)


# ss_clusters_reliable.append(new_region2_ss)  # very large region. НЕ ВЕРЮ!

# ss_clusters_reliable.append(new_region1_ss)
# df_annot = ss_to_df(df_annot, ss_clusters_reliable)
# plot_ss_distribution(df_annot)

# print(df_annot[(df_annot.Pos > 26218) & (df_annot.Pos <= 26271)])


ss_new = read_ss_file('../data/structure_data/SARSCoV2-RNA_from_site_plus_0based.ss', True)

utils.print_ss_fasta_to_draw(ss_new, 72)

df_annot = ss_to_df(df_annot, ss_new)
plot_ss_distribution(df_annot)



# Reliable:
# 69 >motif_25379-25406 ['S-ORF3a_space']
# 86 >motif_29547-29612 ['N-ORF10_space']
# Rough:
# 57 >consensus_25196-25578 ['S-ORF3a_space']
# 59 >consensus_25834-26427 ['ORF3a-E_space'] - this 26221  26244
# 60 >consensus_26463-27101 ['E-M_space'] - this 26473  26522
# 68 >consensus_29396-29547 ['N-ORF10_space']
# 69 >consensus_29548-29870 ['N-ORF10_space']

#                        min    max
# GenName                          
# ORF1ab-S_space       21556  21562
# S-ORF3a_space        25385  25392
# ORF3a-E_space        26221  26244
# E-M_space            26473  26522
# M-ORF6_space         27192  27201
# ORF6-ORF7a_space     27388  27393
# ORF7a-ORF7b_overlay  27756  27759
# ORF7b-ORF8_space     27888  27893
# ORF8-N_space         28260  28273
# N-ORF10_space        29534  29557