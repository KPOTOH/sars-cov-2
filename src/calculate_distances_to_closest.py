from collections import defaultdict
import json
import sys
import time

import numpy as np
import networkx as nx
import pandas as pd
from rna_tools.SecondaryStructure import parse_vienna_to_pairs
import tqdm

sys.path.append('../sourse/')
from mutations_extractor_with_context import get_sequence, read_fasta_generator


REFSEQ_PATH = "../data/covid_ref.fasta"
SUBSTITUTIONS_PATH = "../data/overall_mutations_with_context2.json"
SEC_STRUCT_PATH = '../data/structure_data/SARSCoV2-RNA_from_site_plus_0based.ss'


def dict_from_pairs(pairs: list):
    """there is a list of pairs, that will be converted to
    dict, where each element of pair is key and value at the same time:
    {k: v, v: k} for any (k, v) in pairs
    """
    d1 = {k: v for k, v in pairs}
    d2 = {v: k for k, v in pairs}
    d1.update(d2)
    return d1


def assign_ss_types(start, stop, paired_pos_dict):
    """ return the type of nucleotides in ss;

    type variants:
    0: free
    1: stem
    2: hairpin loop
    3: bulge loop or internal loop

    @param ss_idx, int - index of one secondary cluster from output of `read_ss_file` func
    @param nucl_idx, int - index of nucleotide in the ss
    """
    n = stop - start + 1
    inside_loop = False
    ss_types = [0] * n
    cur_loop_stop = -1

    for pos in range(n):
        if pos in paired_pos_dict:
            if not inside_loop:
                cur_loop_stop = paired_pos_dict[pos]
                inside_loop = True

            cur_ss_type = 1
            last_stem_pos = pos
        else:
            if inside_loop:
                next_stem_pos = pos + 1
                # here we find next stem pos
                while next_stem_pos not in paired_pos_dict:
                    next_stem_pos += 1

                if paired_pos_dict[last_stem_pos] == next_stem_pos:
                    cur_ss_type = 2
                else:
                    cur_ss_type = 3
            else:
                cur_ss_type = 0

        ss_types[pos] = cur_ss_type

        if pos == cur_loop_stop:
            inside_loop = False

    return ss_types


def read_ss_file(filepath, from_site=False):
    """read ss file in multiple fasta format and
    return list of clusters, where there are:
    - cluster header,
    - start_pos (0-based),
    - stop_pos (0-based),
    - ss_pairs of paired nucleotides (0-based)
    - dict of pairs in both direction of pair (very usefull)
    - ss sequence"""

    if "from_site" in filepath:
        from_site = True

    clusters = []
    with open(filepath) as fin:
        for line in fin:
            if line == "\n":
                break
            if line.startswith(">"):
                header = line.strip()
            else:
                ss_seq = line.strip()
                ss_pairs = parse_vienna_to_pairs(ss_seq)
                ss_pairs_0_based = [(x - 1, y - 1) for x, y in ss_pairs[0]]
                paired_pos_dict = dict_from_pairs(ss_pairs_0_based)
                start_pos, stop_pos = map(
                    int, header.lstrip(">consensus_").lstrip(">motif_").split("-"))
                if not from_site:
                    start_pos -= 1
                    stop_pos -= 1
                ss_types = assign_ss_types(
                    start_pos, stop_pos, paired_pos_dict)
                clusters.append((
                    header,
                    start_pos,
                    stop_pos,
                    ss_pairs_0_based,
                    paired_pos_dict,
                    ss_types,
                    ss_seq,
                ))
    return clusters


def ss_graph_constructor(ss_cluster):
    n = ss_cluster[2] - ss_cluster[1] + 1
    G = nx.path_graph(n)
    pairs = ss_cluster[3]
    for v1, v2 in pairs:
        G.add_edge(v1, v2)
    return G


def full_ss_graph_constructor(refseq_len, full_genome_pairs):
    n = refseq_len
    G = nx.path_graph(n)
    for v1, v2 in full_genome_pairs.items():
        G.add_edge(v1, v2)
    return G


def primary_distance(pos1, pos2):
    """@param pos1, pos2: int[0..29902], positions of genome using to calculate"""
    return abs(pos1 - pos2)


def secondary_distance(pos1, pos2):
    """calculate distance using secondary structure

    @param pos1, pos2: int[0..29902], positions of genome using to calculate
    distance between them

    distance calculated by dijkstra_path_length in genome graph
    """
    dist = nx.dijkstra_path_length(full_genome_graph, pos1, pos2)
    return dist


def dist2closest_subst(i, substs, history=None):
    """ посчитает расстояние до ближайшей замены, если дать
    на вход индекс интересующей позиции (i) в общем скопе
    замен (substs)

    return (distance1, distance2)
    """
    MAX_PRIMARY_DISTANCE = 30000
    n_substs = len(substs)
    assert 0 <= i < n_substs, 'substs intex out of range'

    if n_substs == 1:
        return None, None

    cur_pos = substs[i][0]
    around_window = 5
    closest_subst_dist1 = closest_subst_dist2 = MAX_PRIMARY_DISTANCE  # max
    clidex_by_primary = clidex_by_secondary = -1
    cur_d2 = 0

    for j in range(max(i - around_window, 0),
                   min(i + around_window, n_substs - 1) + 1):
        if i == j:
            continue
        maybe_cl_pos = substs[j][0]

        if history is not None:
            _pair_pos = (cur_pos, maybe_cl_pos)
            if _pair_pos in history:
                key = _pair_pos
            elif _pair_pos[::-1] in history:
                key = _pair_pos[::-1]
            else:
                key = None
        else:
            key = None

        if key is not None:
            cur_d1, cur_d2 = history[key]
        else:
            cur_d1 = primary_distance(cur_pos, maybe_cl_pos)
            if cur_d1 > MAX_PRIMARY_DISTANCE:
                continue
            # cur_d2 = secondary_distance(cur_pos, maybe_cl_pos)  # UNCOMMENT TO COUNT DIST2
            history[(cur_pos, maybe_cl_pos)] = (cur_d1, cur_d2)
            history[(maybe_cl_pos, cur_pos)] = (cur_d1, cur_d2)

        if cur_d1 > MAX_PRIMARY_DISTANCE:
            continue
        if cur_d1 < closest_subst_dist1:
            # closest substitution index by primary structure
            clidex_by_primary = j
            closest_subst_dist1 = cur_d1

        if cur_d2 < closest_subst_dist2:
            # closest substitution index by secondary structure
            clidex_by_secondary = j
            closest_subst_dist2 = cur_d2

    return closest_subst_dist1, closest_subst_dist2


def is_paired(nucl1: str, nucl2: str):
    """check if given nucleotides can form Watson-Crick pairing
    
    return:
        - 0 if not paired
        - 1 if weak pair (A:T or G:T)
        - 2 if strong pair
    """
    nset = {nucl1, nucl2}
    if nset == {'A', 'T'} or nset == {'T', 'G'}:
        return 1
    elif nset == {'C', 'G'}:
        return 2
    return 0


def determine_substitution_type(subst: tuple, parent_node, child_node):
    """ looking at secondary structure and return type of subst:
    {complementary, not_complementary, not_in_sh}

    @param subst: tuple(pos, parent_nucl, child_nucl)
    to adress `ss_clusters`

    return stype (int): substitution type, one of:
    {   
        0: "no info",
        1: "free>free",
        2: "free>weak",
        3: "free>strong",
        4: "weak>free",
        5: "weak>weak",
        6: "weak>strong",
        7: "strong>free",
        8: "strong>weak",
        9: "strong>strong",
    }
    """
    pos, parent_nucl, child_nucl = subst[0], subst[1], subst[2]
    stype = 0  # no info by default
    if pos in full_genome_pairs:
        # данный код выполняется, если
        # замутировала та, что находится в спаренном состоянии в рефсеке
        # нужно проверить ее тип, разрушила она парную связь иль нет
        the_pair_pos = full_genome_pairs[pos]

        parent_seq = get_sequence(parent_node)
        child_seq = get_sequence(child_node)

        parent_nucl_of_paired = parent_seq[the_pair_pos]
        child_nucl_of_paired = child_seq[the_pair_pos]

        parent_pos_type = is_paired(parent_nucl, parent_nucl_of_paired)
        child_pos_type = is_paired(child_nucl, child_nucl_of_paired)

        if parent_pos_type == 0:
            if child_pos_type == 0:
                stype = 1  # free2free
            elif child_pos_type == 1:
                stype = 2  #free2weak
            elif child_pos_type == 2:
                stype = 3  #free2strong

        elif parent_pos_type == 1:
            if child_pos_type == 0:
                stype = 4  # weak2free
            elif child_pos_type == 1:
                stype = 5  #weak2weak
            elif child_pos_type == 2:
                stype = 6  #weak2strong

        elif parent_pos_type == 2:
            if child_pos_type == 0:
                stype = 7  # strong2free
            elif child_pos_type == 1:
                stype = 8  #strong2weak
            elif child_pos_type == 2:
                stype = 9  #strong2strong
    return stype


def _is_paired(nucl1, nucl2):
    """check if given nucleotides can form Watson-Crick pairing"""
    nset = {nucl1, nucl2}
    if nset == {'A', 'T'} or nset == {'C', 'G'}:
        return True
    return False


def _determine_substitution_type(subst: tuple, parent_node, child_node):
    """ looking at secondary structure and return type of subst:
    {complementary, not_complementary, not_in_sh}

    @param subst: tuple(pos, parent_nucl, child_nucl)
    to adress `ss_clusters`

    return stype (int): substitution type, one of:
    {
    0: комплементарные замены, приводящие к формированию комплементарных взаимодействий
    1: некомплементарные замены, разрушающие уже имеющиеся комплементарные взаимодействия
    2: замены, которые не восстанавливают комп. взаимодействия, хотя в референсе взаимодействие было
    3: вне вторичных взаимодействий
    4: еще есть комплементарные, которые заново сформированы за 1 итерацию
    }
    """
    pos, parent_nucl, child_nucl = subst[0], subst[1], subst[2]

    if pos in full_genome_pairs:
        # данный код выполняется, если
        # замутировала та, что находится в спаренном состоянии в рефсеке
        # нужно проверить ее тип, разрушила она парную связь иль нет
        the_pair_pos = full_genome_pairs[pos]

        parent_seq = get_sequence(parent_node)
        child_seq = get_sequence(child_node)

        parent_nucl_of_paired = parent_seq[the_pair_pos]
        child_nucl_of_paired = child_seq[the_pair_pos]

        # check if parent pair was existed
        if is_paired(parent_nucl, parent_nucl_of_paired):
            if not is_paired(child_nucl, child_nucl_of_paired):
                stype = 1  # parent yes, child no
            else:
                # здесь присходит парная замена с сохранением вторичной связи
                #                 print('DIN-DON!', file=sys.stderr)
                stype = 4  # самый редкий, должно быть
        else:
            # родительские нуклы не были спарены, спарены ли детские?
            if is_paired(child_nucl, child_nucl_of_paired):
                stype = 0  # parent no, child yes
            else:
                stype = 2  # parent no, child no
    else:
        stype = 3

    return stype


if __name__ == "__main__":
    refseq = next(read_fasta_generator(REFSEQ_PATH))[1].replace('\n', '')
    refseq_len = len(refseq)

    with open(SUBSTITUTIONS_PATH) as fin:
        substitutions = json.load(fin)

    # read ss-file
    ss_clusters = read_ss_file(SEC_STRUCT_PATH)

    full_genome_pairs = defaultdict(lambda: None)  # 0-based
    for ss in ss_clusters:
        start = ss[1]
        pairs = ss[4]
        for k, v in pairs.items():
            full_genome_pairs[k + start] = v + start
    full_genome_pairs = dict(full_genome_pairs)


    ss_clusters_graphs = list(map(ss_graph_constructor, ss_clusters))
    full_genome_graph = full_ss_graph_constructor(refseq_len, full_genome_pairs)


    # collection of substs and distance to closest substs
    final_fantasy = []
    history = dict()

    pair_idx = 0
    for parent_node, child_node, substs in tqdm.tqdm(substitutions):
        # only for SNPs, no indels!!!
        substs = [x for x in substs if x[1] != '-' and x[2] != '-']

        if len(substs) < 2:
            # the pair parent-child dropped,
            # because we need calculate closest subst. for analysis
            pair_idx += 1
            continue

        for j in range(len(substs)):
            cur_subst = substs[j]
            cur_subst_type = determine_substitution_type(
                cur_subst, parent_node, child_node,
            )
            # calculate distance to closest subst. by secondary and primary
            # structure
            nearest_sdist1, nearest_sdist2 = dist2closest_subst(j, substs, history)

            final_fantasy.append((
                pair_idx,
                cur_subst[0],
                cur_subst[1],
                cur_subst[2],
                cur_subst_type,
                nearest_sdist1,
                # nearest_sdist2,
                parent_node,
                child_node,
            ))
        pair_idx += 1


    cols = ["pair_idx", 'pos', "parent_nucl", "child_nucl", 
            'stype', 'primary_dist2nearest',
            'parent_node', 'child_node', ]
    cur_time = time.asctime().replace(" ", "_")

    df = pd.DataFrame(final_fantasy, columns=cols)
    df.to_csv(f"../data/new_final_fantasy_{cur_time}.csv", index=None)
    print(df)
