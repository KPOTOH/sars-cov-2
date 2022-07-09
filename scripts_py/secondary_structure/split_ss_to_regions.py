"""
split regions and controls to subregions
"""

import os
import re
import sys
from typing import List

import pandas as pd

import utils

SEC_STRUCT_PATH_PLUS = '../data/structure_data/SARSCoV2-RNA_from_site_plus_0based.ss'
SEC_STRUCT_PATH_ROUGH = '../data/structure_data/SARSCoV2-RNA.ss'

PATH_TO_REGIONS = '../data/ss_regions_for_IWT.csv'
PATH_TO_CONTROL_REGIONS = '../data/ss_control_regions_for_IWT.csv'
PATH_TO_CONTROL_REGIONS_ROUGH = '../data/ss_control_regions_rough_for_IWT.csv'


def one_control_seq(idx, ndots, ss_start):
    regseq = '.' * ndots
    reg_end = idx  # not including
    reg_start = idx - ndots
    return {
        'ss_seq': regseq,
        'start': ss_start + reg_start,
        'end': ss_start + reg_end,
    }


def release_free_bases(ss_seq, min_len, ss_start, constant_lenght=False):
    regions = []
    ndots = 0
    for i, char in enumerate(ss_seq):
        if char == '.':
            ndots += 1
            if ndots >= min_len and constant_lenght:
                regions.append(one_control_seq(i, ndots, ss_start))
                ndots = 0
        else:
            if ndots >= min_len:
                regions.append(one_control_seq(i, ndots, ss_start))
            ndots = 0

    return regions


def control_patririon(ss_clusters, min_len=5):
    regions = []
    for cur_ss in ss_clusters:
        ss_seq = cur_ss[-1]
        ss_start = cur_ss[1]
        parts = release_free_bases(ss_seq, min_len, ss_start)
        for reg in parts:
            regions.append(reg)
    return regions


def is_in_hairpin(ss_seq, pos, nfront_to_see=15):
    assert ss_seq[pos] == '.'
    last_bracket = ss_seq[pos - 1]
    assert last_bracket in '()'

    is_hairpin = False
    h_end_pos = None
    ndots = 1
    for i in range(pos + 1, pos + nfront_to_see):
        if ss_seq[i] == '.':
            ndots += 1
            continue
        if ss_seq[i] != last_bracket:
            is_hairpin = True
            h_end_pos = i
            break
        else:
            break
    return is_hairpin, h_end_pos


def one_element_reg(ss_seq, start, end, ss_start):
    regseq = ss_seq[start: end]
    return {
        'ss_seq': regseq,
        'start': ss_start + start,
        'end': ss_start + end,
    }


def borders(h_ndots, drop_marginal_dots=True):
    # А ЕСЛИ ВСЕГО 3 НУКЛЕОТИДА В ХЭИРПИНЕ???  
    # ЭТУ ДВОЙКУ НАДО ЗАМЕНИТЬ НА УСЛОВИЕ  
    if drop_marginal_dots:
        return 0, 0

    if h_ndots < 2:
        return 0, 0
    if h_ndots == 2:
        return 1, 1
    if h_ndots == 3:
        return 2, 1
    if h_ndots > 3:
        return 2, 2


def split_ss_part1(ss_seq, ss_start, min_len=20):
    """ split on (...) and )( 
            
    end in ss_clusters is included, 
    but here in regions isn't
    """
    min_len = min(min_len, 20)
    seq_len = len(ss_seq)
    n_dots_before = 0
    regions = []
    start = 0
    pos = 0
    while pos < seq_len:
        char = ss_seq[pos]
        if char == '.':
            n_dots_before += 1
            if n_dots_before == 1 and pos != 0 and ss_seq[pos - 1] != '.':
                is_hairpin, h_end_pos = is_in_hairpin(ss_seq, pos, seq_len - pos)
                if is_hairpin:
                    h_ndots = h_end_pos - pos
                    end_delta, st_delta = borders(h_ndots)
                    end = pos + end_delta  # last excluded
                    if end - start >= min_len:
                        reg = one_element_reg(ss_seq, start, end, ss_start)
                        regions.append(reg)
                    start = h_end_pos - st_delta
                    pos = h_end_pos + 1
                    continue
        else:
            n_dots_before = 0
            if ss_seq[pos] == '(' and ss_seq[pos - 1] == ')' and pos != 0:
                end = pos
                if end - start >= min_len:
                    reg = one_element_reg(ss_seq, start, end, ss_start)
                    regions.append(reg)

                start = pos
        pos += 1

    reg = one_element_reg(ss_seq, start, len(ss_seq), ss_start)
    regions.append(reg)
    return regions


def dot_split(reg, ndots=4, drop_marginal_dots=True):
    splitter = '.' * ndots
    ss_seq = reg['ss_seq']
    parts = ss_seq.split(splitter)
    if len(parts) == 1:
        return [reg]

    rstart = reg['start']
    new_regions = []
    for i, pot_reg_seq in enumerate(parts):
        if pot_reg_seq != '':
            if i == 0 and not drop_marginal_dots:
                pot_reg_seq += '..'
            rend = rstart + len(pot_reg_seq)
            # pot_reg_seq.lstrip('.')
            nreg = {'ss_seq': pot_reg_seq, 'start': rstart, 'end': rend}
            new_regions.append(nreg)
            rstart = rend + ndots
        else:
            rstart += ndots
    return new_regions


def half_splitting(reg):
    ss_seq = reg['ss_seq']
    start = reg['start']
    end = reg['end']
    half = (end - start) // 2
    half_pos = start + half
    reg1 = reg.copy()
    reg1['ss_seq'] = ss_seq[:half]
    reg1['end'] = half_pos
    reg2 = reg.copy()
    reg2['ss_seq'] = ss_seq[half:]
    reg2['start'] = half_pos
    return [reg1, reg2]


def ss_partition(ss_clusters: List[tuple], max_len=30, min_len=5):
    hairpin_splitted_regions = []
    for piece in ss_clusters:
        ss_start, ss_seq = piece[1], piece[-1]
        regs_from_piece = split_ss_part1(ss_seq, ss_start, 5)
        hairpin_splitted_regions.extend(regs_from_piece)
    
    dot_splitted_regions = []
    for reg in hairpin_splitted_regions:
        by_4 = dot_split(reg, 4)
        for sreg in by_4:
            by_3 = dot_split(sreg, 3)
            by_3 = [x for x  in by_3 if len(x['ss_seq']) >= min_len]
            dot_splitted_regions.extend(by_3)

        # for sreg in by_4:
        #     if len(sreg['ss_seq']) > max_len:
        #         by_3 = dot_split(sreg, 3)
        #         dot_splitted_regions.extend(by_3)
        #     else:
        #         dot_splitted_regions.append(sreg)

    ## uncomment if there are need to decrease max len
    # half_splitted_regions = []
    # for reg in dot_splitted_regions:
    #     if len(reg['ss_seq']) > max_len:
    #         half_splitted_regions.extend(half_splitting(reg))
    #     else:
    #         half_splitted_regions.append(reg)
    
    return dot_splitted_regions


def regions2csv(regions: List[dict], filepath):
    df = pd.DataFrame(regions)
    df.to_csv(filepath, index=None)


def main():
    ss_clusters = utils.read_ss_file(SEC_STRUCT_PATH_PLUS)
    ss_clusters_rough = utils.read_ss_file(SEC_STRUCT_PATH_ROUGH)

    regions = ss_partition(ss_clusters)
    control_regions_rough = control_patririon(ss_clusters_rough)
    # control_regions = control_patririon(ss_clusters)

    regions2csv(regions, PATH_TO_REGIONS)
    regions2csv(control_regions_rough, PATH_TO_CONTROL_REGIONS_ROUGH)
    # # regions2csv(control_regions, PATH_TO_CONTROL_REGIONS)


if __name__ == '__main__':
    main()




# ss_seq = ss_clusters[1][-1]
# print(ss_clusters[1][:3])
# print(ss_seq, '\n')


# from Bio import SeqIO
# ref = next(SeqIO.parse('../data/covid_ref.fasta', 'fasta'))
# print(ref.seq[83: 126 + 1])

# regs = split_ss_part1(ss_seq, ss_clusters[1][1])
# print('||'.join([x['ss_seq'] for x in regs]), '\n')
# print(regs)

# print(ref.seq[83: 105])



        # if char != '.' and last_char == '.':
        #     if n_dots_before > max_free_space:
        #         end = pos - n_dots_before // 2
        #         cur_reg = ss_seq[start: end]
        #         if len(cur_reg) >= 10:
        #             regions.append(cur_reg)
        #             splitting_positions.append(start)
        #             start = end

        #     n_dots_before = 0
        # last_char = char
    
    # pos += 1
    # regions.append(ss_seq[start:])
    # splitting_positions.append(start)
    # splitting_positions.append(len(ss_seq))


# print([x[-1] for x in ss_clusters])

# split_ss2(ss_clusters[1][-1])

# for reg in ss_clusters:
#     # split_ss2(reg[-1])
#     print(reg[-1])




# print(max(list(map(len, [max(re.split('[()]', x[-1]), key=len) for x in ss_clusters]))))

# print_ss_fasta_to_draw(ss_clusters, 1)


# for x in ss_clusters_reliable:
#     print(f"{x[0]}\t0")
# break
