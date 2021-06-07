"""
sampling 
"""

from typing import List
import os
from random import randint
import re
import sys

import numpy as np
import pandas as pd

PATH_TO_REGIONS = '../data/ss_regions_for_IWT.csv'
PATH_TO_CONTROL_REGIONS_ROUGH = '../data/ss_control_regions_rough_for_IWT.csv'


# def _release_positions_from_regions(regions: pd.DataFrame):
#     """ create merged list of positions """
#     positions = []
#     for _, start, end in regions.values:
#         for pos in range(start, end):
#             positions.append(pos)
#     return positions


# def _prepare_positions(path_to_regions: str):
#     """ read data and release merged """
#     ss_regions = pd.read_csv(path_to_regions)
#     ss_regions_positions = _release_positions_from_regions(ss_regions)
#     assert ss_regions.ss_seq.apply(len).sum() == len(ss_regions_positions)

#     return ss_regions_positions


# def _sample_from_merged(positions, n_samples=100, rlenght=20):
#     """ sample from merged """
#     samples = []
#     for _ in range(n_samples):
#         pstart = randint(0, len(positions) - 1 - rlenght)
#         reg_pos = positions[pstart: pstart + rlenght]
#         assert len(reg_pos) == rlenght
#         samples.append(reg_pos)
#     return np.array(samples)


def samples_generator(regions: pd.DataFrame, rlenght=20) -> list:
    """ online true random samples generation (no one merged, there are shuffling)"""
    regions_values = regions.values
    while True:
        np.random.shuffle(regions_values)
        output_sample = []
        for _, start, end in regions_values:
            for pos in range(start, end):
                output_sample.append(pos)
                if len(output_sample) == rlenght:
                    yield output_sample
                    output_sample = []


def sample_from_shuffled(regions: pd.DataFrame, n_samples=100, rlenght=20) -> np.ndarray:
    sgen = samples_generator(regions, rlenght)
    samples = np.array([next(sgen) for _ in range(n_samples)])
    return samples


if __name__ == '__main__':
    # ss_reg_pos = prepare_positions(PATH_TO_REGIONS)
    # control_reg_pos = prepare_positions(PATH_TO_CONTROL_REGIONS_ROUGH)
    regions = pd.read_csv(PATH_TO_REGIONS)
    samples = sample_from_shuffled(regions, 10000, 20)
    print(regions)
    print()
    print(samples)
    print(samples.shape)

    
    # ss_samples = sample_from_merged(ss_reg_pos)
    
    # print(ss_samples)
    # print(ss_samples.shape)


# import matplotlib.pyplot as plt

# df = pd.read_csv(PATH_TO_REGIONS)
# print(df.ss_seq.apply(len).quantile(.75))
# # df.ss_seq.apply(len).plot.hist(bins=50)
# # plt.show()

# df = pd.read_csv(PATH_TO_CONTROL_REGIONS_ROUGH)
# print(df.ss_seq.apply(len).quantile(.75))
# # df.ss_seq.apply(len).plot.hist(bins=50)
# # plt.show()
