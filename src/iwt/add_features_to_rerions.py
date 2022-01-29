"""
на вход: csv с регионами
на выход: csv с фичами

classes for fetures construction:
FeatureConstructor
SamplesFeatureConstructor

"""
from functools import partial
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix

from sample_from_merged_regions import sample_from_shuffled

TREE_PATH = "../data/mulal_gisaid_2021-01-22.filtered.twice.fasta.prank.anc.dnd"
PATH_TO_SUBSTITUTIONS = '../data/overall_mutations_with_context2.csv'
# PATH_TO_DISTANCES = '../data/new_final_fantasy_Sat_May__8_13:54:58_2021.csv'
# PATH_TO_DISTANCES = '../data/new_final_fantasy_Sat_Jun__5_20:28:29_2021.csv'
PATH_TO_DISTANCES = '../data/distances_06-06-21.csv'
PATH_TO_REGIONS = '../data/ss_regions_for_IWT.csv'


class FeatureConstructor:
    def __init__(self,
                 regions: pd.DataFrame,
                 substitutions: pd.DataFrame,
                 distances: pd.DataFrame,
                 window=None):
        self.regions = regions.copy()
        self.substitutions = substitutions.copy()

        self.distances = distances.copy()
        self.window = window or self.determine_window()
        print(f"window size: {self.window}", file=sys.stderr)
        self.nsubst = self.substitutions.groupby(
            'pos').apply(len).reset_index()
        self.nsubst.columns = ['pos', 'count']

    def determine_window(self, plot=False):
        sizes = (self.regions.end - self.regions.start + 1).values
        median = int(np.median(sizes))
        if plot:
            plt.hist(sizes, bins=30)
            plt.show()
        return median

    @staticmethod
    def region_split(start, end: int, window: int):
        """ split region to intervals with given window """
        intervals = []
        for i in range(start, end + 1, window):
            intervals.append([i, i + window - 1])
        intervals[-1][1] = end  # костыль
        return intervals

    def feature_base(self, feature: callable):
        """ calculate statistics along regions by given callable `feature` """
        region_stats = []
        for reg in self.regions[['start', 'end']].values:
            start, end = reg[1], reg[2]
            intervals = self.region_split(start, end, self.window)
            istats = [feature(istart, iend) for istart, iend in intervals]
            region_stats.append(istats)
        return pd.Series(region_stats)

    def count_mean_nsubst_in_interval(self, istart: int, iend: int):
        """ count mean number of substitutions per position in interval

        @param istart: interval start (including)
        @param iend: interval end (including)
        """
        nsubst = self.nsubst
        subst_count = nsubst[(nsubst.pos >= istart) &
                             (nsubst.pos <= iend)]['count'].mean()
        return subst_count

    def subst_count_feature(self, add_to_table=False):
        feature_values = self.feature_base(self.count_mean_nsubst_in_interval)
        if add_to_table:
            self.regions['subst_count'] = feature_values
        return feature_values

    def count_mean_distance_to_closest_mut(
            self, istart: int, iend: int, primary: bool):
        """ return mean distance to closest substitution

        @param primary: if True primary rna distance, else secondary distance
        """
        distances = self.distances[self.distances.primary_dist2nearest < 30_000]
        subst_count = distances[(distances.pos >= istart) &
                                (distances.pos <= iend)]
        if primary:
            dist = subst_count.primary_dist2nearest.mean()
        else:
            dist = subst_count.secondary_dist2nearest.mean()
        return dist

    def distance_feature(self, primary=True, add_to_table=False):
        feature_values = self.feature_base(
            partial(self.count_mean_distance_to_closest_mut, primary=primary)
        )
        if add_to_table:
            col_name = f"distance{int(not primary) + 1}"
            self.regions[col_name] = feature_values
        return feature_values


class SamplesFeatureConstructor:
    def __init__(self,
                 regions: np.array,
                 substitutions: pd.DataFrame,
                 distances: pd.DataFrame,
                 window: int,
                 reg_df=None):
        self.regions = regions.copy()
        self.regions_df = reg_df if reg_df is not None else pd.DataFrame({
            'start': np.arange(0, np.prod(regions.shape), regions.shape[1]),
            'end': np.arange(
                regions.shape[1],
                np.prod(regions.shape) + regions.shape[1],
                regions.shape[1]
            ),
        })
        self.substitutions = substitutions.copy()
        self.distances = distances.copy()
        # assert regions.shape[1] % window == 0, 'region len % window != 0'
        self.window = window
        print(f"window size: {self.window}", file=sys.stderr)
        print(f"distances subst size: {distances.shape[0]}", file=sys.stderr)

    
    def calculate_positional_nsubst(self, substitutions=None):
        if substitutions is None:
            substitutions = self.substitutions
        nsubst = substitutions.groupby('pos').apply(len).reset_index()
        nsubst.columns = ['pos', 'count']
        return nsubst

    @staticmethod
    def region_split(reg_pos: np.array, window: int):
        """ split region to intervals with given window """
        intervals = []
        for i in range(0, len(reg_pos), window):
            intervals.append(reg_pos[i: i + window])
        return intervals

    def feature_base(self, feature: callable):
        """ calculate statistics along regions by given callable `feature` """
        region_stats = []
        for reg_pos in self.regions:
            intervals = self.region_split(reg_pos, self.window)
            istats = list(map(feature, intervals))
            region_stats.append(istats)
        return pd.Series(region_stats)

    # def _count_mean_nsubst_in_interval(self, interval_pos: np.array):
    #     """ count mean number of substitutions per position in interval

    #     @param istart: interval start (including)
    #     @param iend: interval end (including)
    #     """
    #     nsubst = self.nsubst
    #     subst_count = nsubst[nsubst.pos.isin(interval_pos)]['count'].mean()
    #     return subst_count

    # def _subst_count_feature(self, add_to_table=False):
    #     self.nsubst = self.calculate_positional_nsubst()
    #     feature_values = self.feature_base(self._count_mean_nsubst_in_interval)
    #     if add_to_table:
    #         self.regions_df['subst_count'] = feature_values
    #     return feature_values

    def count_mean_nsubst_in_interval(self, interval_pos: np.array, pos2nsubst):
        """ count mean number of substitutions per position in interval

        @param 
        """
        subst_count = [pos2nsubst[pos] for pos in interval_pos if pos in pos2nsubst]
        if len(subst_count) == 0:
            return np.nan
        return np.mean(subst_count)

    def subst_count_feature(self, substitutions=None, fname='', add_to_table=False):
        # calculate nsubst on filtered substitutions
        nsubst = self.calculate_positional_nsubst(substitutions)
        pos2nsubst = dict(nsubst.values)
        feature_values = self.feature_base(partial(
            self.count_mean_nsubst_in_interval,
            pos2nsubst=pos2nsubst,
        ))
        if add_to_table:
            self.regions_df[f'subst_count_{fname}'] = feature_values
        return feature_values

    def _count_mean_distance_to_closest_mut(
        self,
        interval_pos: np.ndarray,
        pos2dist: dict
    ):
        """ return mean distance to closest substitution in the interval

        first, calculate mean distance on position - incorrect decision !!!

        @param primary: if True primary rna distance, else secondary distance
        """
        idist = [pos2dist[pos] for pos in interval_pos if pos in pos2dist]
        if len(idist) == 0:
            return np.nan
        return np.mean(idist)

    def count_around_substs(self, interval_pos: np.ndarray, 
                            pos2n_around_subst: dict) -> np.float:
        values = [pos2n_around_subst[pos] for pos in interval_pos \
                                            if pos in pos2n_around_subst]
        return np.mean(values)

    def around_subst_feature(self, add_to_table=False, radius=50):
        substitutions = self.substitutions.copy()
        substitutions['pair'] = substitutions.parent_node + '_>_' + substitutions.child_node
        local_substs = substitutions.groupby(['pair', 'pos']).child_node.count().reset_index()
        pair2pos = dict(local_substs.groupby('pair').pos.apply(np.array).reset_index().values)
        n_around_subst = local_substs.apply(
            lambda x: np.sum(np.abs(x.pos - pair2pos[x.pair]) <= radius) - 1, axis=1)
        local_substs['n_around_subst'] = n_around_subst
        pos2n_around_subst = local_substs.groupby('pos').n_around_subst.mean()
        
        feature_values = self.feature_base(
            partial(
                self.count_around_substs,
                pos2n_around_subst=pos2n_around_subst,
            )
        )
        if add_to_table:
            col_name = f"n_subst_around_{radius}nt"
            self.regions_df[col_name] = feature_values
        return feature_values
    
    def count_mean_distance_to_closest_mut(
        self, interval_pos: np.ndarray, pos2dist: dict):
        """ calculate true interval mean distance to closest subst """
        idist_collection = [pos2dist[pos] for pos in interval_pos if pos in pos2dist]  # list of lists
        if len(idist_collection) == 0:
            return np.nan
        idist = [d for pos_dists in idist_collection for d in pos_dists]
        if len(idist) == 0:
            return np.nan
        return np.mean(idist)

    def default_distance_feature(self, primary=True, add_to_table=False):
        _name = 'primary' if primary else 'secondary'
        fea_name = f'{_name}_dist2nearest'
        # dists = self.distances.groupby('pos')[fea_name].mean().reset_index()
        dists = self.distances.groupby('pos')[fea_name].apply(list)
        pos2dist = dict(dists)  # dict(pos: [d1, d2, d3, ...])
        
        feature_values = self.feature_base(
            partial(
                self.count_mean_distance_to_closest_mut,
                pos2dist=pos2dist,
            )
        )
        if add_to_table:
            col_name = f"distance{int(not primary) + 1}"
            self.regions_df[col_name] = feature_values
        return feature_values

    def count_substs_strenght(self, interval_pos):
        cur_substs = self.position_strenght[
            self.position_strenght['pos'].isin(set(interval_pos))
        ]
        return cur_substs['strenght'].mean()

    def calculate_substitution_strenght(self):
        self.add_strenght_to_dataset(self.substitutions)
        self.position_strenght = self.substitutions.groupby(
            'pos')['strenght'].mean().reset_index()    

    @staticmethod
    def add_strenght_to_dataset(data):
        """ calculate substitution strenght 
        (-1 if number of bounds decreased else 1) 

        add column 'strenght' to data
        """
        data['strenght'] = 0
        # data['strenght'][
        #     (data.parent_nucl == 'A') & (data.child_nucl == 'G') | 
        #     (data.parent_nucl == 'T') & (data.child_nucl == 'C')
        # ] = 1
        # data['strenght'][
        #     (data.parent_nucl == 'G') & (data.child_nucl == 'A') | 
        #     (data.parent_nucl == 'C') & (data.child_nucl == 'T')
        # ] = -1
        data['strenght'][
            (data.parent_nucl.isin({'A', 'T'})) &
            (data.child_nucl.isin({'G', 'C'}))
        ] = 1
        data['strenght'][
            (data.parent_nucl.isin({'G', 'C'})) &
            (data.child_nucl.isin({'A', 'T'}))
        ] = -1

    def subst_strenght_feature(self, add_to_table=False):
        self.calculate_substitution_strenght()
        assert 'strenght' in self.substitutions.columns

        feature_values = self.feature_base(self.count_substs_strenght)
        if add_to_table:
            self.regions_df["strenght"] = feature_values
        return feature_values

    def positional_distance_feature(self, distances, strenght_name, primary=True, add_to_table=False):
        """ добавить позиции, на которых считать фичу, или условия """
        _name = 'primary' if primary else 'secondary'
        fea_name = f'{_name}_dist2nearest'

        # dists = distances.groupby('pos')[fea_name].mean()
        dists = distances.groupby('pos')[fea_name].apply(list)
        pos2dist = dict(dists)
        feature_values = self.feature_base(
            partial(
                self.count_mean_distance_to_closest_mut,
                pos2dist=pos2dist,
            )
        )
        if add_to_table:
            col_name = f"distance{int(not primary) + 1}_{strenght_name}"
            self.regions_df[col_name] = feature_values
        return feature_values
    
    def strenght_ncount_features(self, add_to_table=False):
        """ count distance features with info about substitutions strenght """
        if 'strenght' not in self.substitutions.columns:
            self.add_strenght_to_dataset(self.substitutions)
        
        # substitutions that becoming more hard or weak
        hard_substs = self.substitutions[self.substitutions.strenght == 1]
        weak_substs = self.substitutions[self.substitutions.strenght == -1]

        # count nsubst
        self.subst_count_feature(hard_substs, 'hard', add_to_table)
        self.subst_count_feature(weak_substs, 'weak', add_to_table)
        

    def strenght_distance_features(self, add_to_table=False):
        """ count distance features with info about substitutions strenght """
        # useless features
        if 'strenght' not in self.distances.columns:
            self.add_strenght_to_dataset(self.distances)

        # substitutions that becoming more hard or weak
        hard_dists = self.distances[self.distances.strenght == 1]
        weak_dists = self.distances[self.distances.strenght == -1]

        # count for primary
        self.positional_distance_feature(hard_dists, 'hard', True, add_to_table)
        self.positional_distance_feature(weak_dists, 'weak', True, add_to_table)

        # count for secondary
        self.positional_distance_feature(hard_dists, 'hard', False, add_to_table)
        self.positional_distance_feature(weak_dists, 'weak', False, add_to_table)
        

# print(SamplesFeatureConstructor.region_split(np.random.randint(0, 50, 20), 5))


def add_edge_level_to_table(substitutions: pd.DataFrame):
    """ запускать, только если нужно пересчитать уровни для замещений """
    from ete3 import PhyloTree
    from utils import node_parent

    tree = PhyloTree(TREE_PATH, format=1)
    
    edge_levels = dict()

    for leaf in tree.iter_leaves():
        cur_level = 1
        parent = leaf    
        while parent != tree:
            parent = node_parent(leaf)
            pair = (parent.name, leaf.name)
            edge_levels[pair] = cur_level

            leaf = parent
            cur_level += 1

    assert len(edge_levels) == len(tree.get_descendants())  # E = V - 1
    
    substitutions = substitutions.copy()
    substitutions['edge_level'] = substitutions.apply(
        lambda r: edge_levels[(r.parent_node, r.child_node)], axis=1)
    
    return substitutions


def sample_distances_by_edge_levels(substitutions: pd.DataFrame, level=1, mode='more'):
    """ Отбор лежит в глубине дерева. Мутагенез на терминалях. 
    Будем отбирать замещения начиная с определенного уровня, который назначен каждому ребру. 
    `level = 1` означает брать все замещения; 
    `level = 2` означает удаление замещений из терминальных узлов
    ... максимальное значение - `level = 12`

    @param mode (at each mode level is including):
        - 'more' - means more than level value (default)
        - 'equal' - equal to level value
        - 'less' - less than level value
    """
    if 'edge_level' not in substitutions.columns:
        print('adding edge level to table', file=sys.stderr)
        substitutions = add_edge_level_to_table(substitutions)
    
    if mode == 'more':
        out = substitutions[substitutions['edge_level'] >= level]
    elif mode == 'equal':
        out = substitutions[substitutions['edge_level'] == level]
    elif mode == 'less':
        out = substitutions[substitutions['edge_level'] <= level]
    else:
        Exception('Incorrect `mode` value')
    return out


def filter_input_substitutions(substitutions, distances, max_distance: int, 
                               edge_level: int, filter_mode: str):
    # only substitutions "from root" to catch selection
    # на простых замещениях (без расстояний) обойдемся старой фильтрацией, они не используются сильно
    substitutions = substitutions[substitutions.child_node.str.startswith('#')]
    distances = distances[distances.primary_dist2nearest < max_distance]

    distances = sample_distances_by_edge_levels(distances, edge_level, filter_mode)

    ## old method
    # terminal_nodes = distances[~distances.child_node.str.startswith('#')].child_node.unique()
    # pre_terminal_nodes = distances[~distances.child_node.str.startswith('#')].parent_node.unique()
    # pre_pre_terminal_nodes = distances[~distances.child_node.isin(set(pre_terminal_nodes))].parent_node.unique()

    # distances = distances[~distances.child_node.isin(terminal_nodes)]
    # distances = distances[~distances.child_node.isin(pre_terminal_nodes)]
    # distances = distances[~distances.child_node.isin(pre_pre_terminal_nodes)]

    return substitutions, distances


def build_region_feature(path_to_regions=PATH_TO_REGIONS,
                         path_to_subst=PATH_TO_SUBSTITUTIONS,
                         path_to_dist=PATH_TO_DISTANCES,
                         only_nsubst=True,
                         window_size=10,
                         region_len=20,
                         n_samples=100,
                         max_distance=500,
                         save_positions=False,
                         edge_level=4,
                         filter_mode='more'):
    # если не сохраняем позиции, то загружаем их из предварительно сохраненных
    # ----------------------
    fname = '../data/regions_positions.npy'
    if 'control' in path_to_regions:
        fname = fname.replace('regions', 'control_regions')

    if save_positions:
        raw_regions = pd.read_csv(path_to_regions)
        sample_regions = sample_from_shuffled(raw_regions, n_samples, region_len)
        np.save(fname, sample_regions)
    else:
        try:
            sample_regions = np.load(fname)
        except:
            Exception(f'error with loading sample_regions from {fname}')
    # ----------------------

    substitutions = pd.read_csv(path_to_subst)
    final_fantasy = pd.read_csv(path_to_dist)

    substitutions, distances = filter_input_substitutions(
        substitutions, final_fantasy, max_distance, edge_level, filter_mode,
    )

    constructor = SamplesFeatureConstructor(
        sample_regions, substitutions, distances, window=window_size)

    constructor.subst_count_feature(add_to_table=True)
    if not only_nsubst:
        constructor.default_distance_feature(primary=True, add_to_table=True)
        # constructor.default_distance_feature(primary=False, add_to_table=True)
        # constructor.subst_strenght_feature(add_to_table=True)
        # ------------------
        # constructor.around_subst_feature(add_to_table=True)
        # constructor.strenght_distance_features(add_to_table=True)
        # constructor.strenght_ncount_features(add_to_table=True)  # no need
    return constructor.regions_df


def main():
    regions = build_region_feature(only_nsubst=False)
    print(regions.head(2))

    # dists = pd.read_csv(PATH_TO_DISTANCES)
    # dists = add_edge_level_to_table(dists)
    # dists.to_csv('../data/distances_06-06-21.csv')
    # from IPython import embed; embed()


if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    main()

# substitutions = pd.read_csv(PATH_TO_SUBSTITUTIONS)
# substitutions['mut'] = substitutions.parent_nucl + substitutions.child_nucl

# substitutions['strenght'] = 0
# substitutions['strenght'][
#     (substitutions.parent_nucl.isin({'A', 'T'})) &
#     (substitutions.child_nucl.isin({'G', 'C'}))
# ] = 1
# substitutions['strenght'][
#     (substitutions.parent_nucl.isin({'G', 'C'})) &
#     (substitutions.child_nucl.isin({'A', 'T'}))
# ] = -1
# print(substitutions)
# position_strenght = substitutions.groupby('pos')['strenght'].mean()
# print(position_strenght)
# gr = substitutions.groupby('pos')
# print(gr.get_group(10)['mut'])


# final_fantasy = pd.read_csv(PATH_TO_DISTANCES)
# final_fantasy = final_fantasy[final_fantasy.primary_dist2nearest < 30000]


# print(pos2dist1)

# print(final_fantasy.pos.nunique())
# print(final_fantasy.shape)
# dme = final_fantasy.groupby('pos')['primary_dist2nearest'].mean().reset_index()
# pos2dist1 = dict(dme.values[:20])
# print(dme)
# print(dict(pos2dist1))

# print(pos2dist1[3])

# print(np.mean([]))
