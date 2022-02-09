"""
нужно посчитать фичи на датасетах с различным уровнем фильтрации 
претерминальных узлов. 
удаляться будут постепенно от уровня 1 (все), уровень 3 (без бабушек, родителей и терминалей) 
до 12 уровня, в котором только замещения из корня
"""

import os
import warnings

from add_features_to_rerions import build_region_feature
from regions_csv_to_bed import (
    write_features_datasetsTable,
    PATH_TO_REGIONS,
    PATH_TO_CONTROL_REGIONS_ROUGH,
)

PATH_TO_DISTANCES = '../data/new_final_fantasy_Sat_Jun__5_20:28:29_2021.csv'


if __name__ == '__main__':
    warnings.filterwarnings('ignore')

    window_size = 10
    reglen = 20
    n_samples = 10000
    max_distance = 30000
    save_positions = False
    filter_mode = 'more'

    for edge_level in range(1, 10):
        print(f'edge level: {edge_level}')
        output_path = f'../data/_iwt_dataset{reglen}_{window_size}_{n_samples}_{max_distance}_lvl-{edge_level}_mode-{filter_mode}/'

        # count features
        regions = build_region_feature(
            path_to_regions=PATH_TO_REGIONS,
            only_nsubst=False, 
            window_size=window_size, 
            region_len=reglen,
            n_samples=n_samples,
            max_distance=max_distance,
            save_positions=save_positions,
            edge_level=edge_level,
            filter_mode=filter_mode,
        )
        control_regions = build_region_feature(
            path_to_regions=PATH_TO_CONTROL_REGIONS_ROUGH, 
            only_nsubst=False, 
            window_size=window_size, 
            region_len=reglen,
            n_samples=n_samples,
            max_distance=max_distance,
            save_positions=save_positions,
            edge_level=edge_level,
            filter_mode=filter_mode,
        )
        # prepare directory
        if not os.path.exists(output_path):
            os.mkdir(output_path)
            os.mkdir(os.path.join(output_path, 'files'), )

        write_features_datasetsTable(regions, control_regions, output_path)