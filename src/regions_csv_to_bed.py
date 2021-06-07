import os
import sys
import json
import warnings

import numpy as np
import pandas as pd

from add_features_to_rerions import build_region_feature

PATH_TO_REGIONS = '../data/ss_regions_for_IWT.csv'
# PATH_TO_CONTROL_REGIONS = '../data/ss_control_regions_for_IWT.csv'
PATH_TO_CONTROL_REGIONS_ROUGH = '../data/ss_control_regions_rough_for_IWT.csv'


def write_elements(features_path, regions_path, featurename, regions):
    assert featurename in regions.columns
    chr_id = 'chr1'
    feature_content = []
    bed_content = []
    cols = ['start', 'end', featurename]
    for start, end, feavals in regions[cols].values:
        features = "\t".join(map(str, feavals))
        bed_row = '\t'.join([chr_id, str(start), str(end)])
        feature_row = '\t'.join([chr_id, str(start), str(end), features])
        feature_content.append(feature_row)
        bed_content.append(bed_row)

    with open(features_path, 'w') as fout:
        fout.write('\n'.join(feature_content))

    with open(regions_path, 'w') as fout:
        fout.write('\n'.join(bed_content))


def write_dataset(data_path):
    dataset_path = os.path.join(data_path, 'datasets.txt')
    dataset_content = [
        "id\tname\tregionFile", 
        'elem1\tElements 1\tElements_regions.bed',
        'control\tControls\tControls_regions.bed'
    ]
    with open(dataset_path, 'w') as fout:
        fout.write('\n'.join(dataset_content))


def write_features_datasetsTable(regions, control_regions, data_path):
    print('writing features', file=sys.stderr)
    write_dataset(data_path)

    features_datasetsTable_path = os.path.join(
        data_path, "features_datasetsTable.txt")
    features_datasets_columns = "id	name	elem1	control"
    feature_table_content = [features_datasets_columns]

    feature_names = [x for x in regions.columns if x not in ['start', 'end', 'ss_seq']]
    for i, feaname in enumerate(feature_names):
        idx = f'ftr{i}'
        elements_path = f'{feaname}_elements.txt'
        elements_regions_path = 'Elements_regions.bed'
        control_path = f'{feaname}_control.txt'
        control_regions_path = 'Controls_regions.bed'
        row = '\t'.join([idx, feaname, elements_path, control_path])
        feature_table_content.append(row)
        
        write_elements(os.path.join(data_path, 'files', elements_path), 
                      os.path.join(data_path, 'files', elements_regions_path), 
                      feaname, regions)

        write_elements(os.path.join(data_path, 'files', control_path), 
                       os.path.join(data_path, 'files', control_regions_path), 
                       feaname, control_regions)

    with open(features_datasetsTable_path, 'w') as fout:
        fout.write('\n'.join(feature_table_content))


if __name__ == '__main__':
    warnings.filterwarnings('ignore')

    window_size = 10
    reglen = 20
    n_samples = 10000
    max_distance = 500
    save_positions = False
    edge_level = 4  # grandpa
    filter_mode = 'more'
    output_path = f'../data/__iwt_dataset{reglen}_{window_size}_{n_samples}_{max_distance}_lvl-{edge_level}_mode-{filter_mode}/'

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
