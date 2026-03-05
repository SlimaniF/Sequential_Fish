"""
This script aims at removing spots found in washout. 
If a spot is detected during a washout cycle, all spots detected in succeeding cycles at same location are deleted.
"""
from typing import cast
import pandas as pd
from ..settings import get_settings
from ..tools import safe_merge_no_duplicates
from ..customtypes import PipelineParameters

def main(run_path) :
    
    print(f"washout runing for {run_path}")
    
    pipeline_parameters = get_settings(run_path)
    pipeline_parameters = cast(PipelineParameters, pipeline_parameters)
    WASHOUT_KEY_WORD = pipeline_parameters.WASHOUT_KEY_WORD

    Acquisition = pd.read_feather(run_path + '/result_tables/Acquisition.feather')
    Detection = pd.read_feather(run_path + '/result_tables/Detection.feather')
    Clusters = pd.read_feather(run_path + '/result_tables/Clusters.feather')
    Spots = pd.read_feather(run_path + '/result_tables/Spots.feather')
    Gene_map = pd.read_feather(run_path + '/result_tables/Gene_map.feather')

    Gene_map['is_washout'] = Gene_map['target'].str.contains(WASHOUT_KEY_WORD)


    # Joins
    Detection = safe_merge_no_duplicates(
        Detection,
        Acquisition,
        keys= ['cycle'],
        on= 'acquisition_id'
    )
    Spots = safe_merge_no_duplicates(
        Spots,
        Detection,
        keys=['cycle','color_id'],
        on= 'detection_id'
    )
    print(Gene_map)

    print(Spots['cycle'].unique())
    print(Spots[Spots['color_id'].isna()])

    Spots = safe_merge_no_duplicates(
        Spots,
        Gene_map,
        keys='is_washout',
        on= ['cycle','color_id']
    )

    # putting coordinates in tuples
    Spots['coordinates'] = list(zip(Spots['z'], Spots['y'], Spots['x']))

    #Filtering
    cycle_list = list(Spots['cycle'].unique())
    cycle_list.sort()

    banned_coordinates = []
    for cycle in cycle_list :

        print(f"\nStarting cycle {cycle}")
        print(f"Currently {len(banned_coordinates)} coordinates are banned.")


        new_washout_idx = Spots.loc[
            (Spots['coordinates'].isin(banned_coordinates)) & (Spots['cycle'] == cycle)
            ].index

        Spots.loc[new_washout_idx, ['is_washout']] = True

        new_banned_coordinates = list(
            pd.unique(Spots[
                (Spots['is_washout']) & (Spots['cycle'] == cycle)
                ]['coordinates'])
        )

        #Updating banned coordinates
        ban_len = len(banned_coordinates)
        banned_coordinates += new_banned_coordinates
        banned_coordinates = list(pd.Series(banned_coordinates).unique())
        print(f"{len(banned_coordinates) - ban_len} coordinates were added to ban list.")

    #Merging with clusters and propagating washout
    Spots_groupby_cluster = Spots.groupby('cluster_id')['is_washout'].max().reset_index(drop=False)

    Clusters = safe_merge_no_duplicates(
        Clusters,
        Spots_groupby_cluster,
        keys='is_washout',
        on='cluster_id'
    )

    Spots.loc[Spots['cluster_id'].isin(Clusters[Clusters['is_washout'] == True])]['is_washout'] = True

    print("\nSaving results...")

    Spots.reset_index(drop=True).to_feather(run_path + "/result_tables/Spots.feather")
    Clusters.reset_index(drop=True).to_feather(run_path + "/result_tables/Clusters.feather")
    print("Done")
