"""
This script aims at removing spots found in washout. 
If a spot is detected during a washout cycle, all spots detected in succeeding cycles at same location are deleted.
"""
import warnings, sys
import pandas as pd

from Sequential_Fish.tools import safe_merge_no_duplicates

def main(run_path) :
    
    print(f"washout runing for {run_path}")
    
    if len(sys.argv) == 1:
        from Sequential_Fish.pipeline_parameters import WASHOUT_KEY_WORD
    else :
        from Sequential_Fish.run_saves import get_parameter_dict
        PARAMETERS = ['WASHOUT_KEY_WORD']
        parameters_dict = get_parameter_dict(run_path, parameters=PARAMETERS)
        WASHOUT_KEY_WORD = parameters_dict['WASHOUT_KEY_WORD']

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
            Spots[
                (Spots['is_washout']) & (Spots['cycle'] == cycle)
                ]['coordinates'].unique()
        )

        #Updating banned coordinates
        ban_len = len(banned_coordinates)
        banned_coordinates += new_banned_coordinates
        banned_coordinates = list(pd.unique(banned_coordinates))
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
    print(Spots)


    print("\nSaving results...")

    Spots.reset_index(drop=True).to_feather(run_path + "/result_tables/Spots.feather")
    Clusters.reset_index(drop=True).to_feather(run_path + "/result_tables/Clusters.feather")
    print("Done")

if __name__ == "__main__":
    if len(sys.argv) == 1:
        warnings.warn("Prefer launching this script with command : 'python -m Sequential_Fish pipeline input' or make sure there is no conflict for parameters loading in pipeline_parameters.py")
        from Sequential_Fish.pipeline_parameters import RUN_PATH as run_path
    else :
        run_path = sys.argv[1]
    main(run_path)        
