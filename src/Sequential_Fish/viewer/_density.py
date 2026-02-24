import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import bigfish.detection as detection
from typing import Literal, Iterable

def multichannel_clustering(
        Acquisition : pd.DataFrame, 
        Detection : pd.DataFrame, 
        Spots : pd.DataFrame, 
        Gene_map : pd.DataFrame,
        voxel_size : tuple,
        cluster_radius : int =None,
        nb_min_spots = 4,
        no_filtering = False,
        )-> pd.DataFrame:
    """ 
    Perform a DBSCAN with cluster radius = max(voxel_size) by pulling spots from all channels/cycles together.
    Aggregates results in resulting pd.DataFrame with keys :

        Axis 0 (Index) :
        ----------------
            lvl0 : location
            lvl1 : general_cluster_id : int, unique identifier for multichannel clusters.
        
        Axis 1 (Columns) :
        ------------------
            ['single_molecule_number', 'unique_target_number', 'target_names', 'z','y','x']

    """

    if len(voxel_size) != 3 : raise ValueError("Only 3 dimensional analysis is supported, voxel_size has {0} dimensions".format(len(voxel_size)))

    master_spots = _merge_data_for_Spots_informations(Acquisition, Detection, Spots, Gene_map, no_filtering=no_filtering)
    clustered_spots = _run_DBSCAN(master_spots, voxel_size=voxel_size,cluster_radius= cluster_radius, min_number_cluster=nb_min_spots)
    master_spots = _merge_DBSCAN_results(master_spots, clustered_spots)
    multichannel_clusters = _aggregating_results(master_spots)

    return multichannel_clusters

def _merge_data_for_Spots_informations(
        Acquisition : pd.DataFrame, 
        Detection : pd.DataFrame, 
        Spots : pd.DataFrame, 
        Gene_map : pd.DataFrame,
        no_filtering : bool,
        ) :
    """
    return Master Spots table with added keys : 'acquisition_id', 'target', 'location'
    """

    check_len = len(Detection)
    Detection = pd.merge(
        Detection,
        Acquisition.loc[:,['acquisition_id', 'cycle', 'location']],
        on= 'acquisition_id',
        suffixes=('','_acquisition')
    )
    if no_filtering : assert len(Detection) == check_len

    Detection = pd.merge(
        Detection,
        Gene_map.loc[:,['cycle','color_id','target']],
        on= ['cycle','color_id']
    )
    if no_filtering : assert len(Detection) == check_len

    check_len = len(Spots)
    Spots = pd.merge(
        Spots,
        Detection.loc[:,['detection_id', 'acquisition_id', 'target', 'location']]
    )
    if no_filtering : assert len(Spots) == check_len

    return Spots

def _run_DBSCAN(Spots : pd.DataFrame, voxel_size : tuple, cluster_radius, min_number_cluster = 4) :
    """
    Return Spots DataFrame (duplicates dropped) with keys
        spots coordinates : 'z','y','x'; 
        general_cluster_id : 'cluster_id' (-1 if free); 
        fov : 'location';
        general_cluster_centroid : 'cluster_centroid_z', 'cluster_centroid_y', 'cluster_centroid_x'
    """
    
    if type(cluster_radius) == type(None) : cluster_radius = max(voxel_size)

    #Pooling data    
    spots_coordinates_per_fov = Spots.groupby(['location']).agg({
        'z' : list,
        'y' : list,
        'x' : list,
        })

    #DBSCAN
    Spots_clustered = pd.DataFrame(columns=['cluster_id','z','y','x', 'location'])
    for location in spots_coordinates_per_fov.index : 
        data_selec = spots_coordinates_per_fov.loc[location]

        spots = np.array(
            list(zip(data_selec['z'],data_selec['y'],data_selec['x'],)),
            dtype=int)

        clustered_spot, clusters = detection.detect_clusters(
            spots,
            voxel_size,
            cluster_radius,
            min_number_cluster,
        )

        z,y,x, cluster_id = zip(*clustered_spot)
        cluster_z, cluster_y, cluster_x, spot_number, cluster_index = zip(*clusters)

        new_Spots_clustered = pd.DataFrame({
                'cluster_id' : cluster_id,
                'z' : z,
                'y' : y,
                'x' : x,
                'location' : location,
            })

        new_Clusters = pd.DataFrame({
            'cluster_id' : cluster_index,
            'cluster_centroid_z' : cluster_z, 
            'cluster_centroid_y' : cluster_y, 
            'cluster_centroid_x' : cluster_x, 
        })

        check_len = len(new_Spots_clustered)
        new_Spots_clustered = pd.merge(
            new_Spots_clustered,
            new_Clusters,
            on='cluster_id',
            how='left',
            validate='m:1'
        )
        assert len(new_Spots_clustered) == check_len , "check for duplication"

        Spots_clustered = pd.concat([
            Spots_clustered,
            new_Spots_clustered,
        ], axis= 0)

    Spots_clustered = Spots_clustered.drop_duplicates()
    return Spots_clustered

def _merge_DBSCAN_results(Spots_main : pd.DataFrame, Spots_clustered : pd.DataFrame) :
    check_len = len(Spots_main) #No filtering/duplication expected.

    Spots_clustered = Spots_clustered.rename(columns={'cluster_id' : 'general_cluster_id'})

    Spots_main = pd.merge(
        Spots_main,
        Spots_clustered,
        on=['location','z','y','x'],
        validate= 'm:1'
    )

    assert len(Spots_main) == check_len, "Spots main : {0}\check_len : {1}\n{2}\n{3}".format(len(Spots_main), check_len, Spots_main, Spots_clustered)
    return Spots_main

def _aggregating_results(merged_spots : pd.DataFrame) :
    """
        'spot_id' : 'count',
        'target' : ['nunique','unique'],
        'cluster_centroid_z' : 'first',
        'cluster_centroid_y' : 'first',
        'cluster_centroid_x' : 'first',
    """
    Clustered_spots = merged_spots.loc[merged_spots['general_cluster_id'] != -1] #Filtering free spots

    multichannel_clusters = Clustered_spots.groupby(['location', 'general_cluster_id']).aggregate({
        'spot_id' : 'count',
        'target' : ['nunique','unique'],
        'cluster_centroid_z' : 'first',
        'cluster_centroid_y' : 'first',
        'cluster_centroid_x' : 'first',
    })

    COLUMNS_INDEX = ['single_molecule_number', 'unique_target_number', 'target_names', 'z','y','x']
    multichannel_clusters.columns = pd.Index(COLUMNS_INDEX)

    return multichannel_clusters

def spot_count_map(Acquisition : pd.DataFrame, 
        Detection : pd.DataFrame, 
        Spots : pd.DataFrame, 
        Gene_map : pd.DataFrame,
        no_filtering = False,
        ):

    """
    Return spots count map array.  
    """

    shape = np.max(list(Acquisition['fish_shape']), axis=0)[:3]
    voxel_size = Detection['voxel_size'].iat[0]
    master_spots = _merge_data_for_Spots_informations(Acquisition, Detection, Spots, Gene_map, no_filtering=no_filtering)

    single_counts = master_spots.groupby(['location','z','y','x'])['spot_id'].count().reset_index(level=0,drop=False)
    results = []
    locations = Acquisition['location'].unique()
    for location in  locations :
        data = single_counts.loc[single_counts['location'] == location]
        new_array = np.zeros(shape=shape, dtype=int)
        if len(data) > 0 :
            z = list(data.index.get_level_values(0))
            y = list(data.index.get_level_values(1))
            x = list(data.index.get_level_values(2))
            count = list(data['spot_id'])
            new_array[z,y,x] = count

        results.append(new_array)

    
    count_map = np.stack(results)
    return count_map