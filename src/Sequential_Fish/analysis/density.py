"""
Submodules aiming at de
"""

import os
import logging, traceback
import pandas as pd
import numpy as np
import bigfish.detection as detection
import matplotlib.pyplot as plt

from .utils import _get_min_cluster_radius, _get_voxel_size
from ..tools import safe_merge_no_duplicates

#properties hints
class multi_rna_cluster_DataFrame(pd.DataFrame):
    location : pd.Index
    general_cluster_id : pd.Index
    spot_number: pd.Series 
    rna_number: pd.Series 
    rna_list: pd.Series 
    z: pd.Series 
    y: pd.Series 
    x: pd.Series 

class updated_spots_DataFrame(pd.DataFrame) :
    spot_id : pd.Series
    z: pd.Series 
    y: pd.Series 
    x: pd.Series
    diversity : pd.Series

# Preprocessing

def merge_data_in_Spots(
    Acquisition : pd.DataFrame,
    Detection : pd.DataFrame,
    Spots : pd.DataFrame,
    Gene_map : pd.DataFrame
) :
    """
    Add "target", "location" keys to Spots.
    """
    
    Detection = safe_merge_no_duplicates(
        Detection,
        Acquisition,
        on='acquisition_id',
        keys= ['cycle']
    )
    
    Detection = safe_merge_no_duplicates(
        Detection,
        Gene_map,
        on= ['cycle', 'color_id'],
        keys= 'target'
    )
    
    Spots = safe_merge_no_duplicates(
        Spots,
        Detection,
        on= 'detection_id',
        keys= ["target", 'location']
    )
    
    return Spots

def group_coordinates_per_fov(Spots) :
    
    spots_coordinates_per_fov = Spots.groupby(['location']).agg({
    'z' : list,
    'y' : list,
    'x' : list,
    })
    
    return spots_coordinates_per_fov

def _coloc_clustered_spots(
    spots_coordinates_per_fov : pd.DataFrame,
    voxel_size : tuple,
    cluster_radius : int,
    min_spots_per_cluster : int,
    
) :
    """
    Creates dataframe giving a cluster_id to each pixel where spots where detected.
    
    -1 for unclustered pixel
    
    """
    
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
            min_spots_per_cluster,
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

def multi_rna_clusters(Spots) -> multi_rna_cluster_DataFrame:
    """
    Create df summing up clusters with location, single molecule number, number of different rna (different population) and list of all rna names in cluster.
    """
    
    Clustered_spots = Spots.loc[Spots['general_cluster_id'] != -1]
    
    multi_rna_clusters = Clustered_spots.groupby(['location', 'general_cluster_id']).aggregate({
        'spot_id' : 'count',
        'target' : ['nunique','unique'],
        'cluster_centroid_z' : 'first',
        'cluster_centroid_y' : 'first',
        'cluster_centroid_x' : 'first',
    })

    multi_rna_clusters.columns = pd.Index(['spot_number', 'rna_number', 'rna_list', 'z', 'y', 'x'])
    multi_rna_clusters = multi_rna_clusters.rename(columns={"rna_number" : "diversity"})
    return multi_rna_clusters

def update_and_filter_spots(Spots, multirna_clusters: multi_rna_cluster_DataFrame) -> updated_spots_DataFrame :
    """
    Add multi rna clusters information to single spots df
    """
    updated_spots = pd.merge(
        Spots,
        multirna_clusters,
        on=['location','general_cluster_id'],
        validate= 'm:1'
    )
    
    assert len(Spots[Spots['general_cluster_id'] != -1]) == len(updated_spots), "Non unique (location,general_cluster_id)"
    
    updated_spots = updated_spots.rename(columns={"rna_number" : "diversity"})

    return updated_spots

def create_cluster_unity_df(
    updated_spots,
    min_diversity = 3
) :
    """
    each line is a cluster and each RNA of data set are reprensented in columns with their single molecule number in this cluster; init to 0.
    """
    
    #Filtering spots from cluster with less than MIN_DIVERSITY different rnas.
    data = updated_spots.loc[updated_spots['diversity'] >= min_diversity]

    #Init empty table : 
    cluster_list = pd.Series(zip(data['location'], data['general_cluster_id'])).unique()
    cluster_list = pd.MultiIndex.from_tuples(cluster_list)
    RNA_list = data['target'].unique()
    cluster_unity_counts = pd.DataFrame(index=cluster_list, columns=RNA_list, dtype=int, data=0).sort_index()

    #Creating table with all counts of single molecule per cluster per rna
    agg_table = data.groupby(['location', 'general_cluster_id','target'])['spot_id'].count()
    for location, cluster_id, rna in agg_table.index :
        value = agg_table.at[(location,cluster_id,rna)]
        cluster_unity_counts.loc[(location,cluster_id), rna] += value

    #Assign values from groupby to cluster table
    total = cluster_unity_counts.sum(axis=1).astype(int)
    cluster_unity_counts['total'] = total
    
    return cluster_unity_counts

def create_presence_dict(
    cluster_unity_counts : pd.DataFrame, 
    RNA_list : list
    ) :
    presence_dict = {}
    for rna in RNA_list :
        rna_index = cluster_unity_counts.loc[cluster_unity_counts[rna] > 0].index
        data = cluster_unity_counts.loc[rna_index]
        presence_df = data.apply(lambda x : x.astype(bool))
        presence_dict[rna] = presence_df
        
    return presence_dict
        
def create_affinity_dict(
    cluster_unity_counts : pd.DataFrame, 
    RNA_list : list
    ) :
    affinity_dict = {}
    for rna in RNA_list :
        rna_index = cluster_unity_counts.loc[cluster_unity_counts[rna] > 0].index
        data = cluster_unity_counts.loc[rna_index]
        total = data['total']
        affinity_df = data.apply(lambda x : x/total.loc[rna_index])
        affinity_df = affinity_df.drop('total', axis=1)
        affinity_df_total = affinity_df.sum(axis=1).round(10)
        assert (affinity_df_total == 1).all(), affinity_df[affinity_df_total != 1]
        affinity_dict[rna] = affinity_df

    return affinity_dict
    
# Plots

def single_number_VS_cluster_diversity(
    multirna_clusters : multi_rna_cluster_DataFrame,
    cluster_radius,
    min_number_spots,
    output_path
) :
    fig = plt.figure()
    ax = fig.gca()

    data = multirna_clusters.groupby("diversity")['spot_number'].sum()
    ax.bar(data.index, data, align='center')
    ax.set_xticks(data.index)
    ax.set_ylabel("Total rna number")
    ax.set_xlabel("RNA diversity in cluster")
    ax.set_title(f"cluster radius : {cluster_radius} nm, min_spot : {min_number_spots}")
    
    plt.savefig(output_path + "/Single_number_VS_cluster_diversity.svg")
    plt.close()
    
def single_number_per_rna_per_diversity(
    updated_spots,
    cluster_radius,
    min_number_spots,
    output_path
) :
    cluster_plurality_number = len(updated_spots['diversity'].unique())
    
    fig, axes = plt.subplots(nrows= cluster_plurality_number, ncols=1, figsize = (16,8*cluster_plurality_number))

    for ax, dimension in zip(axes, range(1,cluster_plurality_number + 1)) :
        ax:plt.Axes
        ax.set_title(f"RNA diversity : >={dimension}; cluster radius : {cluster_radius} nm, min_spot : {min_number_spots}")

        data = updated_spots.loc[updated_spots['diversity'] >= dimension]
        data = data.groupby('target')['spot_id'].count()

        X = list(range(len(data)))
        ax.bar(X, data,align='center')
        ax.set_xticks(X, labels=data.index)
        ax.set_ylabel("Total rna number")

    plt.savefig(output_path + "/Single_number_per_rna_per_diversity.svg")
    plt.close()

def presence_plots(
    presence_dict : dict,
    RNA_list : list,
    cluster_radius : int,
    min_nb_cluster : int,
    output_path : str,
) :
    fig, axes = plt.subplots(nrows=len(RNA_list), ncols=1, figsize= (16,8*len(RNA_list)))
    
    for ax, rna in zip(axes, RNA_list) :
        ax : plt.Axes

        data : pd.DataFrame = presence_dict[rna]
        if 'total' in data.columns : data = data.drop(columns="total")
        cluster_number = len(data)
        data = data.sum(axis=0) / cluster_number * 100

        ax.set_title(f"Presence with {rna}; {cluster_number} clusters \nParameters : radius={cluster_radius} nm; min_spots={min_nb_cluster}")
        ax.bar(data.index, data)
        ax.set_ylabel("% presence")

    plt.savefig(output_path + "/presence_plots.svg")
    plt.close()
    
def affinity_plots(
    affinity_dict : dict,
    RNA_list : list,
    cluster_radius : int,
    min_nb_cluster : int,
    output_path : str,
) :
    fig, axes = plt.subplots(nrows=len(RNA_list), ncols=1, figsize= (16,8*len(RNA_list)))

    for ax, rna in zip(axes, RNA_list) :
        ax : plt.Axes

        data : pd.DataFrame = affinity_dict[rna]
        cluster_number = len(data)
        data = data.mean(axis=0)
        data_std = data.std(axis=0)

        ax.set_title(f"Affinity with {rna}; {cluster_number} clusters\nParameters : radius={cluster_radius} nm; min_spots={min_nb_cluster}")
        ax.bar(data.index, data)
        ax.set_ylabel("affinity")
    
    plt.savefig(output_path + "/affinity_plots.svg")
    plt.close()
   
def density_analysis(
    Acquisition : pd.DataFrame,
    Detection : pd.DataFrame,
    Spots : pd.DataFrame,
    Gene_map : pd.DataFrame,
    run_path : str,
    min_number_spots = 3,
    min_diversity = 3,
    cluster_radius = None,
) :
    """
    Main function for density analysis. All new sub function must be added to this function to be called within analysis pipeline.

    ANALYSIS PARAMETERS
    -------------------
        min_number_spots : int; 
        Minimum number of spots to consider when analysis colocalisation clusters.
    
        min_diversity : int, 
        Minimum number of rna population to consider for multi rna colocalisation clusters. In other words dimensionality of colocalization.
    
        cluster_radius : int
            Colocalization radius in **nanometers**. If None max of *voxel_size* is taken.
    """

   
    output_path = run_path + "/analysis/density_analysis/"
    os.makedirs(output_path, exist_ok=True)
    
    log_file = output_path + "/density_analysis_log.log"
    logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    force= True
)
    
    try : 
        print("Starting density analysis...")
        logging.info(f"New density analysis")
        logging.info(f"min_number_spots :\n{min_number_spots}\min_diversity :\n{min_diversity}\cluster_radius :\n{cluster_radius}\n")
        
        #Data processing
        Spots = merge_data_in_Spots(
            Acquisition=Acquisition,
            Detection=Detection,
            Spots=Spots,
            Gene_map=Gene_map,
        )

        voxel_size = _get_voxel_size(Detection)
        if cluster_radius is None :
            cluster_radius = _get_min_cluster_radius(voxel_size)

        spots_coordinates_per_fov = group_coordinates_per_fov(Spots)

        coloc_clustered_spots = _coloc_clustered_spots(
            spots_coordinates_per_fov,
            voxel_size=voxel_size,
            cluster_radius=cluster_radius,
            min_spots_per_cluster=min_number_spots
        )

        coloc_clustered_spots = coloc_clustered_spots.rename(columns={'cluster_id' : 'general_cluster_id'})

        Spots = safe_merge_no_duplicates(
            Spots,
            coloc_clustered_spots,
            on=['location','z','y','x'],
            keys= ['general_cluster_id', 'cluster_centroid_z', 'cluster_centroid_y', 'cluster_centroid_x']
        )

        coloc_clusters = multi_rna_clusters(Spots)
        updated_Spots = update_and_filter_spots(Spots, coloc_clusters)
        cluster_unity_df = create_cluster_unity_df(updated_Spots, min_diversity=min_diversity)

        RNA_list = list(Spots['target'].unique())
        presence_dict = create_presence_dict(cluster_unity_df, RNA_list)
        affinity_dict = create_affinity_dict(cluster_unity_df, RNA_list)

        #Plots
        single_number_VS_cluster_diversity(
            multirna_clusters=coloc_clusters,
            cluster_radius=cluster_radius,
            min_number_spots=min_number_spots,
            output_path=output_path
        )
        single_number_per_rna_per_diversity(
            updated_spots=updated_Spots,
            cluster_radius=cluster_radius,
            min_number_spots=min_number_spots,
            output_path=output_path
        )
        presence_plots(
            presence_dict=presence_dict,
            RNA_list=RNA_list,
            cluster_radius=cluster_radius,
            min_nb_cluster=min_number_spots,
            output_path=output_path
        )
        affinity_plots(
            affinity_dict=affinity_dict,
            RNA_list=RNA_list,
            cluster_radius=cluster_radius,
            min_nb_cluster=min_number_spots,
            output_path=output_path
        )
    
    except Exception as e :
        logging.error(f"analysis failed :\n{traceback.format_exc()}")
        return False
        
    else :
        logging.info(f"analysis succeed")
        
        return True