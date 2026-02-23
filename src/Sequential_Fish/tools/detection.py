import pandas as pd
import numpy as np

def build_Spots_and_Cluster_df(detection_result : dict) :
    """
    Made to build Spots pandas dataframe from result of multi_threaded call to `multi_thread_full_detection`.

    Parameter
    ---------
        detection_result
    """

    SPOTS_COLUMNS = [
        'detection_id',
        'spots',
        'spots_post_decomp',
        'clustered_spots_dataframe',
        'clusters_dataframe',
        'clusters',
        'clustered_spots',
        'free_spots',   
    ]
    
    detection_res = pd.DataFrame(columns=pd.Index(SPOTS_COLUMNS), data=detection_result)
    detection_res = detection_res.set_index('detection_id', verify_integrity=True, drop=False)

    Spots = pd.DataFrame()
    Clusters = pd.DataFrame()

    for detection_id in detection_res.index :
        spots: pd.DataFrame = detection_res.at[detection_id, 'clustered_spots_dataframe']
        clusters: pd.DataFrame = detection_res.at[detection_id, 'clusters_dataframe']

        if isinstance(spots, (int,float)) : 
            spots = pd.DataFrame(columns=pd.Index(['id','cluster_id','population', 'detection_id'])) 
        if isinstance(clusters, (int,float)) : 
            clusters = pd.DataFrame(columns=pd.Index(['id','population', 'detection_id']))
        
        #names
        spots = spots.rename(columns={'id' : 'spot_id'})
        clusters = clusters.rename(columns={'id' : 'cluster_id'})
        
        #is_clustered
        spots['cluster_id'] = spots['cluster_id'].replace(-1,np.nan)
        clusters['cluster_id'] = clusters['cluster_id'].replace(-1,np.nan)
        spots['population'] = spots['cluster_id'].isna()
        spots['population'] = spots['population'].replace({True : 'free', False : 'clustered'})

        #Detection_id
        spots['detection_id'] = detection_id
        clusters['detection_id'] = detection_id

        Spots = pd.concat([
            Spots,
            spots,
        ],axis=0)

        Clusters = pd.concat([
            Clusters,
            clusters,
        ],axis=0)

    Spots = Spots.reset_index(drop=True)
    Clusters = Clusters.reset_index(drop=True)

    return Spots, Clusters