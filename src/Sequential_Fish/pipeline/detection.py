
"""
Spot quantification for Sequential Fish data

This script use results from FishSeq_pipeline_segmentation.py that must be run before
"""

import os
import logging
from typing import cast

import numpy as np
import pandas as pd
from tqdm import tqdm
from pebble import ProcessPool
from smfishtools.detection import multi_thread_full_detection

from ..tools import open_location
from ..settings import get_settings
from ..customtypes.parameters import PipelineParameters

def main(run_path) :

    print(f"detection runing for {run_path}")
    
    pipeline_parameters = get_settings(run_path)
    pipeline_parameters = cast(PipelineParameters, pipeline_parameters)
    
    VOXEL_SIZE = pipeline_parameters.VOXEL_SIZE
    SPOT_SIZE = pipeline_parameters.SPOT_SIZE
    ALPHA = pipeline_parameters.ALPHA
    BETA = pipeline_parameters.BETA
    GAMMA = pipeline_parameters.GAMMA
    CLUSTER_SIZE = pipeline_parameters.CLUSTER_SIZE
    MIN_SPOT_PER_CLUSTER = pipeline_parameters.MIN_SPOT_PER_CLUSTER
    ARTIFACT_RADIUS = pipeline_parameters.ARTIFACT_RADIUS
    DETECTION_SLICE_TO_REMOVE = pipeline_parameters.DETECTION_SLICE_TO_REMOVE
    MAX_WORKERS = pipeline_parameters.detection_MAX_WORKERS
    WAVELENGTH_LIST = pipeline_parameters.WAVELENGTH_LIST
    
    #Loading data
    Acquisition = pd.read_feather(run_path + "/result_tables/Acquisition.feather")

    #preparing folders
    save_path = run_path + "/result_tables/"

    #Main loop
    Detection_save = pd.DataFrame()
    Spots_save = pd.DataFrame()
    Clusters_save = pd.DataFrame()

    os.makedirs(run_path + '/detection_fov/',exist_ok=True)

    max_id = 0
    Detection = pd.DataFrame()
    for location_id, location in enumerate(Acquisition['location'].unique()) :
        print("Starting location {0}...".format(location_id))
        sub_data = Acquisition.loc[Acquisition["location"] == location]

        visual_path = run_path + "/visuals/{0}/".format(location)
        os.makedirs(visual_path, exist_ok=True)

        #Opening images
        multichannel_stack = open_location(Acquisition, location)
        
        dapi_channel = sub_data['dapi_channel'].iat[0]
        bead_channel = sub_data['bead_channel'].iat[0]
        if not bead_channel is None : 
            end_signal = min(dapi_channel, bead_channel)
        else :
            end_signal = dapi_channel

        #Converting na back to None
        bottom_index, top_index = DETECTION_SLICE_TO_REMOVE
        if bottom_index is None :
            pass
        elif not isinstance(bottom_index,int) : 
            bottom_index= int(bottom_index)
        
        if top_index is None :
            pass
        elif np.isnan(top_index) : 
            top_index = None
        elif type(top_index) != int : 
            top_index= int(top_index)
        
        
        #Removing Z slices (USER SETTING)
        if not top_index is None : top_index = -top_index
        
        multichannel_stack = multichannel_stack[:,bottom_index:top_index]
        multichannel_stack = multichannel_stack[...,:end_signal]
        images_list = [np.moveaxis(channel,[3,0,1,2],[0,1,2,3]) for channel in multichannel_stack]
        images_list = [
            [colors for colors in channel]
             for channel in images_list]
        image_number = len(multichannel_stack)
        colors = list(zip(*images_list))

        #Preparing threads arguments
        Detection = pd.DataFrame({
            'acquisition_id' : list(sub_data['acquisition_id'])
            ,'filename' : list(sub_data['full_path'])
            ,'voxel_size' : [tuple(VOXEL_SIZE)] * image_number
            ,'spot_size' : [tuple(SPOT_SIZE)] * image_number
            ,'alpha' : [ALPHA] * image_number
            ,'beta' : [BETA] * image_number
            ,'gamma' : [GAMMA] * image_number
            ,'artifact_radius' : [ARTIFACT_RADIUS] * image_number
            ,'cluster_size' : [CLUSTER_SIZE] * image_number
            ,'min_spot_per_cluster' : [MIN_SPOT_PER_CLUSTER] * image_number
        })

        if 'threshold' in Acquisition.columns : Acquisition = Acquisition.drop(columns='threshold')
        threshold_col_mask = Acquisition.columns.str.contains('Threshold') #looking for Threshold_0, Threshold_1 col
        if threshold_col_mask.any() :
            threshold_col = Acquisition.columns[threshold_col_mask]
            Detection_line_number = len(Detection)
            Acquisition.loc[:,threshold_col] = Acquisition.loc[:,threshold_col].fillna('').replace('',None)
            Detection = pd.merge(
                left= Detection,
                right= Acquisition.loc[:,['acquisition_id'] + list(Acquisition.columns[threshold_col_mask])],
                how='inner'
            )
            assert len(Detection) == Detection_line_number, "Duplicates or missing acquisition id in Acquisition df."

            Detection['threshold'] = None # above added threshold are like Threshold_i with i the color_number (0,1...)

        else :
            print("threshold column not found in Acquisition, automatic threshold will be used.")
            Detection['threshold'] = None

        id_columns = Detection.columns

        colors_columns = []
        for color_num, color in enumerate(colors) :
            colors_columns.append('{0}'.format(color_num))
            Detection['{0}'.format(color_num)] = color
        id_columns = list(id_columns)
        Detection = Detection.melt(
            id_vars= id_columns,
            value_vars= colors_columns,
            var_name= "color_id",
            value_name= "image",
        )
        Detection['visual_name'] = [None]*len(Detection)
        Detection = Detection.reset_index(drop=False, names='detection_id')
        Detection['detection_id'] += max_id +1
        max_id = Detection['detection_id'].max()

        for color_id in Detection['color_id'].unique() :
            target = 'Threshold_{0}'.format(color_id)
            if target in Detection.columns : 
                loc_index = Detection.loc[Detection['color_id'] == color_id].index
                Detection.loc[loc_index,['threshold']] = Detection[target]

        #Launching threads
        futures = []
        args_list = []
        keys = ['spots','spots_post_decomp','clustered_spots_dataframe','clusters_dataframe','clusters','clustered_spots','free_spots','threshold','voxel_size','spot_radius','alpha','beta','gamma','artifact_size','cluster_radius','min_spot_per_cluster',]
        with ProcessPool(max_workers=MAX_WORKERS) as executor:
            for args in zip(
                Detection['image'],
                Detection['voxel_size'],
                Detection['threshold'],
                Detection['spot_size'],
                Detection['alpha'],
                Detection['beta'],
                Detection['gamma'],
                Detection['artifact_radius'],
                Detection['cluster_size'],
                Detection['min_spot_per_cluster'],
                Detection['detection_id'],
            ):
                future = executor.schedule(multi_thread_full_detection, args = list(args), timeout=180)
                futures.append(future)
                args_list.append(args)

            detection_result = []
            for future, args in tqdm(zip(futures, args_list), total=len(futures)):
                try:
                    result = future.result()  # Set your timeout in seconds
                    
                except TimeoutError as e:
                    logging.warning(f"Detection timed out: {e}")
                    detection_id = args[-1]
                    result = dict.fromkeys(keys, np.nan)
                    result['detection_id'] = detection_id
                    logging.warning(f"Thread {detection_id} : Returning nan values")
                else :
                    detection_result.append(result)



        Spots, Clusters = build_Spots_and_Cluster_df(detection_result)

        #Correct coordinates for removed slices
        if type(DETECTION_SLICE_TO_REMOVE[0]) != type(None) :
            Spots['z'] = Spots['z'] + DETECTION_SLICE_TO_REMOVE[0]
            Clusters['z'] = Clusters['z'] + DETECTION_SLICE_TO_REMOVE[0]

        #Saving detection field view
        print("Saving field of views as compressed arrays.")
        save_dict = {
            str(detection_id) : np.max(image,axis=0) for detection_id,image in zip(Detection['detection_id'],Detection['image'])
        }

        np.savez(
            run_path + '/detection_fov/{0}'.format(location), #PATH,
            **save_dict
        )

        Detection['image_path'] = run_path + '/detection_fov/{0}.npz'.format(location)
        Detection['image_key'] = Detection['detection_id'].astype(str)
        Detection = Detection.drop(columns='image')

        #Appending dataframes
        Detection_save = pd.concat([
            Detection_save,
            Detection
            ],ignore_index=True, axis=0)
        
        Spots_save = pd.concat([
            Spots_save,
            Spots
            ], axis=0, ignore_index=True)

        Clusters_save = pd.concat([ 
            Clusters_save,  
            Clusters    
            ], axis=0, ignore_index=True)
        
        print(Detection_save['color_id'].unique())
        ###### End For loop #####

    #Unique Spots_identifier    
    Spots_save = Spots_save.drop(columns='spot_id').reset_index(drop=False, names="spot_id")

    #Explicit cast to int
    Detection_save["color_id"] = Detection_save["color_id"].astype(int)
    print("after cast to int : ",Detection_save['color_id'].unique())
    
    #Setting wavelength
    Detection_save['wavelength'] = None
    color_id_list = list(Detection_save['color_id'].unique())
    color_id_list.sort()

    if WAVELENGTH_LIST is None :
        WAVELENGTH_LIST = [None]*len(color_id_list)
    assert len(WAVELENGTH_LIST) == len(color_id_list)

    for color_id, wv in zip(color_id_list, WAVELENGTH_LIST) :
        Detection_save.loc[Detection_save['color_id'] == color_id, ['wavelength']] = wv
    Detection_save['wavelength'] = Detection_save['wavelength'].astype("Int16")
    print("lastone : ",Detection_save['color_id'].unique())

    #Saving results 
    Detection_save.to_feather(save_path + '/Detection.feather')
    Spots_save.to_feather(save_path + '/Spots.feather') 
    Clusters_save.to_feather(save_path + '/Clusters.feather')


def build_Spots_and_Cluster_df(detection_result : dict | list) :
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