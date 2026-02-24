import numpy as np
import pandas as pd
from bigfish.multistack import match_nuc_cell
from smfishtools.preprocessing import shift_array
from concurrent.futures import ThreadPoolExecutor
from smfishtools.detection.multithread import cell_quantification
from tqdm import tqdm

from ..settings import get_settings
from ..tools import safe_merge_no_duplicates


def main(run_path) :

    print(f"quantification runing for {run_path}")
    
    pipeline_parameters = get_settings(run_path)
    MAX_WORKERS = pipeline_parameters.quantif_MAX_WORKERS
    
    Acquisition = pd.read_feather(run_path + '/result_tables/Acquisition.feather')
    Drift = pd.read_feather(run_path + '/result_tables/Drift.feather')
    Spots = pd.read_feather(run_path + '/result_tables/Spots.feather')
    Clusters = pd.read_feather(run_path + '/result_tables/Clusters.feather')
    Detection = pd.read_feather(run_path + '/result_tables/Detection.feather')

    #Preparing to add columns to Spots
    if not "in_nucleus" in Spots.columns : Spots['in_nucleus'] = False
    if not "in_nucleus" in Clusters.columns : Clusters['in_nucleus'] = False
    if not "cell_label" in Spots.columns : Spots['cell_label'] = np.nan
    if not "cell_label" in Clusters.columns : Clusters['cell_label'] = np.nan

    # Matching location with Drift
    print(Drift['acquisition_id'])
    Drift = safe_merge_no_duplicates(
        Drift,
        Acquisition,
        keys='location',
        on='acquisition_id'
    )

    #Matching location with Detection
    Detection = safe_merge_no_duplicates(
        Detection,
        Acquisition,
        on= 'acquisition_id',
        keys= 'location'
    )

    Cell_save = pd.DataFrame()
    Cell = pd.DataFrame()
    for location in Acquisition['location'].unique() :
        print("Starting {0}".format(location))

        #Opening segmentation
        segmentation_results = np.load(run_path + '/segmentation/{0}_segmentation.npz'.format(location))
        cytoplasm_label = segmentation_results['cytoplasm']
        nucleus_label = segmentation_results['nucleus']
        dapi_signal = segmentation_results['dapi_signal']

        #Correct drift for nucleus
        drift = Drift.loc[(Drift['drift_type'] == 'dapi') & (Drift['location'] == location), ['drift_y', 'drift_x']].to_numpy(dtype=int).squeeze()
        nucleus_label = shift_array(nucleus_label, *drift)

        nucleus_label, cytoplasm_label = match_nuc_cell(nucleus_label, cytoplasm_label, single_nuc=True, cell_alone=False)

        #Getting Detection ids for this fov
        sub_Detection = Detection.loc[Detection['location'] == location]
        selected_detection_id = sub_Detection['detection_id']

        # Adding cell label and if spots are in nuc seg
        sub_Spots = Spots.loc[Spots['detection_id'].isin(selected_detection_id)]
        sub_Spots_index = sub_Spots.index
        y = list(sub_Spots['y'])
        x = list(sub_Spots['x'])
        is_in_nucleus = nucleus_label[y,x].astype(bool) #0 is background
        spots_cell_labels = cytoplasm_label[y,x].astype(int)
        sub_Spots.loc[:,['in_nucleus']] = is_in_nucleus
        sub_Spots.loc[:,['cell_label']] = spots_cell_labels

        sub_Clusters = Clusters.loc[Clusters['detection_id'].isin(selected_detection_id)]
        sub_Clusters_index = sub_Clusters.index

        sub_spots_grouped = sub_Spots.groupby('cluster_id').agg({
            'in_nucleus' : 'max',
            'cell_label' : 'max'
        }).reset_index(drop=False)

        #Remove edge effect : uniformisation of labels within same cluster using max rule
        for index in sub_spots_grouped.index :
            cluster_id = sub_spots_grouped.at[index,'cluster_id']
            in_nucleus = sub_spots_grouped.at[index,'in_nucleus']
            cell_label = sub_spots_grouped.at[index,'cell_label']

            sub_Spots.loc[sub_Spots['cluster_id'] == cluster_id,['in_nucleus']] = in_nucleus
            sub_Spots.loc[sub_Spots['cluster_id'] == cluster_id,['cell_label']] = cell_label
            sub_Clusters.loc[sub_Clusters['cluster_id'] == cluster_id,['in_nucleus']] = in_nucleus
            sub_Clusters.loc[sub_Clusters['cluster_id'] == cluster_id,['cell_label']] = cell_label

        Spots.loc[sub_Spots_index] = sub_Spots
        Clusters.loc[sub_Clusters_index] = sub_Clusters

        sub_Spots = Spots.loc[
            (Spots['detection_id'].isin(selected_detection_id)) & (Spots['is_washout'] == False) & (Spots['cell_label'] > 0)
        ]

        sub_Clusters = Clusters.loc[
            (Clusters['detection_id'].isin(selected_detection_id)) & (Clusters['is_washout'] == False) & (Clusters['cell_label'] > 0)
        ]

        #Select all spots belonging to this fov; one list element per (cycle,color)
        all_fov_spots_lists = [
            np.array(list(
                zip(
                    sub_Spots[sub_Spots['detection_id'] == detection_id]['z'],
                    sub_Spots[sub_Spots['detection_id'] == detection_id]['y'],
                    sub_Spots[sub_Spots['detection_id'] == detection_id]['x'],
                )
            ), dtype=int)
        for detection_id in selected_detection_id]

        #Select all clusters belonging to this fov; one list element per (cycle,color)
        all_fov_clusters_lists = [
            np.array(list(
                zip(
                    sub_Clusters[sub_Clusters['detection_id'] == detection_id]['z'],
                    sub_Clusters[sub_Clusters['detection_id'] == detection_id]['y'],
                    sub_Clusters[sub_Clusters['detection_id'] == detection_id]['x'],
                    sub_Clusters[sub_Clusters['detection_id'] == detection_id]['spot_number'],
                    sub_Clusters[sub_Clusters['detection_id'] == detection_id]['cluster_id'].fillna(-1), #For bigfish compatibility
                )
            ))
        for detection_id in selected_detection_id]
        detection_number = len(selected_detection_id)

        #Launching threads on cell features
        detection_fov = np.load(run_path + '/detection_fov/{0}.npz'.format(location))
        fov_list = [detection_fov[fov_idx] for fov_idx in sub_Detection['image_key']]

        print("Starting individual cell metrics for {0} detections".format(len(sub_Detection)))
        with ThreadPoolExecutor(max_workers= MAX_WORKERS) as executor :
            cell_quantification_result = list(tqdm(executor.map(
                cell_quantification,
                sub_Detection['acquisition_id'],
                selected_detection_id,
                all_fov_spots_lists,
                all_fov_clusters_lists,
                sub_Detection['voxel_size'],
                [cytoplasm_label] * detection_number,
                [nucleus_label] * detection_number,
                fov_list,
                [dapi_signal]*detection_number
            ),total= len(sub_Detection), desc="individual cell metrics"))
        Cell = pd.concat(cell_quantification_result, axis=0) 
        Cell['location'] = location

        ## End of loop
        Cell_save = pd.concat([
            Cell_save,
            Cell
        ],axis=0)


    #Building cell_id (not unique identifier bc Cell actually has one line per detection)
    cell_IDs = Cell_save.groupby(['location','label'])['detection_id'].first().reset_index(drop=False).reset_index(drop=False, names='cell_id')
    cell_len = len(Cell_save)
    Cell_merged = pd.merge(
        left= Cell_save,
        right= cell_IDs.loc[:,['location','label','cell_id']],
        on= ('location','label'),
    )
    if len(Cell_merged) != cell_len : 
        print("\033[33mWARNING : Cell line conservation failed during merge. Saving a copy of the table before merge.\033[00m")
        Cell_save.reset_index(drop=True).to_feather(run_path + "/result_tables/Cell_before_merge.feather")
    else :
        del Cell_save
    Cell['detection_id'] = Cell['detection_id'].astype(int)


    #Save tables
    print("Saving results...")
    Cell_merged = Cell_merged.reset_index(drop=True).reset_index(drop=False, names='quantification_id')
    Cell_merged.to_feather(run_path + "/result_tables/Cell.feather")
    Spots.reset_index(drop=True).to_feather(run_path + "/result_tables/Spots.feather")
    Clusters.reset_index(drop=True).to_feather(run_path + "/result_tables/Clusters.feather")
    Drift.reset_index(drop=True).to_feather(run_path + "/result_tables/Drift.feather")
    Detection.reset_index(drop=True).to_feather(run_path + "/result_tables/Detection.feather")
    print("Quantification finished.")