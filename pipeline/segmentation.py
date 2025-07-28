"""
This script aims at performing segmentation and savings results as .npy format to be used in FishSeq_pipeline_quantification.py
Drift correction is applied in FishSeq_pipeline_drift.py
"""
import os, warnings, sys
import numpy as np
import pandas as pd
import pbwrap.segmentation as segm #TODO integrate in package
from tqdm import tqdm

from Sequential_Fish.tools import open_location, safe_merge_no_duplicates, shift_array
from Sequential_Fish.status import load_pipeline_parameters

#### USER PARAMETERS

def main(run_path) :
    
    print(f"segmentation runing for {run_path}")

    pipeline_parameters = load_pipeline_parameters(run_path)
    MODEL_DICT = pipeline_parameters.MODEL_DICT
    OBJECT_SIZE_DICT = pipeline_parameters.OBJECT_SIZE_DICT

    #Preparing drift correction
    Acquisition = pd.read_feather(run_path + '/result_tables/Acquisition.feather')
    Drift = pd.read_feather(run_path + '/result_tables/Drift.feather')
    # Matching location with Drift
    Drift = safe_merge_no_duplicates(
        Drift,
        Acquisition,
        keys=['location','cycle'],
        on='acquisition_id'
    ).sort_values(["location","cycle"])

    #Reading input folder.
    Acquisition = pd.read_feather(run_path + "/result_tables/Acquisition.feather")
    SAVE_PATH = run_path + "/segmentation/"
    os.makedirs(SAVE_PATH, exist_ok=True)

    locations = Acquisition['location'].unique()
    print(f"Starting segmentation pipeline, {len(locations)} field of view.")

    for location in tqdm(locations, total= len(locations)) :
        sub_data = Acquisition.loc[Acquisition['location'] == location]

        #Setting output folder.
        image = open_location(Acquisition,location)
        nucleus_channel = sub_data['dapi_channel'].iat[0]
        nucleus_image = image[..., nucleus_channel]

        #Correct drift
        drift_array = Drift.loc[Drift['location'] == location, ['drift_z','drift_y', 'drift_x']].to_numpy(dtype=int)
        for cycle, signal in enumerate(nucleus_image) :
            assert signal.ndim == 3, "Uncorrect indexing, signal dimension should be 3"
            nucleus_image[cycle] = shift_array(signal, *drift_array[cycle])

        cytoplasm_image = np.mean(image[...,:nucleus_channel], axis=4)
        print("before drift : ", cytoplasm_image.shape,"\n")
        for cycle, signal in enumerate(cytoplasm_image) :
            assert signal.ndim == 3, "Uncorrect indexing, signal dimension should be 3"
            cytoplasm_image[cycle] = shift_array(signal, *drift_array[cycle])
        
        # #Nucleus_segmentation
        nucleus_image = np.mean(nucleus_image, axis=0)
        nucleus_label = segm.Nucleus_segmentation(
            dapi=nucleus_image,
            diameter=OBJECT_SIZE_DICT['nucleus_size'],
            model_type= MODEL_DICT['nucleus_model'],
            use_gpu= True
        )


        #Cytoplasm segmentation
        print("cytoplasm_image before mean shape : ", cytoplasm_image.shape)
        cytoplasm_image = np.mean(cytoplasm_image, axis=0) #average on cycles
        print("cytoplasm_image shape : ", cytoplasm_image.shape)
        print("nucleus_image shape : ", nucleus_image.shape)
        # cytoplasm_image = np.mean(cytoplasm_image, axis=0) test segmentation 3D

        #Segmentation
        cytoplasm_label = segm.Cytoplasm_segmentation(
            cy3=cytoplasm_image,
            dapi=nucleus_image,
            diameter=OBJECT_SIZE_DICT['cytoplasm_size'],
            model_type=MODEL_DICT['cytoplasm_model'],
            use_gpu=True
        )

        #Saving labels
        np.savez(
            file= SAVE_PATH + "{0}_segmentation".format(location),
            nucleus= nucleus_label,
            cytoplasm= cytoplasm_label,
        )

            
if __name__ == "__main__":
    if len(sys.argv) == 1:
        warnings.warn("Prefer launching this script with command : 'python -m Sequential_Fish pipeline input' or make sure there is no conflict for parameters loading in pipeline_parameters.py") 