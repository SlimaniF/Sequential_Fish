"""
This script aims at performing segmentation and savings results as .npy format to be used in FishSeq_pipeline_quantification.py
Drift correction is applied in FishSeq_pipeline_drift.py
"""
import os, warnings, sys
import numpy as np
import pandas as pd
import pbwrap.segmentation as segm
import bigfish.plot as plot

from tqdm import tqdm
from Sequential_Fish.tools import open_location

#### USER PARAMETERS

def main(run_path) :
    
    print(f"segmentation runing for {run_path}")
    
    if len(sys.argv) == 1:
        from default_pipeline_parameters import MODEL_DICT, OBJECT_SIZE_DICT, PLOT_VISUALS
    else :
        from Sequential_Fish.status import get_parameter_dict
        PARAMETERS = ['nucleus_model', 'cytoplasm_model', 'nucleus_size', 'cytoplasm_size', 'PLOT_VISUALS']
        parameters_dict = get_parameter_dict(run_path, parameters=PARAMETERS)        
        PLOT_VISUALS = parameters_dict['PLOT_VISUALS']
        nucleus_model = parameters_dict['nucleus_model']
        cytoplasm_model = parameters_dict['cytoplasm_model']
        nucleus_size = parameters_dict['nucleus_size']
        cytoplasm_size = parameters_dict['cytoplasm_size']

        MODEL_DICT = {
            'cytoplasm_model' : cytoplasm_model,
            'nucleus_model' : nucleus_model,
        }
        
        OBJECT_SIZE_DICT = {
            'nucleus_size' : nucleus_size,
            'cytoplasm_size' : cytoplasm_size
        }
        


    #Reading input folder.
    Acquisition = pd.read_feather(run_path + "/result_tables/Acquisition.feather")
    nuc_number = len(Acquisition['dapi_full_path'].unique())
    SAVE_PATH = run_path + "/segmentation/"
    VISUAL_PATH = run_path + "/visuals/segmentation/"
    os.makedirs(SAVE_PATH, exist_ok=True)
    os.makedirs(VISUAL_PATH, exist_ok=True)

    print("Starting segmentation pipeline, {0} nucleus images and {0} cytoplasm images to segment".format(nuc_number))

    for location in tqdm(Acquisition['location'].unique()) :
        sub_data = Acquisition.loc[Acquisition['location'] == location]

        #Setting output folder.
        image = open_location(Acquisition,location)
        nucleus_channel = sub_data['nucleus_channel'].iat[0]
        nucleus_image = image[..., nucleus_channel]

        #Nucleus_segmentation
        nucleus_image_save = nucleus_image.copy()
        nucleus_image = np.mean(nucleus_image, axis=0)
        nucleus_label = segm.Nucleus_segmentation(
            dapi=nucleus_image,
            diameter=OBJECT_SIZE_DICT['nucleus_size'],
            model_type= MODEL_DICT['nucleus_model'],
            use_gpu= True
        )

        #Cytoplasm segmentation
        cytoplasm_image = np.mean(image[...,:nucleus_channel], axis=4)
        cytoplasm_image = np.mean(cytoplasm_image, axis=0)

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
            dapi_signal = nucleus_image_save,
        )

        if PLOT_VISUALS : 
            plot.plot_segmentation_boundary(
                image=cytoplasm_image,
                cell_label=cytoplasm_label,
                nuc_label=nucleus_label,
                boundary_size=3,
                contrast=True,
                path_output=VISUAL_PATH + "/{0}_segmentation_cyto_view.png".format(location),
                show=False
            )
            plot.plot_segmentation_boundary(
                image=nucleus_image,
                cell_label=cytoplasm_label,
                nuc_label=nucleus_label,
                boundary_size=3,
                contrast=True,
                path_output=VISUAL_PATH + "/{0}_segmentation_nuc_view.png".format(location),
                show=False
            )
            
if __name__ == "__main__":
    if len(sys.argv) == 1:
        warnings.warn("Prefer launching this script with command : 'python -m Sequential_Fish pipeline input' or make sure there is no conflict for parameters loading in pipeline_parameters.py")
        from default_pipeline_parameters import RUN_PATH as run_path
    else :
        run_path = sys.argv[1]
    main(run_path)    