"""
This script aims at performing segmentation and savings results as .npy format to be used in FishSeq_pipeline_quantification.py
Drift correction is applied in FishSeq_pipeline_drift.py
"""
import os, warnings, sys
import numpy as np
import pandas as pd
import smfishtools.segmentation as segm
import bigfish.plot as plot

from tqdm import tqdm
from Sequential_Fish.tools.utils import open_image, reorder_image_stack

#### USER PARAMETERS

def main(run_path) :
    
    print(f"segmentation runing for {run_path}")
    
    if len(sys.argv) == 1:
        from Sequential_Fish.pipeline_parameters import MODEL_DICT, OBJECT_SIZE_DICT, PLOT_VISUALS
    else :
        from Sequential_Fish.run_saves import get_parameter_dict
        PARAMETERS = ['nucleus_model', 'cytoplasm_model', 'nucleus_size', 'cytoplasm_size', 'PLOT_VISUALS']
        pipeline_parameters = get_parameter_dict(run_path, parameters=PARAMETERS)        
        PLOT_VISUALS = pipeline_parameters['PLOT_VISUALS']
        nucleus_model = pipeline_parameters['nucleus_model']
        cytoplasm_model = pipeline_parameters['cytoplasm_model']
        nucleus_size = pipeline_parameters['nucleus_size']
        cytoplasm_size = pipeline_parameters['cytoplasm_size']

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

        #Nucleus_segmentation
        nucleus_path = sub_data['dapi_full_path'].unique()
        nucleus_map = sub_data['dapi_map'].iat[0]
        assert len(nucleus_path) == 1, '{}'.format(nucleus_path)
        nucleus_path = nucleus_path[0]
        nucleus_image = open_image(nucleus_path)
        nucleus_image = reorder_image_stack(nucleus_image, nucleus_map)
        assert nucleus_image.ndim == 4, nucleus_image.shape
        nucleus_image = nucleus_image[:,:,:,0]
        nucleus_image_save = nucleus_image.copy()
        nucleus_image = np.mean(nucleus_image, axis=0)
        nucleus_label = segm.Nucleus_segmentation(
            dapi=nucleus_image,
            diameter=OBJECT_SIZE_DICT['nucleus_size'],
            model_type= MODEL_DICT['nucleus_model'],
            use_gpu= True
        )

        #Cytoplasm segmentation
        cytoplasm_path = sub_data['full_path'].iat[0] #First washout, also avoid opening all images together.
        cytoplasm_map = sub_data['fish_map'].iat[0] #First washout, also avoid opening all images together.
        cytoplasm_image = open_image(cytoplasm_path)
        cytoplasm_image = reorder_image_stack(cytoplasm_image, cytoplasm_map)
        cytoplasm_image = cytoplasm_image[0,:,:,:,-1]
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
        from Sequential_Fish.pipeline_parameters import RUN_PATH as run_path
    else :
        run_path = sys.argv[1]
    main(run_path)    