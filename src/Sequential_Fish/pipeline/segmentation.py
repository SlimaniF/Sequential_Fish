"""
This script aims at performing segmentation and savings results as .npy format to be used in FishSeq_pipeline_quantification.py
Drift correction is applied in FishSeq_pipeline_drift.py
"""
import os
import numpy as np
import pandas as pd
import cellpose.models as models
import bigfish.plot as plot

from tqdm import tqdm
from ..tools.utils import open_image, reorder_image_stack
from ..settings import get_settings

#### USER PARAMETERS

def main(run_path) :
    
    print(f"segmentation runing for {run_path}")
    
    pipeline_parameters = get_settings(run_path)

    PLOT_VISUALS = pipeline_parameters.PLOT_VISUALS
    MODEL_DICT = pipeline_parameters.MODEL_DICT
    OBJECT_SIZE_DICT = pipeline_parameters.OBJECT_SIZE_DICT
    DO_3D_SEGMENTATION = pipeline_parameters.DO_3D_SEGMENTATION
    VOXEL_SIZE = pipeline_parameters.VOXEL_SIZE
    DAPI_CHANNEl = pipeline_parameters.DAPI_CHANNEl


    if DO_3D_SEGMENTATION :
        if len(VOXEL_SIZE) < 3 :
            raise ValueError(f"For 3D segmentation expected 3 dimensions in voxelsize : {VOXEL_SIZE}") 
        anisotropy = float(VOXEL_SIZE[0] / VOXEL_SIZE[1])
    else :
        anisotropy = 1.

    #Reading input folder.
    Acquisition = pd.read_feather(run_path + "/result_tables/Acquisition.feather")
    save_path = run_path + "/segmentation/"
    visual_path = run_path + "/visuals/segmentation/"
    os.makedirs(save_path, exist_ok=True)
    os.makedirs(visual_path, exist_ok=True)

    print("Starting segmentation pipeline")

    # Init cellpose models
    nucleus_model = models.CellposeModel(gpu= True, model_type = MODEL_DICT['nucleus_model'])
    cytoplasm_model = models.CellposeModel(gpu= True, model_type = MODEL_DICT['cytoplasm_model'])


    for location in tqdm(Acquisition['location'].unique()) :
        sub_data = Acquisition.loc[Acquisition['location'] == location]

        image_path = sub_data['full_path'].iat[0] #First washout, also avoid opening all images together.
        image_map = sub_data['fish_map'].iat[0] #First washout, also avoid opening all images together.
        image = open_image(image_path)
        image = reorder_image_stack(image, image_map)
        
        image = np.mean(image, axis=0)# mean on cycles
        if not DO_3D_SEGMENTATION :
            image = np.mean(image, axis=0)# mean on z

        #Nucleus_segmentation
        nucleus_image = image[...,DAPI_CHANNEl]
        nucleus_image_save = nucleus_image.copy()
        
        nucleus_label,*_ = nucleus_model.eval(
            nucleus_image, 
            diameter= OBJECT_SIZE_DICT['nucleus_size'], 
            do_3D= DO_3D_SEGMENTATION, 
            anisotropy=anisotropy
            )

        #Cytoplasm segmentation
        if image.shape[-1] > 2 :
            cytoplasm_image = np.mean(image[...,:DAPI_CHANNEl], axis=-1)
        else :
            cytoplasm_image = image[...,DAPI_CHANNEl-1]
        cytoplasm_label, *_ = cytoplasm_model.eval(
            cytoplasm_image, 
            diameter= OBJECT_SIZE_DICT['cytoplasm_size'], 
            do_3D= DO_3D_SEGMENTATION, 
            anisotropy=anisotropy
            )


        #Saving labels
        np.savez(
            file= save_path + "{0}_segmentation".format(location),
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
                path_output=visual_path + "/{0}_segmentation_cyto_view.png".format(location),
                show=False
            )
            plot.plot_segmentation_boundary(
                image=nucleus_image,
                cell_label=cytoplasm_label,
                nuc_label=nucleus_label,
                boundary_size=3,
                contrast=True,
                path_output=visual_path + "/{0}_segmentation_nuc_view.png".format(location),
                show=False
            )