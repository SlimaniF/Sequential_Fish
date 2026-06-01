"""
This script aims at performing segmentation and savings results as .npy format to be used in FishSeq_pipeline_quantification.py
Drift correction is applied in FishSeq_pipeline_drift.py
"""
import os
import logging
from typing import cast

import numpy as np
import pandas as pd
import cellpose.models as models
import bigfish.plot as plot
from tqdm import tqdm

from ..tools.utils import open_image, reorder_image_stack
from ..settings import get_settings
from ..settings import PipelineParameters

#### USER PARAMETERS

def main(run_path) :
    
    print(f"segmentation runing for {run_path}")
    
    pipeline_parameters = get_settings(run_path)
    pipeline_parameters = cast(PipelineParameters, pipeline_parameters)

    PLOT_VISUALS = pipeline_parameters.PLOT_VISUALS
    MODEL_DICT = pipeline_parameters.MODEL_DICT
    OBJECT_SIZE_DICT = pipeline_parameters.OBJECT_SIZE_DICT
    DO_3D_SEGMENTATION_NUCLEUS = pipeline_parameters.DO_3D_SEGMENTATION_NUCLEUS
    DO_3D_SEGMENTATION_CYTOPLASM = pipeline_parameters.DO_3D_SEGMENTATION_CYTOPLASM
    FLOW_3D_SMOOTH = pipeline_parameters.FLOW_3D_SMOOTH
    VOXEL_SIZE = pipeline_parameters.VOXEL_SIZE
    DAPI_CHANNEl = pipeline_parameters.DAPI_CHANNEl
    REFERENCE_CYCLE = pipeline_parameters.REFERENCE_CYCLE
    FLOW_THRESHOLD = pipeline_parameters.FLOW_THRESHOLD
    CELLPROB_THRESHOLD = pipeline_parameters.CELLPROB_THRESHOLD
    MIN_SIZE = pipeline_parameters.MIN_SIZE


    if DO_3D_SEGMENTATION_NUCLEUS or DO_3D_SEGMENTATION_CYTOPLASM :
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
    nucleus_model = models.CellposeModel(gpu= True, pretrained_model= MODEL_DICT['nucleus_model'])
    cytoplasm_model = models.CellposeModel(gpu= True, pretrained_model= MODEL_DICT['cytoplasm_model'])


    for location in tqdm(Acquisition['location'].unique()) :
        logging.info(f"Starting location : {location}")
        sub_data = Acquisition.loc[Acquisition['location'] == location]

        image_path = sub_data['full_path'].iat[0] #First washout, also avoid opening all images together.
        image_map = sub_data['fish_map'].iat[0] #First washout, also avoid opening all images together.
        image = open_image(image_path)
        image = reorder_image_stack(image, image_map)
        
        image = image[REFERENCE_CYCLE]

        #Nucleus_segmentation
        if not DO_3D_SEGMENTATION_NUCLEUS :
            nucleus_image = np.mean(image,axis=0)
        else :
            nucleus_image = image
        nucleus_image = nucleus_image[...,DAPI_CHANNEl]
        nucleus_image_save = nucleus_image.copy()
        
        nucleus_label,*_ = nucleus_model.eval(
            nucleus_image, 
            diameter= OBJECT_SIZE_DICT['nucleus_size'], 
            do_3D= DO_3D_SEGMENTATION_NUCLEUS,
            z_axis=0,
            flow3D_smooth=FLOW_3D_SMOOTH["nucleus"],
            anisotropy=anisotropy,
            cellprob_threshold=CELLPROB_THRESHOLD["nucleus"],
            flow_threshold=FLOW_THRESHOLD["nucleus"],
            min_size=MIN_SIZE["nucleus"]
            )

        #Cytoplasm segmentation
        if not DO_3D_SEGMENTATION_NUCLEUS :
            cytoplasm_image = np.mean(image,axis=0)
        else :
            cytoplasm_image = image

        if image.shape[-1] > 2 :
            cytoplasm_image = np.mean(cytoplasm_image[...,:DAPI_CHANNEl], axis=-1)
        else :
            cytoplasm_image = image[...,DAPI_CHANNEl-1]
        cytoplasm_label, *_ = cytoplasm_model.eval(
            cytoplasm_image, 
            diameter= OBJECT_SIZE_DICT['cytoplasm_size'], 
            do_3D= DO_3D_SEGMENTATION_CYTOPLASM, 
            anisotropy=anisotropy,
            flow_threshold=FLOW_THRESHOLD["cytoplasm"], #not used for 3D
            cellprob_threshold=CELLPROB_THRESHOLD["cytoplasm"],
            min_size=MIN_SIZE["cytoplasm"]
            )


        #Saving labels
        np.savez(
            file= save_path + "{0}_segmentation".format(location),
            nucleus= nucleus_label,
            cytoplasm= cytoplasm_label,
            dapi_signal = nucleus_image_save,
        )
        logging.info("Segmentation labels saved.")

        if PLOT_VISUALS :
            if DO_3D_SEGMENTATION_CYTOPLASM :
                logging.warning("Cannot represent 3D segmentation in 2D png file. To check cytoplasm segmentation results use viewer module.") 
            else :
                plot.plot_segmentation_boundary(
                    image=cytoplasm_image,
                    cell_label=cytoplasm_label,
                    nuc_label=nucleus_label,
                    boundary_size=3,
                    contrast=True,
                    path_output=visual_path + "/{0}_segmentation_cyto_view.png".format(location),
                    show=False
                )
            if DO_3D_SEGMENTATION_NUCLEUS :
                logging.warning("Cannot represent 3D segmentation in 2D png file. To check nucleus segmentation results use viewer module.") 
            else :
                plot.plot_segmentation_boundary(
                    image=nucleus_image,
                    cell_label=cytoplasm_label,
                    nuc_label=nucleus_label,
                    boundary_size=3,
                    contrast=True,
                    path_output=visual_path + "/{0}_segmentation_nuc_view.png".format(location),
                    show=False
                )
