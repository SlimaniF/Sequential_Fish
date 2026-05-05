"""
Aim : Compute correlation matrix (Pearson coefficient) between all signals
Note correcting drift is necessary to ensure good chromatic abberations correction, as the model should be trained after alignement. On the contraty if single color no need to
correct drift as correlation coef is invariant to translation.
"""

from ..settings import get_settings, AnalysisParameters, PipelineParameters
from ..chromatic_abberrations import apply_polynomial_transform_to_signal, load_calibration
from ..tools import open_location, shift_array

import os
import logging
from typing import cast
from datetime import datetime
from time import time

import pandas as pd
import numpy as np
from tqdm import tqdm
from pebble import ThreadPool
from multiprocessing import shared_memory, Pool

def main(run_path) :

    log_file = run_path + "/signal_correlation_log.log"
    script_start = datetime.now()
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
    )

    parameters = cast(AnalysisParameters, get_settings(run_path, "analysis"))
    pipeline_parameters = cast(PipelineParameters,get_settings(run_path, "pipeline"))
    voxel_size = pipeline_parameters.VOXEL_SIZE
    washout_keyword = pipeline_parameters.WASHOUT_KEY_WORD
    
    #1.Check data tables are found
    tables_names = ["Acquisition","Drift","Detection","Gene_map",]
    result_path = os.path.join(run_path,"result_tables")
    os.makedirs(result_path, exist_ok=True)
    list_tables = os.listdir(result_path)

    tables_found = [not table + ".feather" in list_tables for table in tables_names]
    if any(tables_found) :
        raise FileNotFoundError(f"Couldn't found datatables {np.array(tables_names)[tables_found]} at {result_path}")


    #2.Load data table
    Acquisition = pd.read_feather(os.path.join(result_path,"Acquisition.feather"))
    Drift = pd.read_feather(os.path.join(result_path,"Drift.feather"))
    Detection = pd.read_feather(os.path.join(result_path,"Detection.feather"))
    Gene_map = pd.read_feather(os.path.join(result_path,"Gene_map.feather"))


    #3.
    logging.info("\n\nSignal correlation analysis is starting\n")
    mean_correlation_matrix = np.mean([_process_location(
        Acquisition,
        Drift=Drift,
        Gene_map=Gene_map,
        Detection=Detection,
        washout_keyword=washout_keyword,
        location=location,
        parameters=parameters,
        voxel_size=voxel_size,
        dapi_channel= pipeline_parameters.DAPI_CHANNEl
        ) for location in tqdm(list(Acquisition["location"].unique()), desc= "Processing location")], axis=0)

    logging.info("\nComputing correlation matrix end. Saving results\n")


    Gene_map = Gene_map.sort_values(["cycle","color_id"])
    columns_name = list(Gene_map[~Gene_map["target"].str.contains(washout_keyword)]["target"])
    correlation_df = pd.DataFrame(mean_correlation_matrix, columns= pd.Index(columns_name), index=pd.Index(columns_name))
    correlation_df = correlation_df.to_feather(os.path.join(run_path,"result_tables","signal_correlations.feather"))
    logging.info("\nDone.\n")

    return True

def _process_location(
    Acquisition : pd.DataFrame,
    Drift : pd.DataFrame,
    Detection : pd.DataFrame,
    Gene_map : pd.DataFrame,
    washout_keyword : str | None,
    location : str,
    parameters : AnalysisParameters,
    voxel_size : tuple[int,int,int],
    dapi_channel : int
    ) :
    
    clock = time()
    location_stack = open_location(Acquisition, location) #shaped (cycle,z,y,x,color)

    #Remove dapi
    location_stack = np.delete(location_stack, dapi_channel, axis=4)

    #Remove washout
    if not washout_keyword is None :
        clock = time()
        washout_slices = list(Gene_map.loc[Gene_map["target"].str.contains(washout_keyword)]["cycle"].unique())
        location_stack = np.delete(location_stack,washout_slices,axis=0)
    else :
        washout_slices = []
    shape = location_stack.shape
    location_stack = location_stack.reshape(shape[0]*shape[4],*shape[1:4]) # one line per target (cycle*color) flatten arrays

    if not parameters.reference_wavelength is None :

        def _correct_slice(
        reference_wavelength : int,
        index : int,
        drift : tuple[int,int,int],
        color_id_number : int,
        reference_color_id : int,
        wavelength_series : pd.Series,
        voxel_size: tuple[int,int,int],
        ) :
            if index % color_id_number == reference_color_id :
                return None
            else : wavelength = wavelength_series.iat[index % color_id_number]

            #1. Correct Drift
            location_slice = shift_array(location_stack[index], *drift)

            #2. Correct Chromatic abberations
            calibration = load_calibration(
                reference_wavelength=reference_wavelength, 
                corrected_wavelength=wavelength
                )

            location_slice = apply_polynomial_transform_to_signal(
                image= location_slice,
                poly= calibration["polynomial_features"],
                model_x=calibration["x_fit"],
                model_y=calibration["y_fit"],
                model_z=calibration["z_fit"],
                voxel_size= np.array(voxel_size, dtype=int)
            )

            return location_slice

        reference_color_id = Detection.loc[Detection["wavelength"] == parameters.reference_wavelength].iloc[0].at["color_id"]
        wavelength_series = Detection.groupby(["color_id","wavelength"],as_index=False)["wavelength"].first()["wavelength"] #ordered by color_id
        wavelength_series = cast(pd.Series, wavelength_series)
        color_id_number = len(Detection["color_id"].unique())

        location_acquistion_ids = list(Acquisition.loc[(Acquisition["location"] == location) & (~Acquisition["cycle"].isin(washout_slices))]["acquisition_id"])
        drift_list = Drift.set_index("acquisition_id",verify_integrity=True).sort_values("cycle").loc[location_acquistion_ids,["drift_z","drift_y","drift_x"]]
        drift_list = np.repeat(drift_list.to_numpy(), color_id_number, axis=0)

        with ThreadPool(max_workers=2, ) as pool :
            task_number = len(location_stack)
            assert len(drift_list) == len(location_stack), f"drift_list : {len(drift_list)} vs location_stack : {len(location_stack)} \n{drift_list}"
            
            futures = pool.map(
                _correct_slice,
                [parameters.reference_wavelength]* task_number,
                list(range(len(location_stack))),
                drift_list,
                [color_id_number]* task_number,
                [reference_color_id]* task_number,
                [wavelength_series]* task_number,
                [voxel_size]* task_number,
            )
            for index, res in tqdm(enumerate(futures.result()), desc="Correcting chromatic abberations", total=task_number) :
                if not res is None :
                    location_stack[index] = res
    
    location_stack = location_stack.reshape(shape[0]*shape[4],np.prod(shape[1:4])) # one line per target (cycle*color) flatten arrays
    correlation_matrix = np.corrcoef(location_stack)
    
    return correlation_matrix

    