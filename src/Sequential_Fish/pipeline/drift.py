"""
Aims at finding drift value for each field of view and store it into a dataframe : Drift
"""

import os
import multiprocessing
from typing import cast

import pandas as pd
import numpy as np
import smfishtools.preprocessing.alignement as prepro
from tqdm import tqdm
from pebble import ProcessPool
from multiprocessing.shared_memory import SharedMemory
import tifffile

from ..settings import PipelineParameters
from ..tools import reorder_image_stack
from ..settings import get_settings

def main(run_path) :


    pipeline_parameters = get_settings(run_path)
    pipeline_parameters = cast(PipelineParameters, pipeline_parameters)
    DRIFT_SLICE_TO_REMOVE = pipeline_parameters.DRIFT_SLICE_TO_REMOVE
    VOXEL_SIZE = pipeline_parameters.VOXEL_SIZE
    BEAD_SIZE = pipeline_parameters.BEAD_SIZE
    DO_HIGHPASS_FILTER = pipeline_parameters.DO_HIGHPASS_FILTER
    DAPI_CHANNEL = pipeline_parameters.DAPI_CHANNEl
    BEAD_CHANNEL = pipeline_parameters.BEAD_CHANNEl
    REFERENCE_CYCLE = pipeline_parameters.REFERENCE_CYCLE
    MAX_WORKERS = min(4, multiprocessing.cpu_count())
        
    print(f"Drift runing for {run_path}\nBEAD CHANNEL : {BEAD_CHANNEL}\nDAPI CHANNE : {DAPI_CHANNEL}\nREFERENCE_CYCLE : {REFERENCE_CYCLE}")

    save_path = run_path + '/visuals/'
    Drift_columns = [
        'acquisition_id',
        'drift_type',
        'drift_z',
        'drift_y',
        'drift_x',
        'voxel_size',
        'bead_size',
        'removed_slices',
        'highpass_filter',
        'max_projection',

    ]
    Drift_save = pd.DataFrame(columns=pd.Index(Drift_columns))
    Drift_save['max_projection'] = Drift_save['max_projection'].astype(bool)
    Drift_save['highpass_filter'] = Drift_save['highpass_filter'].astype(bool)

    ### MAIN ###

    Acquisition = pd.read_feather(run_path + "/result_tables/Acquisition.feather")

    if not DRIFT_SLICE_TO_REMOVE[1] is None : DRIFT_SLICE_TO_REMOVE[1] *= -1

    for location in Acquisition['location'].unique() :
        print(f"Starting location : {location}") 
        Drift = process_location(
            location,
            Acquisition=Acquisition,
            save_path=save_path,
            DRIFT_SLICE_TO_REMOVE=DRIFT_SLICE_TO_REMOVE,
            REFERENCE_CYCLE=REFERENCE_CYCLE,
            DAPI_CHANNEL=DAPI_CHANNEL,
            VOXEL_SIZE=VOXEL_SIZE,
            MAX_WORKERS=MAX_WORKERS
        )

        Drift_save = pd.concat([
            Drift_save,
            Drift
        ], axis=0)


    print("All locations computed. Saving results...")
    Drift_save['voxel_size'] = [VOXEL_SIZE] * len(Drift_save)
    Drift_save['bead_size'] = [BEAD_SIZE] * len(Drift_save)
    Drift_save['highpass_filter'] = [DO_HIGHPASS_FILTER] * len(Drift_save)
    Drift_save = Drift_save.reset_index(drop=True).reset_index(drop=False, names= 'drift_id')
    Drift_save.to_feather(run_path + '/result_tables/Drift.feather')
    Drift_save.to_excel(run_path + '/result_tables/Drift.xlsx')

    print("Done")

class DimensionError(ValueError) :
    """Exception raised when image has not expected dimension"""

def process_location(
    location : str,
    Acquisition : pd.DataFrame,
    save_path : str,
    DRIFT_SLICE_TO_REMOVE : list,
    REFERENCE_CYCLE : int,
    DAPI_CHANNEL : int,
    VOXEL_SIZE : tuple[int,int,int],
    MAX_WORKERS : int,
    ) :
    plot_path = os.path.join(save_path,"drift",location)
    os.makedirs(plot_path,exist_ok=True)
    
    sub_acq = Acquisition.loc[Acquisition["location"] == location].sort_values('cycle')
    path = sub_acq['full_path'].iat[0]
    image_map = sub_acq['fish_map'].iat[0]

    with tifffile.TiffFile(path) as tif :
        image = tif.asarray() 
    image = reorder_image_stack(image, image_map)
    if len(image.shape) != 5 : raise DimensionError(f"Expected 5 dimensions for image stack got {image.ndim}")
    image = image[:,DRIFT_SLICE_TO_REMOVE[0]:DRIFT_SLICE_TO_REMOVE[1]]
    dapi_image_stack = image[...,DAPI_CHANNEL]
    
    #Preparing memory share
    dapi_reference_shape = dapi_image_stack.shape
    dapi_reference_dtype = dapi_image_stack.dtype
    dapi_image_stack = dapi_image_stack.ravel()
    del image

    shared_dapi_mem = SharedMemory(
        create=True,
        size=dapi_image_stack.nbytes
    )
    shared_dapi = np.ndarray(
    dapi_image_stack.shape,
    dtype=dapi_image_stack.dtype,
    buffer=shared_dapi_mem.buf
    )
    try :
        np.copyto(shared_dapi, dapi_image_stack)

        ref_acquisition_id = sub_acq[sub_acq['cycle'] == REFERENCE_CYCLE]['acquisition_id'].iat[0]
        Drift = pd.DataFrame({
            'acquisition_id' : [ref_acquisition_id],
            'drift_type' : ['dapi'],
            'drift_z' : [0],
            'drift_y' : [0],
            'drift_x' : [0],
            'error' : [np.nan],
            'phasediff' : [np.nan],
            'max_projection' : [None],
            'cycle' : [REFERENCE_CYCLE]
        })


        location_loc = sub_acq.loc[sub_acq['cycle'] != REFERENCE_CYCLE,['acquisition_id', 'cycle']].sort_values("cycle")
        cycle_number = len(location_loc)
        with ProcessPool(max_workers=MAX_WORKERS) as task_pool :
            futures = task_pool.map(
                    process_cycle,
                    [(shared_dapi_mem.name, dapi_reference_shape, dapi_reference_dtype)]*cycle_number,
                    [REFERENCE_CYCLE] * cycle_number,
                    location_loc["cycle"].to_list(),
                    [VOXEL_SIZE] * cycle_number,
                    location_loc["acquisition_id"].to_list()
            )

            results = [res for res in tqdm(futures.result(), desc="processing cycles", total=cycle_number)]

    finally :
        shared_array = None  # Dereference
        shared_dapi_mem.close()
        shared_dapi_mem.unlink()
    
    Drift = pd.concat([
        Drift,
        pd.DataFrame(data=results, columns= Drift.columns)
    ], axis=0, ignore_index=True
    )

    Drift = Drift.sort_values("cycle")

    return Drift


def process_cycle(
    dapi_image_stack_tuple : tuple,
    reference_index : int,
    cycle_index : int,
    VOXEL_SIZE : tuple[int,int,int],
    acquisition_id : int,
    ) -> dict :

    shared_dapi_name, dapi_reference_shape, dapi_reference_dtype = dapi_image_stack_tuple
    shared_mem = SharedMemory(name=shared_dapi_name)
    dapi_stack = np.ndarray(
        dapi_reference_shape,
        dtype=dapi_reference_dtype,
        buffer=shared_mem.buf
    )

    # Finding drift with dapi
    dapi_results = prepro.fft_phase_correlation_drift(
        dapi_stack[reference_index],
        dapi_stack[cycle_index],
        voxel_size=VOXEL_SIZE
    )

    if (dapi_results['drift_z'], dapi_results['drift_y'], dapi_results['drift_x']) == (0,0,0) : #No drift found in 3D, try in 2D
        max_proj = True
        dapi_results = prepro.fft_phase_correlation_drift(
            reference_image= np.max(dapi_stack[reference_index],axis=0),
            drifted_image= np.max(dapi_stack[cycle_index], axis=0),
            voxel_size=VOXEL_SIZE,
        )
        dapi_results["drift_z"] = 0
    else : 
        max_proj = False

    dapi_results['drift_type'] = 'dapi'
    dapi_results['max_projection'] = max_proj
    dapi_results['acquisition_id'] = acquisition_id
    dapi_results['cycle'] = cycle_index

    return dapi_results