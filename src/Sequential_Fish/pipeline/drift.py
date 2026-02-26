"""
Aims at finding drift value for each field of view and store it into a dataframe : Drift
"""

import os
from typing import cast

import pandas as pd
import numpy as np
import smfishtools.preprocessing.alignement as prepro
from tqdm import tqdm

from ..customtypes import PipelineParameters
from ..tools import open_image, reorder_image_stack
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
    has_beads = not BEAD_CHANNEL is None 
        
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

    for location in Acquisition['location'].unique() : 

        print('Starting ',location)
        plot_path = save_path + '/drift/{0}/'.format(location)
        os.makedirs(plot_path,exist_ok=True)
        #
        print("opening images...")
        sub_acq = Acquisition.loc[Acquisition["location"] == location].sort_values('cycle')
        path = sub_acq['full_path'].iat[0]
        image_map = sub_acq['fish_map'].iat[0]


        image = open_image(path)
        image = reorder_image_stack(image, image_map)
        assert len(image.shape) == 5
        if not DRIFT_SLICE_TO_REMOVE[1] is None : DRIFT_SLICE_TO_REMOVE[1] *= -1
        image = image[:,DRIFT_SLICE_TO_REMOVE[0]:DRIFT_SLICE_TO_REMOVE[1]]


        ref_acquisition_id = sub_acq[sub_acq['cycle'] == REFERENCE_CYCLE]['acquisition_id'].iat[0]
        stack_index = 0
        Drift = pd.DataFrame({
            'acquisition_id' : [ref_acquisition_id],
            'drift_type' : ['dapi'],
            'drift_z' : [0],
            'drift_y' : [0],
            'drift_x' : [0],
        })
        

        #Selecting images
        if has_beads :
            fish_image_stack = image[...,BEAD_CHANNEL] # Selecting beads channel
            assert len(sub_acq) == len(fish_image_stack) == len(sub_acq['acquisition_id'])
            fish_reference_image = fish_image_stack[REFERENCE_CYCLE]
            Drift = pd.concat([
                Drift,
                pd.DataFrame({
                    'acquisition_id' : [ref_acquisition_id],
                    'drift_type' : ['fish'],
                    'drift_z' : [0],
                    'drift_y' : [0],
                    'drift_x' : [0],
                })
            ])
        else : fish_reference_image = None

        dapi_image_stack = image[...,DAPI_CHANNEL]
        dapi_reference = dapi_image_stack[1]

        for acquisition_id in tqdm(sub_acq[sub_acq['cycle'] != REFERENCE_CYCLE]['acquisition_id']) : #is ordered by cycle which is image stack ordered.

            # Finding drift for dapi
            dapi_results = prepro.fft_phase_correlation_drift(
                dapi_reference,
                dapi_image_stack[stack_index],
                voxel_size=VOXEL_SIZE
            )
            dapi_results = pd.DataFrame(dapi_results)

            if (dapi_results.at[0,'drift_z'], dapi_results.at[0,'drift_y'], dapi_results.at[0,'drift_x'],) == (0,0,0) : #No drift found in 3D, try in 2D
                max_proj = True
                dapi_results = prepro.fft_phase_correlation_drift(
                    reference_image= np.max(dapi_reference,axis=0),
                    drifted_image= np.max(dapi_image_stack[stack_index], axis=0),
                    voxel_size=VOXEL_SIZE,
                )
                dapi_results = pd.DataFrame(dapi_results)
                dapi_results["drift_z"] = 0
            else : 
                max_proj = False

            dapi_results['drift_type'] = 'dapi'
            dapi_results['max_projection'] = max_proj
            dapi_results['acquisition_id'] = acquisition_id
            Drift = pd.concat([
                Drift,
                dapi_results
            ], axis=0)

            if has_beads :
                assert not fish_reference_image is None
                max_proj = False
                fish_result = prepro.fft_phase_correlation_drift(
                    reference_image= fish_reference_image,
                    drifted_image= fish_reference_image[stack_index],
                    bead_size=BEAD_SIZE,
                    voxel_size=VOXEL_SIZE,
                )

                if (fish_result['drift_z'], fish_result['drift_y'], fish_result['drift_x'],) == (0,0,0) : #No drift found in 3D, try in 2D
                    max_proj = True
                    fish_result = prepro.fft_phase_correlation_drift(
                        reference_image= np.max(fish_reference_image,axis=0),
                        drifted_image= np.max(fish_reference_image[stack_index], axis=0),
                        bead_size=BEAD_SIZE,
                        voxel_size=VOXEL_SIZE,
                    )

                fish_result = pd.DataFrame(fish_result)
                if max_proj : 
                    fish_result["drift_z"] = 0
                fish_result['acquisition_id'] = acquisition_id
                fish_result['drift_type'] = 'fish'
                fish_result['max_projection'] = max_proj

                Drift = pd.concat([
                    Drift,
                    fish_result,
                ], axis=0)
            stack_index +=1

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