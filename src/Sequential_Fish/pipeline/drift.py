"""
Aims at finding drift value for each field of view and store it into a dataframe : Drift
"""

import os, sys
import pandas as pd
import numpy as np
import smfishtools.preprocessing.alignement as prepro
from Sequential_Fish.tools import open_image, reorder_image_stack
from tqdm import tqdm

def main(run_path) :

    print(f"drift runing for {run_path}")

    if len(sys.argv) == 1:
        from Sequential_Fish.pipeline_parameters import DRIFT_SLICE_TO_REMOVE, VOXEL_SIZE, BEAD_SIZE, DO_HIGHPASS_FILTER
    else :
        from Sequential_Fish.run_saves import get_parameter_dict
        PARAMETERS = ['DRIFT_SLICE_TO_REMOVE', 'VOXEL_SIZE', 'BEAD_SIZE', 'DO_HIGHPASS_FILTER']
        pipeline_parameters= get_parameter_dict(run_path, PARAMETERS)
        DRIFT_SLICE_TO_REMOVE = pipeline_parameters['DRIFT_SLICE_TO_REMOVE']
        VOXEL_SIZE = pipeline_parameters['VOXEL_SIZE']
        BEAD_SIZE = pipeline_parameters['BEAD_SIZE']
        DO_HIGHPASS_FILTER = pipeline_parameters['DO_HIGHPASS_FILTER']
        

    SAVE_PATH = run_path + '/visuals/'
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
    Drift_save = pd.DataFrame(columns=Drift_columns)
    Drift_save['max_projection'] = Drift_save['max_projection'].astype(bool)
    Drift_save['highpass_filter'] = Drift_save['highpass_filter'].astype(bool)

    ### MAIN ###

    Acquisition = pd.read_feather(run_path + "/result_tables/Acquisition.feather")

    for location in Acquisition['location'].unique() : 

        print('Starting ',location)
        plot_path = SAVE_PATH + '/drift/{0}/'.format(location)
        os.makedirs(plot_path,exist_ok=True)
        #
        print("opening images...")
        sub_acq = Acquisition.loc[Acquisition["location"] == location].sort_values('cycle')
        dapi_path = sub_acq['dapi_full_path'].unique()
        dapi_map = sub_acq['dapi_map'].iat[0]
        assert len(dapi_path) == 1, 'multiple files for dapi found : {0}'.format(len(dapi_path))
        dapi_path = dapi_path[0]

        fish_path = sub_acq.loc[sub_acq['cycle']==0]['full_path']
        fish_map = sub_acq.loc[sub_acq['cycle']==0]['fish_map']
        assert len(fish_path) == 1, 'multiple files for fish found : {0}'.format(len(fish_path))
        fish_path = fish_path.iat[0]
        fish_map = fish_map.iat[0]

        dapi_image = open_image(dapi_path)
        dapi_image = reorder_image_stack(dapi_image, dapi_map)
        fish_image_stack = open_image(fish_path)
        fish_image_stack = reorder_image_stack(fish_image_stack, fish_map)

        max_z = max(
            dapi_image.shape[0],
            fish_image_stack.shape[1]
        )

        #Selecting images
        fish_image_stack = fish_image_stack[...,-1] # Selecting beads channel
        assert len(sub_acq) == len(fish_image_stack) == len(sub_acq['acquisition_id'])
        fish_image_stack = fish_image_stack[:,DRIFT_SLICE_TO_REMOVE[0]:-DRIFT_SLICE_TO_REMOVE[1],...] # removing first slice that is noisy
        fish_reference_image = fish_image_stack[0]
        fish_other_images = fish_image_stack[1:]
        dapi_image = dapi_image[DRIFT_SLICE_TO_REMOVE[0]:-DRIFT_SLICE_TO_REMOVE[1],:,:,-1]

        #Reference fov has no drift
        ref_acquisition_id = sub_acq[sub_acq['cycle'] == 0]['acquisition_id'].iat[0]
        Drift = pd.DataFrame({
            'acquisition_id' : [ref_acquisition_id],
            'drift_type' : ['fish'],
            'drift_z' : [0],
            'drift_y' : [0],
            'drift_x' : [0],
        })

        #preparing fish loop
        stack_index = 0

        #dapi drift
        if fish_reference_image.shape != dapi_image.shape :
            matching_shape = tuple(
                np.min([fish_reference_image.shape, dapi_image.shape], axis=0)
                )
        else : matching_shape = fish_reference_image.shape
        indexer = tuple([slice(axis) for axis in matching_shape])
        dapi_results = prepro.fft_phase_correlation_drift(
            fish_reference_image[indexer],
            dapi_image[indexer],
            bead_size=BEAD_SIZE,
            voxel_size=VOXEL_SIZE
        )
        dapi_results = pd.DataFrame(dapi_results)
        dapi_results['acquisition_id'] = ref_acquisition_id
        dapi_results['drift_type'] = 'dapi'

        Drift = pd.concat([
            Drift,
            dapi_results
        ], axis=0)

        print("Computing drift values for drifted fish stack...")
        for acquisition_id in tqdm(sub_acq[sub_acq['cycle'] != 0]['acquisition_id']) : #is ordered by cycle which is image stack ordered.
            stack_index +=1
            drifted_fish = fish_image_stack[stack_index]

            max_proj = False
            fish_result = prepro.fft_phase_correlation_drift(
                reference_image= fish_reference_image,
                drifted_image= drifted_fish,
                bead_size=BEAD_SIZE,
                voxel_size=VOXEL_SIZE,
            )

            if (fish_result['drift_z'], fish_result['drift_y'], fish_result['drift_x'],) == (0,0,0) :
                max_proj = True
                fish_result = prepro.fft_phase_correlation_drift(
                    reference_image= np.max(fish_reference_image,axis=0),
                    drifted_image= np.max(drifted_fish, axis=0),
                    bead_size=BEAD_SIZE,
                    voxel_size=VOXEL_SIZE,
                )

            fish_result = pd.DataFrame(fish_result)
            fish_result['acquisition_id'] = acquisition_id
            fish_result['drift_type'] = 'fish'
            fish_result['max_projection'] = max_proj

            Drift = pd.concat([
                Drift,
                fish_result,
            ], axis=0)

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