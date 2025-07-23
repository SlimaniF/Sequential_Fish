"""
Aims at finding drift value for each field of view and store it into a dataframe : Drift
"""

import sys
import pandas as pd
import warnings
import pbwrap.preprocessing.alignement as prepro
from tqdm import tqdm

from Sequential_Fish.status import load_pipeline_parameters
from Sequential_Fish.tools import open_location

def main(run_path) :

    print(f"drift runing for {run_path}")

    pipeline_parameters= load_pipeline_parameters(run_path)
    VOXEL_SIZE = pipeline_parameters.VOXEL_SIZE
    BEAD_SIZE = pipeline_parameters.BEAD_SIZE

    Drift_columns = [
        'acquisition_id',
        'drift_type',
        'drift_z',
        'drift_y',
        'drift_x',
        'voxel_size',
        'removed_slices',
        'max_projection',

    ]
    Drift_save = pd.DataFrame(columns=Drift_columns)
    Drift_save['max_projection'] = Drift_save['max_projection'].astype(bool)

    ### MAIN ###
    Acquisition = pd.read_feather(run_path + "/result_tables/Acquisition.feather")

    for location in Acquisition['location'].unique() : 

        print('Starting ',location)
        sub_acq = Acquisition.loc[Acquisition['location'] == location]
        
        print("opening images...")
        image = open_location(Acquisition, location)
        dapi_channel = sub_acq['dapi_channel'].iat[0]

        #Selecting images
        image = image[...,dapi_channel]
        reference_image = image[0]

        #Reference fov has no drift
        ref_acquisition_id = sub_acq[sub_acq['cycle'] == 0]['acquisition_id'].iat[0]
        Drift = pd.DataFrame({
            'acquisition_id' : [ref_acquisition_id],
            'drift_z' : [0],
            'drift_y' : [0],
            'drift_x' : [0],
        })

        #preparing cycle loop
        stack_index = 0

        print("Computing drift values for drifted fish stack...")
        for acquisition_id in tqdm(sub_acq[sub_acq['cycle'] != 0]['acquisition_id']) : #is ordered by cycle which is image stack ordered.
            stack_index +=1
            drifted_fish = image[stack_index]

            max_proj = False
            fish_result = prepro.fft_phase_correlation_drift(
                reference_image= reference_image,
                drifted_image= drifted_fish,
            )

            if (fish_result['drift_z'], fish_result['drift_y'], fish_result['drift_x'],) == (0,0,0) :
                max_proj = True
                fish_result = prepro.fft_phase_correlation_drift(
                    reference_image= reference_image,
                    drifted_image= drifted_fish,
                )

            fish_result = pd.DataFrame(fish_result)
            fish_result['acquisition_id'] = acquisition_id
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
    Drift_save = Drift_save.reset_index(drop=True).reset_index(drop=False, names= 'drift_id')
    Drift_save.to_feather(run_path + '/result_tables/Drift.feather')
    Drift_save.to_excel(run_path + '/result_tables/Drift.xlsx')

    print("Done")
    
if __name__ == "__main__":
    if len(sys.argv) == 1:
        warnings.warn("Prefer launching this script with command : 'python -m Sequential_Fish pipeline drift'")
        from default_pipeline_parameters import RUN_PATH as run_path
    else :
        run_path = sys.argv[1]
    main(run_path)    