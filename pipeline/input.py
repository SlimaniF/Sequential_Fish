"""
This script aims at reading the input folder and preparing data folders and locations for next scripts.
"""
import sys
import Sequential_Fish.tools._folder_integrity as prepro
import pandas as pd
import os
import warnings
import numpy as np
from Sequential_Fish.tools.utils import auto_map_channels, _find_one_or_NaN, reorder_image_stack
from Sequential_Fish.tools import open_location
from Sequential_Fish.status import load_pipeline_parameters

def infer_channel(Cycle_map : pd.DataFrame, keyword= 'DAPI') :
    index = np.argmax(Cycle_map.eq(keyword).all(0)) - 1 # -1 since first columns is cycle number
    
    if index < 0 : raise ValueError(f'Could find {keyword} keyword in Cycle map.')
    
    return index

def main(run_path) :
    
    from Sequential_Fish import __version__

    pipeline_parameters = load_pipeline_parameters(run_path)
    MAP_FILENAME = pipeline_parameters.MAP_FILENAME
    cycle_regex = pipeline_parameters.cycle_regex
    CYCLE_KEY = pipeline_parameters.CYCLE_KEY
    GENES_NAMES_KEY = pipeline_parameters.GENES_NAMES_KEY
    WASHOUT_KEY_WORD = pipeline_parameters.WASHOUT_KEY_WORD
    FOLDER_KEYS = pipeline_parameters.FOLDER_KEYS 
    
    #Reading input folder.
    file_dict = prepro.assert_run_folder_integrity(
        run_path=run_path,
        fish_folder=FOLDER_KEYS.get('fish_folder'),
        )
    location_list = list(file_dict.keys())
    location_list.sort()
    location_number = len(location_list)
    print("{0} locations found.".format(location_number))

    #Init pandas DF
    COLUMNS = [
        "acquisition_id",
        "location",
        "cycle",
        "full_path",
        "fish_shape",
        "fish_map",
        "bead_channel", #TODO Add a way to infer bead channel or None if not found in cycle map
        "dapi_channel", #TODO Add a way to infer dapi channel or raise Error if not found
        "pipeline_version"
        ]
    Acquisition = pd.DataFrame(columns=COLUMNS)
    cycle_map = pd.read_excel(run_path + '/' + MAP_FILENAME)
    color_number = len(GENES_NAMES_KEY)
    cycle_number = len(cycle_map)
    print("Expected {0} colors.".format(color_number))
    print("Expected {0} cycles.".format(cycle_number))

    bead_channel = infer_channel(cycle_map, keyword= 'beads')
    dapi_channel = infer_channel(cycle_map, keyword= 'DAPI')
    has_bead = not bead_channel is None

    Acquisition['acquisition_id'] = np.arange(len(location_list)*cycle_number)
    Acquisition['location'] = location_list * cycle_number
    cycles_list = list(cycle_map[CYCLE_KEY])*location_number
    cycles_list.sort()
    Acquisition['cycle'] = cycles_list
    for location_index, location in enumerate(location_list) :
        index = Acquisition[Acquisition['location'] == location].index

        #Setting fish full path
        fish_path = run_path + "/{0}/{1}/".format(FOLDER_KEYS.get('fish_folder'), location)
        fish_path_list = os.listdir(fish_path)
        fish_path_list.sort() # THIS MUST GIVE CYCLE ORDERED LIST
        full_path_list = [fish_path + file for file in fish_path_list]
        Acquisition.loc[index, "full_path"] = full_path_list
        
        while len(full_path_list) < len(index) :
            full_path_list.append(np.NaN)

        fish_im = open_location(
            Acquisition=Acquisition,
            location=location
            )
        
        fish_map = auto_map_channels(fish_im, color_number=color_number, cycle_number=cycle_number, has_bead_channel= has_bead)
        fish_shape = fish_im.shape[:fish_map['cycles']] + fish_im.shape[(fish_map['cycles'] + 1):] #1cycle per acquisition
        reoderdered_shape = reorder_image_stack(fish_im, fish_map).shape
        fish_reodered_shape = reoderdered_shape[1:]

        Acquisition.loc[index, "fish_shape"] = pd.Series((fish_shape,)*cycle_number, index=index)
        Acquisition.loc[index, "fish_map"] = pd.Series((fish_map,)*cycle_number, index=index)
        Acquisition.loc[index, "fish_reodered_shape"] = pd.Series((fish_reodered_shape,)*cycle_number, index=index)

        cycle_regex_result = Acquisition.loc[:, 'full_path'].apply(_find_one_or_NaN, regex=cycle_regex)

    #Integrity checks
    assert all(Acquisition['cycle'].isin(cycle_map[CYCLE_KEY])), "Some cycle are not found in map"
    assert len(cycle_map[CYCLE_KEY] == len(Acquisition['cycle'])), "{0} column length doesn't match cycle number ({1})".format(CYCLE_KEY, len(Acquisition['cycle']))
    for key in GENES_NAMES_KEY : 
        assert len(cycle_map[key] == len(Acquisition['cycle'])), "{0} column length doesn't match cycle number ({1})".format(key, len(Acquisition['cycle']))

    cycle_regex_result = Acquisition.loc[:, 'full_path'].apply(_find_one_or_NaN, regex=cycle_regex)
    cycles_match = all(Acquisition.loc[~Acquisition['full_path'].isna(),"cycle"] == cycle_regex_result[~cycle_regex_result.isna()])
    if not cycles_match : raise ValueError("Missmatch between cycles assigned and cycles found in filenames. Maybe filenames could not be used to sort on cycles.")
    if any(Acquisition['full_path'].isna()) : warnings.warn("Warning : Some images registered in metadata were not found in folder. Ignore this message if some files were deleted after acquisition, in such a case pipeline should return as well 'OME series failed to read [...]. Missing data are zeroed' warning. ")

    Acquisition = pd.merge(
        left=Acquisition,
        right=cycle_map,
        left_on='cycle',
        right_on=CYCLE_KEY
    ).sort_values('acquisition_id').reset_index(drop=True)

    map_dict ={"cycle" : list(cycle_map[CYCLE_KEY])}
    map_dict.update({
        "{0}".format(gene_number) : list(cycle_map[gene_key]) for gene_number, gene_key in enumerate(GENES_NAMES_KEY)
    })

    color_columns = ["{0}".format(gene_number) for gene_number, gene_key in enumerate(GENES_NAMES_KEY)]
    Gene_map = pd.DataFrame(map_dict)
    Gene_map = Gene_map.melt(
        id_vars=['cycle'],
        value_vars=color_columns,
        value_name= "target",
        var_name="color_id"
    )
    Gene_map =Gene_map.reset_index(drop=False, names="map_id")
    washout_index = Gene_map[Gene_map['target'] == WASHOUT_KEY_WORD].index
    Gene_map.loc[washout_index, ['target']] = Gene_map.loc[washout_index]['target'] + '_' + Gene_map.loc[washout_index]['cycle'].astype(str) + '_' + Gene_map.loc[washout_index]['color_id'].astype(str)
    assert len(Gene_map['target']) == len(Gene_map['target'].unique()), "{1} duplicates found in Gene map even after washout renaming... If several cycle targets same RNA please add suffix in Gene map to differenciate.\nFound genes : \n{0}".format(Gene_map['target'], len(Gene_map['target']) - len(Gene_map['target'].unique()))

    #Set constant
    Acquisition['bead_channel'] = bead_channel
    Acquisition['dapi_channel'] = dapi_channel
    Acquisition['pipeline_version'] = __version__

    #Set index
    Gene_map = Gene_map.set_index('cycle', verify_integrity=True, drop=False)
    Acquisition = Acquisition.set_index(['location','cycle'], verify_integrity=True, drop=False)

    #Explicit dtype cast
    Gene_map['color_id'] = Gene_map['color_id'].astype(int)
    
    #Output
    save_path = run_path + '/result_tables/'
    os.makedirs(save_path, exist_ok=True)
    Acquisition.to_excel(save_path + '/Acquisition.xlsx')
    Acquisition.to_feather(save_path + '/Acquisition.feather')
    Gene_map.to_excel(save_path + 'Gene_map.xlsx')
    Gene_map.to_feather(save_path + 'Gene_map.feather')
    print("Done")
    
if __name__ == "__main__":
    warnings.warn("Prefer launching this script with command : 'python -m Sequential_Fish pipeline input'")
