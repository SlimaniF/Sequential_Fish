import os
import numpy as np

"""
This submodule checks that folder put as input contains file architecture compatible with pipeline.

in version 0.2 : nucleus signal was merged with fish signal so only fish folder remains to be checked.

"""

def _lvl1(RUN_PATH:str, fish_folder) :
    """
    returns True if ok else raise FileNotFoundError
    """

    dirlist = os.listdir(RUN_PATH)
    if not fish_folder in dirlist : raise FileNotFoundError("{0} folder not found in run folder.".format(fish_folder))

    return True

def _lvl2(RUN_PATH:str, fish_folder) :
    """
    returns locations list of ok else raise ValueError
    """
    Cy3_dirlist = os.listdir(RUN_PATH + "/{0}".format(fish_folder))

    locations_Cy3 = []

    for file in Cy3_dirlist :
        if 'Location-' in file : 
            locations_Cy3.append(file)

    locations_Cy3.sort()

    return locations_Cy3

def _lvl3(RUN_PATH, locations, fish_folder) :

    #Fish
    file_number = []
    file_dict = {}
    for location in locations :
        dirlist = os.listdir(RUN_PATH + '/{0}/'.format(fish_folder) + location)
        file_dict[location] = dirlist.copy()
        file_dict[location].sort()
        for file in dirlist : 
            if not file.endswith(".ome.tif") or file.startswith("._") : file_dict[location].remove(file)
        file_number.append(len(file_dict[location]))
    assert len(np.unique(file_number)) == 1, "Different file numbers found for fish Z-stacks amongst locations : {0}".format(np.unique(file_number))

    return file_dict

def assert_run_folder_integrity(run_path, fish_folder) :
    _lvl1(run_path, fish_folder=fish_folder)
    locations_list = _lvl2(run_path, fish_folder=fish_folder)
    file_dict = _lvl3(run_path, locations_list, fish_folder=fish_folder)

    return file_dict