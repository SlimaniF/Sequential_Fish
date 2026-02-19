import os
import numpy as np

def _lvl1(RUN_PATH:str, nucleus_folder, fish_folder) :
    """
    returns True if ok else raise FileNotFoundError
    """

    dirlist = os.listdir(RUN_PATH)
    if not fish_folder in dirlist : raise FileNotFoundError("{0} folder not found in run folder.".format(fish_folder))
    if not nucleus_folder in dirlist : raise FileNotFoundError("{0} folder not found in run folder.".format(nucleus_folder))

    return True

def _lvl2(RUN_PATH:str, nucleus_folder, fish_folder) :
    """
    returns locations list of ok else raise ValueError
    """
    Cy3_dirlist = os.listdir(RUN_PATH + "/{0}".format(fish_folder))
    Dapi_dirlist = os.listdir(RUN_PATH + "/{0}".format(nucleus_folder))

    locations_Cy3 = []
    locations_Dapi = []

    for file in Cy3_dirlist :
        if 'Location-' in file : 
            locations_Cy3.append(file)

    for file in Dapi_dirlist :
        if 'Location-' in file : 
            locations_Dapi.append(file)

    locations_Cy3.sort()
    locations_Dapi.sort()

    if locations_Cy3 != locations_Dapi :
        raise ValueError("Missmatch between locations found in Cy3 folder and in Dapi folder")
    else :
        return locations_Cy3

def _lvl3(RUN_PATH, locations, fish_folder, nucleus_folder) :

    #Dapi
    for location in locations :
        dirlist = os.listdir(RUN_PATH + '/{0}/'.format(nucleus_folder) + location)
        if len(dirlist) == 0 : raise FileNotFoundError("Dapi acquisition not found for location : {0}".format(location))
        elif len(dirlist) > 1 : raise FileNotFoundError("More than 1 dapi stack for location : {0}".format(location))
    
    #Cy3
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

def assert_run_folder_integrity(run_path, fish_folder, nucleus_folder) :
    _lvl1(run_path, fish_folder=fish_folder, nucleus_folder=nucleus_folder)
    locations_list = _lvl2(run_path, fish_folder=fish_folder, nucleus_folder=nucleus_folder)
    file_dict = _lvl3(run_path, locations_list, fish_folder=fish_folder, nucleus_folder=nucleus_folder)

    return file_dict