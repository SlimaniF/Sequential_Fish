"""
This script aims at reading the input folder and preparing data folders and locations for next scripts.

Note : Initially this script is made to handle .ome.tif data. This format corresponds to tiff series with the particularity that metadata are stored exclusively in the cycle 0 file.
In case your system uses a different mechanism you will need to adapt this script to make a compatible output, by doing so you ensure downstream pipeline need not being change.

"""
import pandas as pd
import os
import warnings
from typing import TypedDict, cast
import tifffile
from ome_types import from_xml
from ome_types import OME
from tqdm import tqdm

import numpy as np

from ..settings import PipelineParameters
from ..tools.utils import auto_map_channels, _find_one_or_NaN, reorder_image_stack
from ..settings import get_settings

def main(run_path) :
    
    from Sequential_Fish import __version__

    pipeline_parameters = get_settings(run_path, settings_name="pipeline")
    pipeline_parameters = cast(PipelineParameters, pipeline_parameters)
    MAP_FILENAME = pipeline_parameters.MAP_FILENAME
    cycle_regex = pipeline_parameters.cycle_regex
    CYCLE_KEY = pipeline_parameters.CYCLE_KEY
    GENES_NAMES_KEY = pipeline_parameters.GENES_NAMES_KEY
    WASHOUT_KEY_WORD = pipeline_parameters.WASHOUT_KEY_WORD
    DAPI_CHANNEL = pipeline_parameters.DAPI_CHANNEl
    BEAD_CHANNEL = pipeline_parameters.BEAD_CHANNEl
    FISH_FOLDER = pipeline_parameters.FISH_FOLDER
    LOCATION_KEYWORD = pipeline_parameters.LOCATION_KEYWORD
    has_bead_channel = not BEAD_CHANNEL is None
    WAVELENGTH_LIST = pipeline_parameters.WAVELENGTH_LIST
    

    if WAVELENGTH_LIST is None :
        warnings.warn("WAVELENGTH were not set in parameters, pipeline will not be able to perform chromatic abberations corrections.")
    
    #Reading input folder.
    file_dict = assert_run_folder_integrity(
        run_path=run_path,
        fish_folder=FISH_FOLDER,
        location_key_word=LOCATION_KEYWORD
        )
    location_list = list(file_dict.keys())
    location_list.sort()
    location_number = len(location_list)
    print("{0} locations found.".format(location_number))

    #Init pandas DF
    COLUMNS = pd.Index([
        "acquisition_id",
        "location",
        "cycle",
        "full_path",
        "fish_shape",
        "fish_map",
        "bead_channel",
        "dapi_channel",
        "pipeline_version"
        ])
    Acquisition = pd.DataFrame(columns=COLUMNS)
    cycle_map = pd.read_excel(run_path + '/' + MAP_FILENAME)
    color_number = len(GENES_NAMES_KEY)
    cycle_number = len(cycle_map)
    print("Expected {0} colors.".format(color_number))
    print("Expected {0} cycles.".format(cycle_number))

    Acquisition['acquisition_id'] = np.arange(len(location_list)*cycle_number)
    Acquisition['location'] = location_list * cycle_number
    cycles_list = list(cycle_map[CYCLE_KEY])*location_number
    cycles_list.sort()
    Acquisition['cycle'] = cycles_list
    for location in tqdm(location_list) :

        #Setting dapi informations
        index = Acquisition[Acquisition['location'] == location].index

        #Setting general fish informations
        fish_path = run_path + "/{0}/{1}/".format(FISH_FOLDER, location)
        fish_path_list = os.listdir(fish_path)
        fish_path_list.sort() # THIS MUST GIVE CYCLE ORDERED LIST ie : filename cycle matches map cycles and rest of filename doesn't change list order.
        
        fish_im = None
        found_cycle_number = None
        if ".ome." in fish_path_list[0] :
            with tifffile.TiffFile(os.path.join(fish_path, fish_path_list[0])) as main_tif :
                #Axes map
                fish_map = map_from_metadata(main_tif.series[0].axes)
                if fish_map is None : #Could not infer from metadata
                    warnings.warn("Could not infer axes map from metadata opening image.")
                    fish_im = main_tif.asarray()
                    fish_map = cast(AxisMap, auto_map_channels(fish_im, color_number=color_number, cycle_number=cycle_number, has_bead_channel=has_bead_channel))
                
                #Shape
                fish_shape = None
                fish_reodered_shape = None
                #Infer from metadata
                if not main_tif.ome_metadata is None :
                    ome_metadata = from_xml(main_tif.ome_metadata)
                    fish_shape = get_memory_shape_from_metadata(ome_metadata, fish_map)
                    
                    if not fish_shape is None :
                        found_cycle_number = fish_shape[fish_map['cycles']]
                        fish_shape = fish_shape[:fish_map['cycles']] + fish_shape[(fish_map['cycles'] + 1):] #1cycle per acquisition

                    fish_reodered_shape = get_ordered_shape_from_metadata(ome_metadata)
                
                #Infer from open image : longer
                if main_tif.ome_metadata is None or fish_shape is None or fish_reodered_shape is None :
                    warnings.warn("Could not infer shape from metadata opening image.")
                    if fish_im is None : 
                        fish_im = main_tif.asarray()
                    fish_im = cast(np.ndarray, fish_im)
                    
                    fish_shape = fish_im.shape[:fish_map['cycles']] + fish_im.shape[(fish_map['cycles'] + 1):] #1cycle per acquisition
                    found_cycle_number = fish_im.shape[fish_map['cycles']]
                    fish_reodered_shape = reorder_image_stack(fish_im, fish_map).shape
                
                if found_cycle_number != cycle_number : raise UnevenCycleNumber(
                    cycle_number, 
                    found_cycle_number, 
                    f"Cycle file has {cycle_number} entries but only {found_cycle_number} were found in metadata."
                    )

                fish_reodered_shape = cast(tuple, fish_reodered_shape[1:])

        else :
            raise NotImplementedError("Initially this script is made to handle .ome.tif data. This format corresponds to tiff series with the particularity that metadata are stored exclusively in the cycle 0 file.\nIn case your system uses a different mechanism you will need to adapt this script to make a compatible output, by doing so you ensure downstream pipeline need not being change.")
        

        full_path_list = [fish_path + file for file in fish_path_list]
        while len(full_path_list) < len(index) :
            full_path_list.append(cast(str,np.nan))

        Acquisition.loc[index, "fish_shape"] = pd.Series((fish_shape,)*cycle_number, index=index)
        Acquisition.loc[index, "fish_map"] = pd.Series((fish_map,)*cycle_number, index=index)
        Acquisition.loc[index, "fish_reodered_shape"] = pd.Series((fish_reodered_shape,)*cycle_number, index=index)
        Acquisition.loc[index, "full_path"] = pd.Series(full_path_list, index=index, dtype="string")

    #Integrity checks
    assert all(Acquisition['cycle'].isin(cycle_map[CYCLE_KEY])), "Some cycle are not found in map"
    assert len(cycle_map) == len(Acquisition['cycle'].unique()), "{0} column length doesn't match cycle number ({1})".format(len(cycle_map), len(Acquisition['cycle']))

    cycle_regex_result = Acquisition.loc[:, 'full_path'].apply(_find_one_or_NaN, regex=cycle_regex)
    cycles_match = all(Acquisition.loc[~Acquisition['full_path'].isna(),"cycle"] == cycle_regex_result[~cycle_regex_result.isna()])
    if not cycles_match : raise ValueError("Missmatch between cycles assigned and cycles found in filenames. Maybe filenames could not be used to sort on cycles.")
    if any(Acquisition['full_path'].isna()) : 
        warnings.warn("Warning : Some images registered in metadata were not found in folder. Ignore this message if some files were deleted after acquisition, in such a case pipeline should return as well 'OME series failed to read [...]. Missing data are zeroed' warning. ")

    Acquisition = pd.merge(
        left=Acquisition,
        right=cycle_map,
        left_on='cycle',
        right_on=CYCLE_KEY
    ).sort_values('acquisition_id').reset_index(drop=True)
    Acquisition['dapi_channel'] = DAPI_CHANNEL
    Acquisition['bead_channel'] = BEAD_CHANNEL

    threshold_columns = Acquisition.columns[Acquisition.columns.str.contains("Threshold")]
    if Acquisition.loc[:, threshold_columns].isna().any().any() :
        raise ValueError("Found empty values in thresholds.")

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
    Acquisition['bead_channel'] = BEAD_CHANNEL
    Acquisition['dapi_channel'] = DAPI_CHANNEL
    Acquisition['pipeline_version'] = __version__

    #Explicit dtype cast
    Gene_map['color_id'] = Gene_map['color_id'].astype(int)
    
    #Set index
    Gene_map = Gene_map.reset_index(drop=True)
    Acquisition = Acquisition.reset_index(drop=True)
    
    #Output
    save_path = run_path + '/result_tables/'
    os.makedirs(save_path, exist_ok=True)
    Acquisition.to_excel(save_path + '/Acquisition.xlsx')
    Acquisition.to_feather(save_path + '/Acquisition.feather')
    Gene_map.to_excel(save_path + 'Gene_map.xlsx')
    Gene_map.to_feather(save_path + 'Gene_map.feather')
    print("Done")


class AxisMap(TypedDict):
    x : int
    y : int
    z : int
    c : int
    cycles : int

def map_from_metadata(axes_str : str | None) -> AxisMap | None :

    axis_map = {}
    if axes_str is None or axes_str == "" :
        return None
    elif "Z" in axes_str.upper() and "C" in axes_str.upper() :
        axis_map["z"] =  axes_str.upper().index("Z")
        axis_map["c"] =  axes_str.upper().index("C")
    elif "I" in axes_str.upper() : #I = Intensity usually Z and C are mixed -> Z*C
        axis_map["z"] = axes_str.upper().index("I")
        axis_map["c"] = axis_map["z"] + 1
    else :
        return None #could not infer

    for i in ["x","y"] :
        if i.upper() in axes_str.upper() :
            axis_map[i] = axes_str.upper().index(i.upper())
        else :
            return None
    
    if "T" in axes_str.upper() : #cycles
        axis_map["cycles"] = axes_str.upper().index("T")    

    return cast(AxisMap, axis_map)

def get_memory_shape_from_metadata(ome_metadata : OME, im_map : AxisMap) -> tuple | None :
    pixels_informations = ome_metadata.images[0].pixels

    shape = []
    items = cast(list, im_map.items())
    for axe, _ in sorted(items, key= lambda t : t[1]) :
        match axe :
            case "cycles" :
                shape.append(pixels_informations.size_t)
            case "c" :
                shape.append(pixels_informations.size_c)
            case "z" :
                shape.append(pixels_informations.size_z)
            case "y" :
                shape.append(pixels_informations.size_y)
            case "x" :
                shape.append(pixels_informations.size_x)
            case _ :
                raise AssertionError("Uncorrect map keys")
    
    if not all([isinstance(i,int) for i in shape]) : 
        warnings.warn("Error in getting shape from metadata, metadata containing not int axes")
        return None
    if not len(shape) == 5 : 
        warnings.warn("Expected 5 dimensions in shape")
        return None
    

    return tuple(shape)

def get_ordered_shape_from_metadata(ome_metadata : OME) -> tuple | None :
    pixels_informations = ome_metadata.images[0].pixels

    shape = (pixels_informations.size_t, pixels_informations.size_z,pixels_informations.size_y,pixels_informations.size_x,pixels_informations.size_c)
    
    if not all([isinstance(i,int) for i in shape]) : 
        warnings.warn("Error in getting shape from metadata, metadata containing not int axes")
        return None

    return tuple(shape)

class InputIntegrityError(ValueError) :
    """
    Raised if files and metadata are not consistent between themself or with cycle mapping.
    """
    pass


class UnevenCycleNumber(InputIntegrityError) :
    """
    Raised when metadata and expected file number are not consistent
    """

    def __init__(self, expected_cycle, found_cycle, *args: object) -> None:
        super().__init__(*args)
        self.expected_cycle = expected_cycle
        self.found_cycle = found_cycle

class UnevenFileNumber(InputIntegrityError) :
    """
    Raise when file number are not consistent amongst locations.
    """
    pass

def _lvl1(RUN_PATH:str, fish_folder) :
    """
    returns True if ok else raise FileNotFoundError
    """

    dirlist = os.listdir(RUN_PATH)
    if not fish_folder in dirlist : raise FileNotFoundError("{0} folder not found in run folder.".format(fish_folder))

    return True

def _lvl2(
    RUN_PATH:str, 
    location_key_word : str,
    fish_folder
    ) :
    """
    returns locations list of ok else raise ValueError
    """
    Cy3_dirlist = os.listdir(RUN_PATH + "/{0}".format(fish_folder))

    locations_Cy3 = []

    for file in Cy3_dirlist :
        if location_key_word in file : 
            locations_Cy3.append(file)
    locations_Cy3.sort()

    return locations_Cy3

def _lvl3(RUN_PATH, locations, fish_folder) :

    file_number = []
    file_dict = {}
    for location in locations :
        dirlist = os.listdir(RUN_PATH + '/{0}/'.format(fish_folder) + location)
        file_dict[location] = dirlist.copy()
        file_dict[location].sort()
        for file in dirlist : 
            if not file.endswith(".ome.tif") or file.startswith("._") : file_dict[location].remove(file)
        file_number.append(len(file_dict[location]))
    if len(np.unique(file_number)) != 1 : 
        raise UnevenFileNumber("Different file numbers found for fish Z-stacks amongst locations : {0}".format(np.unique(file_number)))

    return file_dict

def assert_run_folder_integrity(run_path, fish_folder, location_key_word) :
    _lvl1(run_path, fish_folder=fish_folder)
    locations_list = _lvl2(run_path, fish_folder=fish_folder, location_key_word=location_key_word)
    file_dict = _lvl3(run_path, locations_list, fish_folder=fish_folder)

    return file_dict