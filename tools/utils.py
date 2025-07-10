import pandas as pd
import numpy as np
import datetime as dt
import re, os, warnings
from skimage import io
from czifile import imread as _imread
from bigfish.stack import read_image as _read_image
from datetime import datetime

import warnings
from aicsimageio import AICSImage
from typing import Optional, Tuple


class MappingError(Exception) :
    """
    Raised when user inputs an incorrect image mapping.
    """
    pass

def auto_map_channels(image: np.ndarray, color_number: int, cycle_number: int, has_bead_channel = True) :
    """
    Assume z is the smallest spatial dimension
    """

    dim = image.ndim
    shape = image.shape
    reducing_list = list(shape)
    map_ = dict()

    try :
        c_idx = shape.index(color_number + has_bead_channel + 1 )# +1 for DAPI
    except ValueError :
        if has_bead_channel :
            error = "{0} colors channels (+2 for beads and dapi) are expected from experimental file but no matching axis was found in shape {1}.".format(color_number, shape)
        else :
            error = "{0} colors channels are expected from experimental file but no matching axis was found in shape {1}.".format(color_number, shape)
        raise MappingError(error)
    else :
        map_['c'] = c_idx
        reducing_list[c_idx] = max(reducing_list) + 1

    if dim > 4 : #Else dapi is being mapped and has only one cycle.
        try :
            cycles_idx = shape.index(cycle_number)
        except ValueError :
            raise MappingError("{0} cycles are expected from experimental file but no matching axis was found in shape {1}.".format(cycle_number, shape))
        else :
            map_['cycles'] = cycles_idx
            reducing_list[cycles_idx] = max(reducing_list) + 1

    #Set the smallest dimension to z
    z_val = min(reducing_list)
    z_idx = shape.index(z_val)
    map_['z'] = z_idx

    for index, value in enumerate(reducing_list) :
        if index in map_.values() : continue
        else : 
            if not 'y' in map_.keys() : map_['y'] = index
            else : map_['x'] = index

    return map_


def reorder_image_stack(image, channel_map, is_3D = True) :
    """
    will order image to cycles-zyxc
    """
    
    dim = image.ndim
    if is_3D :
        if dim == 5 :
            new_order = (channel_map['cycles'], channel_map['z'], channel_map['y'], channel_map['x'], channel_map['c'])
            ref_order = [0,1,2,3,4]
        elif dim == 4 :
            new_order = (channel_map['z'], channel_map['y'], channel_map['x'], channel_map['c'])
            ref_order = [0,1,2,3]
        else :
            raise AssertionError()
    else :
        if dim == 4 :
            new_order = (channel_map['cycles'], channel_map['y'], channel_map['x'], channel_map['c'])
            ref_order = [0,1,2,3]
        elif dim == 3 :
            new_order = (channel_map['y'], channel_map['x'], channel_map['c'])
            ref_order = [0,1,2]
        else :
            raise AssertionError()



    image = np.moveaxis(image, new_order, ref_order)
    return image


def open_image(path:str, image_number=None) :
    """
    Supports czi, png, jpg, jpeg, tif or tiff extensions.  
    If `image_number` is provided, extension must be .tif or .tiff
    """

    SUPPORTED_TYPES = ('.png', '.jpg', '.jpeg','.tif', '.tiff')

    if not image_number is None :
        if path.endswith('.tif') or path.endswith('.tiff') :
            im = [io.imread(path, plugin="tifffile", img_num=im_index) for im_index in range(image_number)]
            im = np.stack(im)
            return im
        else :
            warnings.warn("'image_number' is provided with a non tiff extension, ignoring 'image_number' argument")


    if path.endswith('.czi') :
        im = _imread(path)
    elif path.endswith(SUPPORTED_TYPES) :
        im = _read_image(path)
    else :
        raise ValueError("Unsupported type. Currently supported types are {0}".format(SUPPORTED_TYPES))

    return im

def get_datetime():
    return datetime.now().strftime("%d_%m_%Y_%H_%M_%S")


def _find_one_or_NaN(path, regex) :

    if isinstance(path, str) :
        res = re.findall(regex, path)
        if len(res) > 1 : raise ValueError("cycle regex yields ambiguous results.")
        elif len(res) == 0 : raise ValueError("cycle regex yields no results for {0}.".format(path))
        res = int(res[0])
    
    elif np.isnan(path) :
        res = np.NaN
    
    else :
        raise AssertionError("Unexpected element in Acquisition['full_path'] : {0}".format(path))
    
    return res


def safe_merge_no_duplicates(
        left : pd.DataFrame,
        right : pd.DataFrame,
        keys : 'list[str]',
        on : str = None,
        left_on : str = None,
        right_on : str = None,
        warn = False,
) :
    """
    Always perform left merge, aimed for 1:1 or m:1 merges. (error if duplicating or removing lines).
    Sole purpose is to safely add columns from right to left.
    """

    if type(keys) == str : keys  = [keys]

    if type(on) != type(None) :
        if type(on) == str : loc = [on]
        elif isinstance(on, (list,tuple)) : loc = on
        else : raise TypeError(f"Wrong type for 'on' parameter : {type(on)}")
    elif type(right_on) != type(None) : 
        if type(right_on) == str : loc = [right_on]
        elif isinstance(right_on, (list,tuple)) : loc = right_on
        else : raise TypeError(f"Wrong type for 'on' parameter : {type(right_on)}")
    else :
        raise ValueError("'on' parameter must be passed or couple ('left_on','right_on')")

    keys_to_merge = []
    for key in keys :
        if key not in left.columns : 
            keys_to_merge.append(key)
        else :
            if warn : warnings.warn(f"{key} already in left dataframe columns, {key} was removed from columns to merge.")

    if len(keys_to_merge) == 0 : 
        if warn : warnings.warn("No column to merge.")
        return left
    
    check_len = len(left)
    left = pd.merge(
        left,
        right.loc[:,loc + keys_to_merge],
        on=on,
        left_on=left_on,
        right_on=right_on,
    )

    if len(left) != check_len : raise ValueError(f"Lines were duplicated or removed during safe merge.\nPrevious count {check_len}; count after merge : {len(left)}")

    return left


def open_location(
        Acquisition : pd.DataFrame,
        location : str,
) :
    """
    Open all cycles of a location and reorder stacks in order (cycle,z,y,x,channel)
    """
    loc_Acquisition = Acquisition.loc[Acquisition['location'] == location].index
    assert len(loc_Acquisition) == 1, "Duplicates locations or no location found"

    fish_path = Acquisition.at[loc_Acquisition[0], 'full_path']
    fish_path_list = os.listdir(fish_path)
    fish_path_list.sort() # THIS MUST GIVE CYCLE ORDERED LIST ie : filename cycle matches map cycles and rest of filename doesn't change list order.
    fish_im = open_image(fish_path + fish_path_list[0]) #Opening first tiff file will open all tiff files of this location (multitif_file) with correct reshaping. Ignoring first dim which will be the cycles gives us image dimension

    stack_map = Acquisition.loc[Acquisition['location'] == location]['fish_map'].iat[0]   
    fish_im = reorder_image_stack(fish_im, channel_map=stack_map)

    return fish_im

def open_cycle(
        Acquisition : pd.DataFrame,
        location : str,
        cycle : int,
) :
    """
    Open specific cycle of a location and reorder stacks in order (z,y,x,channel)
    """
    loc_Acquisition = Acquisition.loc[Acquisition['location'] == location].index
    assert len(loc_Acquisition) == 1, "Duplicates locations or no location found"
    fish_path = Acquisition.at[loc_Acquisition[0], 'full_path']

    #Getting image informations
    fish_path_list = os.listdir(fish_path)
    fish_path_list.sort() # THIS MUST GIVE CYCLE ORDERED LIST ie : filename cycle matches map cycles and rest of filename doesn't change list order.
    stack_map = Acquisition.iat[(location,cycle), "fish_map"]
    stack_shape = Acquisition.iat[(location, cycle), "fish_shape"] 
    fullpath = Acquisition.iat[(location,cycle), "full_path"]
    
    #Preparing image shape
    z = stack_map['z']
    c = stack_map['c']
    image_number = stack_shape[z] * stack_shape[c]
    image_stack = open_image(fullpath, image_number= image_number)

    image_stack = reorder_image_stack(image_stack, channel_map=stack_map)

    return image_stack

def get_voxel_size_from_metadata(filepath: str) -> Optional[Tuple[Optional[int], Optional[int], Optional[int]]]:
    """
    Returns voxel size in nanometers (nm) as a tuple (X, Y, Z).
    Any of the dimensions may be None if not available.
    /WARINING\ : the unit might not be nm
    """
    try:
        img = AICSImage(filepath)
        voxel_sizes = img.physical_pixel_sizes  # values in meters
        if voxel_sizes is None:
            return None
        x = voxel_sizes.X * 1e3 if voxel_sizes.X else None
        y = voxel_sizes.Y * 1e3 if voxel_sizes.Y else None
        z = voxel_sizes.Z * 1e3 if voxel_sizes.Z else None
        return (z, y, x)
    except Exception as e:
        raise ValueError(f"Failed to read voxel size from {filepath}: {e}")
        return None
    
def get_min_cluster_radius(voxel_size) :
    return max(voxel_size)
    
    
def get_voxel_size(Detection : pd.DataFrame) :
    voxel_size = tuple(Detection['voxel_size'].iat[0])
    voxel_size = [int(i) for i in voxel_size]
    
    return voxel_size