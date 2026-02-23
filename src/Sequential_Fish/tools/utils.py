import pandas as pd
import numpy as np
import datetime as dt
import re
from czifile import imread as _imread
from bigfish.stack import read_image as _read_image
from typing import cast

import warnings


class MappingError(Exception) :
    """
    Raised when user inputs an incorrect image mapping.
    """
    

def auto_map_channels(
    image: np.ndarray, 
    color_number: int, 
    cycle_number: int, 
    has_bead_channel : bool = False,
    ) :
    """
    Assume z is the smallest dimension
    """

    dim = image.ndim
    shape = image.shape
    reducing_list = list(shape)
    map_ = dict()

    try :
        c_idx = shape.index(color_number + has_bead_channel + 1) #+1 for dapi
    except ValueError as e:
        
        error = "{0} colors channels (+1 for dapi and optionally +1 for beads; ) are expected from experimental file but no matching axis was found in shape {1}.".format(color_number, shape)
        if str(e) == "tuple.index(x): x not in tuple" :
            raise MappingError(error) from e
        else : raise e
    else :
        map_['c'] = c_idx
        reducing_list[c_idx] = max(reducing_list) + 1

    if dim > 4 : #Else dapi is being mapped and has only one cycle.
        try :
            cycles_idx = shape.index(cycle_number)
        except ValueError as e:
            if str(e) == "tuple.index(x): x not in tuple" :
                raise MappingError("{0} cycles are expected from experimental file but no matching axis was found in shape {1}.".format(cycle_number, shape)) from e
            else : raise e
        else :
            map_['cycles'] = cycles_idx
            reducing_list[cycles_idx] = max(reducing_list) + 1

    #Set the smallest dimension to z
    z_val = min(reducing_list)
    z_idx = shape.index(z_val)
    map_['z'] = z_idx

    for index, _ in enumerate(reducing_list) :
        if index in map_.values() : continue
        else : 
            if not 'y' in map_.keys() : map_['y'] = index
            else : map_['x'] = index

    return map_


def reorder_image_stack(image, _map) :
    """
    will order image to cycles-zyxc
    """
    
    dim = image.ndim
    if dim == 5 :
        new_order = (_map['cycles'], _map['z'], _map['y'], _map['x'], _map['c'])
        ref_order = [0,1,2,3,4]
    elif dim == 4 :
        new_order = (_map['z'], _map['y'], _map['x'], _map['c'])
        ref_order = [0,1,2,3]
    else :
        raise AssertionError()

    image = np.moveaxis(image, new_order, ref_order)
    return image

# def reorder_shape_tuple(shape, map) :


def open_image(path:str, _map=None) :
    """
    Supports czi, png, jpg, jpeg, tif or tiff extensions.
    """

    SUPPORTED_TYPES = ('.png', '.jpg', '.jpeg','.tif', '.tiff')

    if path.endswith('.czi') :
        im = _imread(path)
    elif path.endswith(SUPPORTED_TYPES) :
        im = _read_image(path)
    else :
        raise ValueError("Unsupported type. Currently supported types are {0}".format(SUPPORTED_TYPES))
    
    if type(_map) != type(None) :
        im =reorder_image_stack(im, _map)

    return im

def get_datetime():
    return dt.datetime.now().strftime("%Y%m%d %H-%M-%S")


def _find_one_or_NaN(path, regex) :

    if isinstance(path, str) :
        res = re.findall(regex, path)
        if len(res) > 1 : raise ValueError("cycle regex yields ambiguous results.")
        elif len(res) == 0 : raise ValueError("cycle regex yields no results for {0}.".format(path))
        res = int(res[0])
    
    elif np.isnan(path) :
        res = np.nan
    
    else :
        raise AssertionError("Unexpected element in Acquisition['full_path'] : {0}".format(path))
    
    return res


def safe_merge_no_duplicates(
        left : pd.DataFrame,
        right : pd.DataFrame,
        keys : list[str],
        on : str | None = None,
        left_on : str | None = None,
        right_on : str | None = None,
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

    keys_to_merge : list[str] = []
    loc = cast(list[str],loc)
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

def inv_FWHM(FWHM) :
    """
    Returns gaussian variance for gaussian defined with Full Width at Half Maximum.
    """
    return 1/(2*np.sqrt(2* np.log(2))) * FWHM

def inv_FWTM(FWTM) :
    """
    Returns gaussian variance for gaussian defined with Full Width at Tenth Maximum.
    """
    return 1/(2*np.sqrt(2* np.log(2))) * FWTM

def gaussian_kernel_size(object_size_nm, voxel_size, width= 'FWHM') :
    """
    Computes the kernel size of a Gaussian function so that either Full Width at Half Maximum (width = 'FWHM') or Full Width at Tenth Maximum (width = 'FWTM') corresponds to the object dimension.
    Object_size_nm should be coherent in type (and length if tuples or lists are passed).
    Always returns a list with object dim number of element.
    """
    
    if width == 'FWHM' : variance_func = inv_FWHM
    elif width == 'FWTM' : variance_func = inv_FWTM
    else : raise ValueError("with should either be 'FWHM' or 'FWTM'. It is {0}".format(width))

    if isinstance(object_size_nm, (tuple, list)) and isinstance(voxel_size, (tuple, list)):
        if len(object_size_nm) != len(voxel_size) : raise ValueError("Length of object_size_nm and voxel_size parameters should be coherent.")
        return [variance_func(obj_size_nm / scale_factor) for obj_size_nm, scale_factor in zip(object_size_nm, voxel_size)]
    elif isinstance(object_size_nm, (float,int)) and isinstance(voxel_size, (tuple, list)): 
        return [variance_func(object_size_nm / scale_factor for scale_factor in voxel_size)]
    else : raise TypeError("object size and voxel_size parameters should be tuple or float like object and be coherent.")


