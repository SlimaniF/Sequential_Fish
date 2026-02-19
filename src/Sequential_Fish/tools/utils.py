import pandas as pd
import numpy as np
import datetime as dt
import re
from czifile import imread as _imread
from bigfish.stack import read_image as _read_image

import warnings


class MappingError(Exception) :
    """
    Raised when user inputs an incorrect image mapping.
    """
    pass

def auto_map_channels(image: np.ndarray, color_number: int, cycle_number: int, bead_channel = True) :
    """
    Assume z is the smallest dimension
    """

    dim = image.ndim
    shape = image.shape
    reducing_list = list(shape)
    map_ = dict()

    try :
        c_idx = shape.index(color_number + bead_channel)
    except ValueError :
        if bead_channel :
            error = "{0} colors channels (+1 for beads) are expected from experimental file but no matching axis was found in shape {1}.".format(color_number, shape)
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


def reorder_image_stack(image, map) :
    """
    will order image to cycles-zyxc
    """
    
    dim = image.ndim
    if dim == 5 :
        new_order = (map['cycles'], map['z'], map['y'], map['x'], map['c'])
        ref_order = [0,1,2,3,4]
    elif dim == 4 :
        new_order = (map['z'], map['y'], map['x'], map['c'])
        ref_order = [0,1,2,3]
    else :
        raise AssertionError()

    image = np.moveaxis(image, new_order, ref_order)
    return image

# def reorder_shape_tuple(shape, map) :


def open_image(path:str, map=None) :
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
    
    if type(map) != type(None) :
        im =reorder_image_stack(im, map)

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
