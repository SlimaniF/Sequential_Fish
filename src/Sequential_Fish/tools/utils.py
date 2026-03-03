import re
import warnings
from typing import cast, Optional, Tuple
import datetime as dt

import pandas as pd
import numpy as np
import tifffile
from czifile import imread as _imread
from czifile import CziFile
from bigfish.stack import read_image as _read_image

from scipy.ndimage import distance_transform_edt

class MappingError(Exception) :
    """
    Raised when user inputs an incorrect image mapping.
    """

class NoVoxelSizeFound(Exception) :
    """
    Raised when could not infer voxel size from metadata.
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
            raise AssertionError(f"Uncorrect dimension : {dim} ; expected 5 or 4")
    else :
        if dim == 4 :
            new_order = (channel_map['cycles'], channel_map['y'], channel_map['x'], channel_map['c'])
            ref_order = [0,1,2,3]
        elif dim == 3 :
            new_order = (channel_map['y'], channel_map['x'], channel_map['c'])
            ref_order = [0,1,2]
        else :
            raise AssertionError(f"Uncorrect dimension : {dim} ; expected 4 or 3")



    image = np.moveaxis(image, new_order, ref_order)
    return image

def pad_to_shape(array : np.ndarray, new_shape) :
    shape = array.shape

    if len(array.shape) != len(new_shape) : raise ValueError("dimensions of array and new_shape don't match")

    pad_width_list = []

    for axis, axis_size in enumerate(shape) :
        target_size = new_shape[axis]
        
        pad_width = int(target_size - axis_size)
        if pad_width >= 0 :
            pad_width_list.append([0,pad_width])
        else :
            raise ValueError("Can't pad to new size {0} on axis {1} because current size {2} is bigger.".format(target_size, axis, axis_size))
    
    array = np.pad(array, pad_width_list)

    return array

def open_image(path:str, _map=None) :
    """
    Supports czi, png, jpg, jpeg, tif or tiff extensions.  
    If `image_number` is provided, extension must be .tif or .tiff
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
        keys : list[str] | str,
        on : str |list[str] | None = None,
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



def open_location(
        Acquisition : pd.DataFrame,
        location : str,
) :
    """
    Open all cycles of a location and reorder stacks in order (cycle,z,y,x,channel)
    """

    if not ('location' in Acquisition.index.names and 'cycle' in Acquisition.index.names) :
        Acquisition = Acquisition.set_index(['location','cycle'], verify_integrity=True)
    
    fish_path = Acquisition.at[(location,0), 'full_path']
    
    with tifffile.TiffFile(fish_path) as tif :
        location_stack = tif.asarray()

    stack_map = Acquisition.at[(location,0), 'fish_map']
    location_stack = reorder_image_stack(location_stack, channel_map=stack_map)

    return location_stack

def open_all_locations_one_cycle(
    Acquisition : pd.DataFrame,
    cycle: int,
) :

    if not ('location' in Acquisition.index.names and 'cycle' in Acquisition.index.names) :
        Acquisition = Acquisition.set_index(['location','cycle'], verify_integrity=True)

    location_list = list(Acquisition.index.get_level_values(0).unique())
    location_list.sort()

    max_shape_no_channel = np.max(Acquisition["fish_reodered_shape"].to_list(), axis=0)
    # max_shape_no_channel = max_shape_no_channel[:-1]
    location_stack = []

    for location in location_list :
        location = open_cycle(Acquisition,location,cycle)

        if not np.equal(location.shape, max_shape_no_channel).all() :
            location = pad_to_shape(location, new_shape=max_shape_no_channel)
        location_stack.append(location)

    if len(location_stack) > 1 :
        location_stack = np.stack(location_stack)
    else :
        location_stack = location_stack[0].reshape((1,) + location_stack[0].shape)
    return location_stack

def open_cycle(
        Acquisition : pd.DataFrame,
        location : str,
        cycle : int,
) :
    """
    Open specific cycle of a location and reorder stacks in order (z,y,x,channel)
    """
    if not ('location' in Acquisition.index.names and 'cycle' in Acquisition.index.names) :
        Acquisition = Acquisition.set_index(['location','cycle'], verify_integrity=True)

    #Getting image informations
    stack_map = Acquisition.at[(location,cycle), "fish_map"]
    stack_shape = Acquisition.at[(location, cycle), "fish_shape"] 
    fullpath = Acquisition.at[(location,cycle), "full_path"]
    stack_map = correct_map(stack_map)
    
    #Preparing image shape
    z = stack_map['z']
    c = stack_map['c']
    image_number = stack_shape[z] * stack_shape[c]

    print(stack_shape)
    print("image number : ", image_number)

    with tifffile.TiffFile(fullpath) as tif :
        cycle_stack = tif.asarray(key=range(0, image_number)).reshape(*stack_shape)

    image_stack = reorder_image_stack(cycle_stack, channel_map=stack_map)

    return image_stack

def get_voxel_size_from_metadata(filepath: str) -> Optional[Tuple[Optional[float], Optional[float], Optional[float]]]:
    """
    Returns voxel size in nanometers (nm) as a tuple (X, Y, Z).
    Any of the dimensions may be None if not available.
    /WARINING\\ : the unit might not be nm
    """
    try:
        if filepath.endswith('.czi'):
            with CziFile(filepath) as czi:
                metadata = cast(str,czi.metadata())  # returns XML metadata
                # try to parse voxel sizes from XML
                import xml.etree.ElementTree as ET
                root = ET.fromstring(metadata)
                scaling_distance = root.findall('.//Scaling//Items//Distance//Value')
                if len(scaling_distance) in [2,3] :
                    res = []
                    for scale in scaling_distance :
                        res = [float(cast(str,scale.text)) * 1e9 for scale in scaling_distance] #m to nm
                    res.reverse()
                    return tuple(res)
                else :
                    raise NoVoxelSizeFound("Couln't find voxel size on xml metadata")

        elif filepath.endswith(('.tif', '.tiff')):
            with tifffile.TiffFile(filepath) as tif:
                ij_meta = tif.imagej_metadata
                page = tif.pages[0]  # first image page
                # X/Y resolution as (numerator, denominator)
                xres = page.tags['XResolution'].value
                _ = page.tags['YResolution'].value
                # ResolutionUnit: must be 'nm' for this calculation
                res_unit = ij_meta.get("unit")

                if res_unit and str(res_unit) != 'nm':
                    xy_size = 1 / (xres[0] / xres[1]) * 1e3 #um to nm
                elif res_unit and str(res_unit) != 'NONE':
                    xy_size = 1 / (xres[0] / xres[1]) #um to nm
                else:
                    xy_size = None

                # Z spacing from ImageJ metadata
                if res_unit and str(res_unit) != 'nm':
                    z_size = ij_meta.get('spacing', None) * 1e3 
                else :
                    z_size = ij_meta.get('spacing', None)

                return (z_size,xy_size, xy_size )
    except NoVoxelSizeFound as e:
        print(f"Failed to read voxel size from {filepath}: {e}")
        return None
    
def get_min_cluster_radius(voxel_size) :
    return max(voxel_size)

def get_centroids_list(clusters_df) :

    """
    clusters_list should be a pd.DataFrame with ['z', 'y', 'x'] or ['y', 'x'] keys.
    """

    if 'y' in clusters_df.columns and 'x' in clusters_df.columns :
        if 'z' in clusters_df.columns : keys = [clusters_df['z'], clusters_df['y'], clusters_df['x']]
        else : keys = [clusters_df['y'], clusters_df['x']]
    else : raise ValueError("Expected keys : ['z', 'y', 'x'] or ['y', 'x']")

    return list(zip(*keys))
    
    
def get_voxel_size(Detection : pd.DataFrame) :
    voxel_size = tuple(Detection['voxel_size'].iat[0])
    voxel_size = [int(i) for i in voxel_size]
    
    return voxel_size

def get_centroids_array(cluster_df) :

    if len(cluster_df) == 0 :
        return np.empty(shape=(0,0), dtype=int)

    else : return np.array(get_centroids_list(cluster_df), dtype= int)

def _compute_critical_spot_number(radius_nm, voxel_size, density) :
    
    max_pixel_distance = int(max(nanometer_to_pixel(radius_nm, voxel_size)))
    kernel = np.ones(shape=(2*max_pixel_distance+1 ,2*max_pixel_distance+1, 2*max_pixel_distance+1)) #always odd number so middle is always at [pixel_radius-1, pixel_radius-1, pixel_radius-1]
    kernel[max_pixel_distance, max_pixel_distance, max_pixel_distance] = 0
    kernel = distance_transform_edt(kernel, sampling= voxel_size) <= radius_nm

    return int(round(kernel.sum() * density/100))


def nanometer_to_pixel(value, scale) :
    if isinstance(scale, (float,int)) : scale = [scale]
    if isinstance(value, (float,int)) : value = [value]*len(scale)
    if len(value) != len(scale) : raise ValueError("value and scale must have the same dimensionality")

    return list(np.array(value) / np.array(scale))

def compute_anisotropy_coef(voxel_size) :
    """
    voxel_size : tuple (z,y,x).
    """

    if not isinstance(voxel_size, (tuple, list)) : raise TypeError("Expected voxel_size tuple or list")
    if len(voxel_size) == 2 : is_3D = False
    elif len(voxel_size) == 3 : is_3D = True
    else : raise ValueError("Expected 2D or 3D voxel, {0} element(s) found".format(len(voxel_size)))

    if is_3D :
        z_anisotropy = voxel_size[0] / voxel_size [2]
        xy_anisotropy = voxel_size[1] / voxel_size [2]
        return (z_anisotropy, xy_anisotropy, 1)

    else :
        return (voxel_size[0] / voxel_size[1], 1)
    

def shift_array(arr : np.ndarray,*args) :
    indexer_new_array = []
    indexer_old_array = []
    for delta in args :
        if delta == 0 : 
            indexer_new_array.append(slice(None))
            indexer_old_array.append(slice(None))
        elif delta > 0 :
            indexer_new_array.append(slice(delta,None))
            indexer_old_array.append(slice(None,-delta))

        else :
            indexer_new_array.append(slice(None,delta))
            indexer_old_array.append(slice(-delta,None))

    if len(args) < arr.ndim :
        indexer_old_array.append(...)
        indexer_new_array.append(...)

    indexer_new_array = tuple(indexer_new_array)
    indexer_old_array = tuple(indexer_old_array)
    new_arr = np.zeros_like(arr)

    if len(args) > arr.ndim :
        raise ValueError("too many axis to shift; dim : {0}, shift : {1}".format(arr.ndim, args))
    else :
        new_arr[indexer_new_array] = arr[indexer_old_array]
        new_arr[indexer_new_array] = arr[indexer_old_array]

    return new_arr

def correct_map(_map:dict) : 
    """
    Maps needs to be corrected when used on image where cycles axis has been removed.
    """

    cycles_axis = _map['cycles']
    res = _map.copy()
    for key, value in res.items() :
        if value < cycles_axis :
            pass
        else :
            res[key] = value - 1
    
    res.pop('cycles')

    return res
