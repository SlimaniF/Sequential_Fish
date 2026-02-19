import numpy as np
import os
from typing import Literal
from skimage import io
from tqdm import tqdm

def open_image(fullpath, image_number=1) :
    arrays = [io.imread(fullpath, plugin="tifffile", img_num=im_index) for im_index in range(image_number)]



    arrays = np.stack(arrays)
    return arrays

def open_segmentation(segmentation_folder_fullpath: str, locations : 'list[str]', object : Literal['nucleus','cytoplasm'], z_repeat=None) :
    """
    Open with sorting on 'location'
    """

    masks = []

    #Opening masks
    for location in tqdm(locations, desc = "opening {0} masks".format(object)) :
        new_mask = np.load(segmentation_folder_fullpath + '/{0}_segmentation.npz'.format(location))[object]

        if type(z_repeat) != type(None) :
            new_mask = np.repeat(
                new_mask[np.newaxis],
                repeats= z_repeat,
                axis=0
            )
        masks.append(new_mask)

    masks = np.stack(masks)
    return masks

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

def correct_map(map:dict) : 
    """
    Maps needs to be corrected when used on image where cycles axis has been removed.
    """

    cycles_axis = map['cycles']
    for key, value in map.items() :
        if value < cycles_axis :
            pass
        else :
            map[key] = value - 1
    
    map.pop('cycles')

    return map

def reshape_stack(stack : np.ndarray, image_map : dict, im_shape : tuple) :
    """
    Reshape and reorder image stack safely.
    """
    
    map_ = image_map.copy()

    if map_['c'] == 1 : #Nothing to do
        
        stack = stack.reshape(*im_shape)
        stack = reorder_image_stack(stack, map_)
        return stack
    
    #.reshape method has to be called with z,c,y,x order independently of order that was yield in input pipeline.
    else :
        target_shape = (im_shape[map_['z']],
                        im_shape[map_['c']],
                        im_shape[map_['y']],
                        im_shape[map_['x']]
                        )

        stack = stack.reshape(*target_shape)

    # now we want to reoder to usual zyxc
    for key in ['z','y','x'] :
        if map_[key] < map_['c'] and  map_[key] >= 1:
            map_[key] +=1

    map_['c'] = 1
    stack = reorder_image_stack(stack, map_)

    return stack