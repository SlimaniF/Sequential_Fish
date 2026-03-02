import numpy as np
from typing import Literal
from tqdm import tqdm
from itertools import zip_longest
from math import ceil
from typing import cast

def open_segmentation(
        segmentation_folder_fullpath: str, 
        locations : 'list[str]', 
        target : Literal['nucleus','cytoplasm'], 
        z_repeat=None
        ) :
    """
    Open with sorting on 'location'
    """

    masks = []

    #Opening masks
    for location in tqdm(locations, desc = "opening {0} masks".format(target)) :
        new_mask = np.load(segmentation_folder_fullpath + '/{0}_segmentation.npz'.format(location))[target]
        new_mask = cast(np.ndarray, new_mask)

        if type(z_repeat) != type(None) and new_mask.ndim == 2:
            new_mask = np.repeat(
                new_mask[np.newaxis],
                repeats= z_repeat,
                axis=0
            )
        masks.append(new_mask)

    masks = np.stack(masks)
    return masks


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

def correct_map(_map:dict) : 
    """
    Maps needs to be corrected when used on image where cycles axis has been removed.
    """

    cycles_axis = _map['cycles']
    for key, value in _map.items() :
        if value < cycles_axis :
            pass
        else :
            _map[key] = value - 1
    
    _map.pop('cycles')

    return _map

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

def _get_red_colors() :
    return ["#D0312D", "#990F02", "#60100B", "#7E2811", "#4E0707", "#BC544B", "#680C07"]

def _get_orange_colors():
    return ["#ED7014", "#FCAE1E", "#B56727", "#BE5504", "#D67229", "#E34A27", "#FF8C00"]

def _get_yellow_colors():
    return ["#D6B85A", "#DFC98A", "#C8A951", "#E7C27D", "#BDA55D", "#E4D00A", "#FFEF00"]

def _get_green_colors():
    return ["#3CB043", "#3A5311", "#728C69", "#AEF359", "#5DBB63", "#028A0F", "#234F1E", "#568203", "#4CBB17", "#487800"]

def _get_blue_colors():
    return ["#3944BC", "#63C5DA", "#0A1172", "#281E5D", "#1338BE", "#48AAAD", "#016064", "#2832C2", "#1F456E", "#4682B4"]

def _get_purple_colors():
    return ["#A32CC4", "#7A4988", "#601A35", "#A1045A", "#663046", "#311432", "#9867C5", "#880085"]

def _get_pink_colors():
    return ["#FC94AF", "#F25278", "#FE7D6A", "#FD5DA8", "#E3256B", "#FF6EC7"]

def _get_brown_colors():
    return ["#4B371C", "#231709", "#795C34", "#CC7722", "#65350F", "#652A0E", "#8A3324"]

def _get_black_colors():
    return ["#000000"]

def _get_grey_colors():
    return ["#808080", "#373737", "#594D5B", "#3E3D53", "#9897A9", "#63645E"]

def get_colors_list(
    size:int = 100, 
    remove_black= False, 
    remove_grey = False, 
    remove_brown= False
    ) -> list[str]:
    """
    Get a list of color from matplotlib.colors of length 'size'.
    100 different shade in the library
    """
    if not isinstance(size, int) : raise TypeError("size should be an int, it is a {0}".format(type(size)))
    if size < 1 : raise ValueError("Size should be >= 1")

    red = _get_red_colors()
    yellow = _get_yellow_colors()
    green = _get_green_colors()
    blue = _get_blue_colors()
    purple = _get_purple_colors()
    brown = _get_brown_colors()
    pink = _get_pink_colors()
    orange = _get_orange_colors()
    black = _get_black_colors()
    grey = _get_grey_colors()

    color_list = list(sum([*zip_longest(red, green, blue, black, orange, purple, grey, yellow, brown, pink)],()))
    length = len(color_list)
    while None in color_list : color_list.remove(None)
    if size > length :
        iteration = ceil(size / length)
        return (color_list * iteration)[:size]

    if remove_black : 
        for color in _get_black_colors() : color_list.remove(color)
    if remove_grey :
        for color in _get_grey_colors() : color_list.remove(color)
    if remove_brown :
        for color in _get_brown_colors() : color_list.remove(color)

    return (color_list)[:size]