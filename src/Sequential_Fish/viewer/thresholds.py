"""
Widget to handle spot detection related tests in napari viewer.
"""

import os
import pandas as pd

from napari.layers import Image
import numpy as np
from tqdm import tqdm

from magicgui import magicgui
from magicgui.widgets import SpinBox
from napari.types import LayerDataTuple
from napari.components import ViewerModel
from .types import NapariWidget, UserInputError
from .load import LoadWidget

import bigfish.stack as stack
import bigfish.detection as detection
from bigfish.detection.utils import get_object_radius_pixel
from bigfish.detection import spots_thresholding, automated_threshold_setting

from typing import cast

class ThresholdsWidget(NapariWidget) :
    pass


_THRESHOLDS_WIDGETS : list = []
def register_thresholds_widget(cls) :
    _THRESHOLDS_WIDGETS.append(cls)
    return cls

def initiate_thresholds_widgets(
    **kwargs
) -> list:
    widgets_list = []
    for cls in _THRESHOLDS_WIDGETS :
        instance = cls(**kwargs)
        if hasattr(instance,"enabled") :
            if instance.enabled :
                widgets_list.extend(instance.widget)
        else :
            widgets_list.extend(instance.widget)

    return widgets_list

@register_thresholds_widget
class ThresholdSelector(LoadWidget) :
    """
    Widget aimed at helping user to set detection parameters : threshold, spot radius and so on...
    """

    def __init__(
        self,
        viewer: ViewerModel,
        default_threshold : int,
        default_spot_size : tuple,
        voxel_size : tuple,
        **_
        ) :
        
        self.viewer = viewer
        self.voxel_size = voxel_size
        self.dim = len(voxel_size)
        self.default_threshold = default_threshold
        self.spot_radius = default_spot_size
        self.kernel_size = None
        self.min_distance = None
        self.image = None
        self.filtered_image = None
        self.local_maxima = None
        self.layer_name = ""
        self.do_update = False
        
        super().__init__(viewer=viewer,table_dict=None)

    def _update_filtered_image(self) :

        if self.image is None : return None
        print("Computing filtered image.", end="", flush=True)
        self.filtered_image = _apply_log_filter(
            image=self.image,
            voxel_size=self.voxel_size,
            spot_radius=self.spot_radius,
            log_kernel_size= cast(tuple, self.kernel_size)
        )
        self.local_maxima = _local_maxima_mask(
            image_filtered=self.filtered_image,
            voxel_size=self.voxel_size,
            spot_radius=self.spot_radius,
            minimum_distance=cast(tuple, self.min_distance)
        )

    def _create_widget(self) :
        
        if not self.default_threshold is None and not self.filtered_image is None: 
            default_threshold = min(self.default_threshold, self.filtered_image.max())
        else :
            default_threshold = None

        @magicgui(
            threshold = {"widget_type" : SpinBox, "min" : 0, "value" : default_threshold},
            spot_radius = {"label" : "spot radius(zyx)", "value" : self.spot_radius},
        )
        def find_spots(
            threshold : int | None,
            spot_radius : tuple[int,int,int],
        ) :

            if (np.array(spot_radius) < 0).any() :
                raise ValueError("Spot radius : set value > 0 (0 to ignore argument)")
            
            if isinstance(spot_radius, tuple) :
                if not all(spot_radius) : spot_radius = cast(tuple[int,int,int], None) #any value set to 0

            if spot_radius != self.spot_radius :
                self.spot_radius = spot_radius
                self.do_update = True

            if not isinstance(self.viewer.layers.selection.active, Image) :
                self.do_update = False
                raise UserInputError("Selected layer must be an image.")
            
            elif "filtered image" in self.viewer.layers.selection.active.name :
                selected_name = self.viewer.layers.selection.active.name
                selected_name = selected_name.replace("filtered image","").strip()
                if selected_name in self.viewer.layers :
                    self.image = self.viewer.layers[selected_name]
                    if self.layer_name == selected_name :
                        pass
                    else :
                        self.layer_name = selected_name
                        self.do_update = True
                else :
                    self.do_update = False
                    raise UserInputError("Selected layer must not be a filtered image, select a signal image.")
                
            elif self.viewer.layers.selection.active.name != self.layer_name :
                self.image = self.viewer.layers.selection.active.data
                self.layer_name = self.viewer.layers.selection.active.name
                self.do_update = True
            
            if self.do_update :
                self._update_filtered_image()
                self.do_update = False
                self.widget.threshold.max = self.filtered_image.max() + 1 if not self.filtered_image is None else None
            
            self.filtered_image = cast(np.ndarray, self.filtered_image)
            self.local_maxima = cast(np.ndarray, self.local_maxima)
            ndim = self.filtered_image.ndim
            if threshold == 0 :
                print("Computing automated threshold : ...", end="", flush=True)
                if ndim == 4 :
                    shape = self.filtered_image.shape
                    threshold = cast(int, automated_threshold_setting(
                        self.filtered_image.reshape(shape[0]*shape[1],shape[2],shape[3]),
                        mask_local_max=self.local_maxima
                    ))

                else :
                    threshold = cast(int, automated_threshold_setting(
                        self.filtered_image,
                        mask_local_max=self.local_maxima
                    ))
                threshold = round(threshold)
                self.widget.threshold.value = threshold
                print("\rComputing automated threshold : done.")

            if ndim == 4 :
                location_index = 0
                spots = np.empty(shape=(0,4))
                for filtered_im, maxim_im in zip(self.filtered_image, self.local_maxima) :
                    new_spots = spots_thresholding(
                            filtered_im,
                            maxim_im,
                            threshold
                        )[0]
                    location_coords = np.array([location_index]*len(new_spots), dtype=int).reshape(-1,1)
                    new_spots = np.concatenate([location_coords, new_spots], axis=1)
                    location_index +=1
                    spots = np.concatenate([spots,new_spots],axis=0)

            else :

                spots = spots_thresholding(
                image=self.filtered_image,
                mask_local_max=self.local_maxima,
                threshold=threshold
                )[0]

            scale = self.voxel_size
    
            spot_layer_args = {
                'name' : f"{self.layer_name} detection",
                'size': 10, 
                'scale' : (1,) + scale if ndim == 4 else scale, 
                'face_color' : 'transparent', 
                'border_color' : 'red', 
                'symbol' : 'disc', 
                'opacity' : 0.7, 
                'blending' : 'translucent', 
                'units' : "nm",
                'visible' : True,
                }

            filtered_image_layer_args = {
                "colormap" :  'gray',
                "scale" : scale,
                "blending" : 'additive',
                "name" : f"{self.layer_name} filtered image",
                "projection_mode" : "max",
                'units' : "nm",
            }

            print(f"Thresholding done ({threshold})")

            return [
                    LayerDataTuple((self.filtered_image, filtered_image_layer_args, 'image')),
                    LayerDataTuple((spots, spot_layer_args, 'points'))
                    ]

        return find_spots

#Not registering to remove label
class ThreholdsFileEditor(ThresholdsWidget) :
    def __init__(self, cyclefile_path : str, color_number : int,  **_):

        if cyclefile_path.endswith("csv") :
            self.reader = pd.read_csv
            self.writer = pd.DataFrame.to_csv
            self.file_extension = ".csv"
        elif cyclefile_path.endswith("xlsx") :
            self.reader = pd.read_excel
            self.writer = pd.DataFrame.to_excel
            self.file_extension = ".xlsx"
        else :
            raise ValueError(f"unsupported extension for file : {cyclefile_path}.\nUse .csv or .x")

        if os.path.isfile(cyclefile_path) :
            self.enabled = True
            self.file = self.reader(cyclefile_path)
            for color in range(color_number) :
                target_col = f'Threshold_{color}'
                if not target_col in self.file.columns :
                    self.file[target_col] = pd.Series([pd.NA]*len(self.file), dtype=pd.Int64Dtype)
                else :
                    self.file[target_col] = self.file[target_col].astype(pd.Int64Dtype())

        else :
            print(f"{cyclefile_path} is not a valid file creating empty file.")
            self.file = pd.DataFrame(columns=pd.Index(["cycle", "gene", "Threshold_0"]))
            self.writer(self.file, cyclefile_path)
        
        self.path = cyclefile_path
        super().__init__()
        self.writer(self.file, self.path.replace(self.file_extension, f"_save{self.file_extension}"))

    def _create_widget(self) :

        @magicgui(threshold_table = {
            "widget_type" : "Table", 
            "value" : self.file,
            },
            call_button="Save")
        def edit_thresholds(threshold_table : pd.DataFrame)-> None :
            self.file = cast(dict[str,pd.Series],threshold_table)
            self.writer(pd.DataFrame(data= self.file["data"], columns= self.file["columns"]), self.path, index=False)
            print("Thresholds sucessfully saved.")
        
        return edit_thresholds




def _apply_log_filter(
        image: np.ndarray,
        voxel_size : tuple,
        spot_radius : tuple,
        log_kernel_size : tuple[float,float,float] | tuple[float,float],

) :
    """
    Apply spot detection steps until local maxima step (just before final threshold).
    Return filtered image.
    """
    
    ndim = image.ndim
    if type(log_kernel_size) == type(None) :
        log_kernel_size = get_object_radius_pixel(
                voxel_size_nm=voxel_size,
                object_radius_nm=spot_radius,
                ndim=3
                )
    

    if ndim == 4 :
        image_filtered = np.stack([stack.log_filter(im, log_kernel_size) for im in tqdm(image, desc= "Computing filtered image")])
    elif ndim == 3 :
        image_filtered = stack.log_filter(image, log_kernel_size)
    else :
        raise AssertionError(f"Unforseen dimension {ndim}")

    return image_filtered
    
def _local_maxima_mask(
    image_filtered: np.ndarray,
    voxel_size : tuple,
    spot_radius : tuple,
    minimum_distance : tuple,
    ) : 

    ndim = image_filtered.ndim

    if type(minimum_distance) == type(None) :
        minimum_distance = get_object_radius_pixel(
            voxel_size_nm=voxel_size,
            object_radius_nm=spot_radius,
            ndim=3)

    if ndim == 4 :
        mask_local_max = np.stack([detection.local_maximum_detection(im, minimum_distance) for im in tqdm(image_filtered, desc= "Computing local maxima")])
    elif ndim == 3 :
        mask_local_max = detection.local_maximum_detection(image_filtered, minimum_distance)
    else :
        raise AssertionError(f"Unforseen dimension {ndim}")
    
    return mask_local_max.astype(bool)