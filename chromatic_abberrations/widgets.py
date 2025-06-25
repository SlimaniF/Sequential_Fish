"""
Widgets for chromatic abberrations correction calibration.
"""

import os, joblib
from napari.layers import Points, Image
from typing import Tuple
from magicgui import magicgui
from bigfish.detection import detect_spots
from sklearn.linear_model import LinearRegression

from Sequential_Fish.chromatic_abberrations import CALIBRATION_FOLDER
from ..tools import get_datetime
from ..customtypes import NapariWidget
from .calibration import match_beads
from .calibration import fit_polynomial_transform_3d
from .calibration import update_calibration_index
from .correction import apply_polynomial_transform_3d_to_signal

_calibration_widgets = []

def register_calibration_widget(cls) :
    _calibration_widgets.append(cls)
    return cls

def initiate_all_calibration_widgets() :
    return {type(cls()) : cls() for cls in _calibration_widgets}
        


@register_calibration_widget
class BeadsDetector(NapariWidget) :
    """
    This widget allow user to detect beads centers using LoG filter + threshold.
    """
    def __init__(
            self,
            voxel_size = (200,97,97),
            beads_radius = (300,150,150)
            ):
        
        self.default_voxel_size = voxel_size
        self.default_beads_radius = beads_radius
        super().__init__()

    def _create_widget(self):

        @magicgui(
            beads_image = {'label' : 'Image'},
            threshold = {'label' : 'Threshold', 'min': 1, 'max': 2**64-1},
            voxel_size = {'annotation' : Tuple[int,int,int], 'label' : 'Voxel size (zyx)'},
            beads_radius = {'annotation' : Tuple[int,int,int], 'label' : 'Beads radius (zyx)'},
            auto_call=False,
            call_button="Detect beads"
                
        )
        def detect_beads(
            beads_image : Image,
            threshold : int = 490,
            beads_radius = self.default_beads_radius,
            voxel_size = self.default_voxel_size,
        ) -> Points :
            
            print("Detecting beads...")
            coordinates = detect_spots(
                beads_image.data,
                threshold=threshold,
                voxel_size=voxel_size,
                spot_radius=beads_radius
            )

            detected_beads = Points(
                data=coordinates, 
                ndim=len(voxel_size),
                name = f"{beads_image.name}_beads",
                blending='additive',
                scale=voxel_size,
                face_color='transparent',
                symbol='disc',
                size=8
                )
            
            print(f"Found {len(detect_beads.data)} beads for {beads_image.name}.")

            return detected_beads
        
        
        return detect_beads

@register_calibration_widget
class ChromaticAberrationCorector(NapariWidget) :
    def __init__(self, degree = 2):
        super().__init__()

        self.model_x = LinearRegression()
        self.model_y = LinearRegression()
        self.model_z = LinearRegression()
        self.inv_model_x = LinearRegression()
        self.inv_model_y = LinearRegression()
        self.inv_model_z = LinearRegression()
        self.calibration_folder = CALIBRATION_FOLDER
        self.voxel_size = (1,1,1)
        self.degree = degree
        self.timestamp = get_datetime()

    def _create_widget(self):
        """
        Perform calibration for chromatic abberration correction and create a layer with corrected signal to evaluate quality of fit.
        """

        @magicgui(
                image_abberation={'label' : 'Image to correct :'},
                spatial_reference_shifted={'label' : 'Points with aberrations'},
                spatial_reference={'label' : 'Points reference'},
                degree={'label' : 'Degree'},
                auto_call=False,
                call_button= "Correct chromatic aberrations",
        )
        def create_corrected_layer(
            image_abberation : Image,
            spatial_reference_shifted : Points,
            spatial_reference : Points,
            degree : int = self.degree,
        ) ->  Image :
            
            voxel_size = spatial_reference.scale
            self.voxel_size = tuple([int(v) for v in voxel_size]) # save as reference if user save calibration
            
            #Convert pixel coordinates to nm to account for anisotropy
            coords1 = spatial_reference.data * voxel_size
            coords2 = spatial_reference_shifted.data * voxel_size

            beads, dist = match_beads(
                coords1= coords1,
                coords2= coords2,
                max_dist= int(max(voxel_size) * 4)
            )

            poly, self.model_x, self.model_y, self.model_z = fit_polynomial_transform_3d(
                                                beads,
                                                dist, 
                                                degree=degree
                                                )
            
            poly, self.inv_model_x, self.inv_model_y, self.inv_model_z = fit_polynomial_transform_3d(
                                                beads,
                                                dist, 
                                                degree=degree
                                                )
            

            image_corrected = apply_polynomial_transform_3d_to_signal(
                image_abberation.data,
                poly=poly,
                model_x=self.model_x,
                model_y=self.model_y,
                model_z=self.model_z,
                voxel_size=voxel_size
            )

            return Image(
                data= image_corrected,
                name= f"{image_abberation.name}_corrected",
                scale= image_abberation.scale,
                blending='additive',
                colormap=image_abberation.colormap,
                interpolation2d= image_abberation.interpolation2d
            )

        self.timestamp = get_datetime()

        return create_corrected_layer
    
    def _create_save_widget(self) :
        """
        This widget allow user to save previously performed calibration.
        """

        @magicgui(
                auto_call=False, 
                call_button= "Save calibration"
                )
        def save_fit_model(
            reference_wavelength : int,
            corrected_wavelength : int
        ) :
            
            if not os.path.isdir(self.calibration_folder) : os.makedirs(self.calibration_folder)
            filename = self.calibration_folder + f"/{reference_wavelength}_{corrected_wavelength}_{self.timestamp}" 

            joblib.dump({
                'x_fit' : self.model_x,
                'y_fit' : self.model_y,
                'z_fit' : self.model_z,
                'x_inv_fit' : self.inv_model_x,
                'y_inv_fit' : self.inv_model_y,
                'z_inv_fit' : self.inv_model_z,
                'voxel_size' : self.voxel_size,
                'degree' : self.degree,
                'reference_wavelength' : reference_wavelength,
                'corrected_wavelength' : corrected_wavelength,
                'timestamp' : self.timestamp,
            },
            filename
            )

            update_calibration_index(
                calib_folder=self.calibration_folder,
                reference_wavelength=reference_wavelength,
                corrected_wavelength=corrected_wavelength,
                filename= filename
            )
        
        return save_fit_model