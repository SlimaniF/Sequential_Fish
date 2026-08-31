from typing import cast


from .types import ThreadedWidget, NapariWidget, UserInputError
from ..tools import get_datetime
from ..chromatic_abberrations import calibration_exist, load_calibration
from ..chromatic_abberrations.calibration import match_beads, fit_polynomial_transform_3d, save_fit_model
from ..chromatic_abberrations import apply_polynomial_transform_spots, apply_polynomial_transform_to_signal
from ..chromatic_abberrations import CALIBRATION_FOLDER


from napari import Viewer
from napari.types import LayerDataTuple
from napari.layers import Points, Image
from magicgui import magicgui

import numpy as np
from tqdm import tqdm
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures


class ChromaticWidget(ThreadedWidget) :

    def __init__(self, *,voxel_size : tuple[int,int,int], wavelength_list : list[int], viewer: Viewer):
        self.voxel_size =voxel_size
        self.wavelength_list = wavelength_list
        super().__init__(viewer=viewer)



_CHROMATIC_WIDGETS : 'list[NapariWidget]' = []
def register_chromatic_widget(cls) :
    _CHROMATIC_WIDGETS.append(cls)
    return cls

def initiate_chromatic_widgets(
        viewer : Viewer,
        wavelength_list : list[int],
        voxel_size : tuple,
) -> tuple[list, list[NapariWidget]]:
    widget_list = []
    linked_widgets = []
    for cls in _CHROMATIC_WIDGETS :
        instance = cls(voxel_size=voxel_size, wavelength_list = wavelength_list, viewer=viewer)
        if hasattr(instance,"enabled") :
            if instance.enabled :
                widget_list.extend(instance.get_widgets())
        else :
            widget_list.extend(instance.get_widgets())

        if hasattr(instance, "update") and callable(getattr(instance, "update")):
            linked_widgets.append(instance)

    return widget_list, linked_widgets

@register_chromatic_widget
class SpotCorrector(ChromaticWidget) :
    def __init__(self, *, viewer: Viewer, wavelength_list : list[int], voxel_size : tuple[int,int,int]):
        super().__init__(viewer=viewer, wavelength_list = wavelength_list, voxel_size = voxel_size)

    def _create_widget(self):
        
        @magicgui(
            reference_wavelength = {"choices" : self.wavelength_list},
            layer_wavelenth = {"choices" : self.wavelength_list}
        )
        def correct_spots(
            Spots : Points,
            reference_wavelength : int,
            layer_wavelenth : int
        ) :

            if not calibration_exist(reference_wavelength, corrected_wavelength=layer_wavelenth) :
                raise UserInputError(f"Not calibration was found for reference wavelength : {reference_wavelength}nm and layer wavelength : {layer_wavelenth}")

            calibration = load_calibration(reference_wavelength=reference_wavelength, corrected_wavelength=layer_wavelenth)
            new_coordinates = apply_polynomial_transform_spots(
                coords=Spots.data,
                poly=calibration['polynomial_features_inv'],
                model_x = calibration['x_inv_fit'],            
                model_y = calibration['y_inv_fit'],
                model_z = calibration['z_inv_fit'],       
                voxel_size= np.array(self.voxel_size, dtype=int) 
            ).round().astype(int)

            res = LayerDataTuple((
                new_coordinates,
                {
                    "name" : Spots.name,
                    "size" : Spots.size,
                    "face_color" : Spots.face_color,
                    "opacity" : Spots.opacity,
                    "blending" : Spots.blending,
                    "border_color" : Spots.border_color,
                    "symbol" : Spots.symbol,
                    "scale" : Spots.scale,
                    'units' : "nm",
                },
                'Points'
            ))

            return res

        return correct_spots

@register_chromatic_widget
class SignalCorrector(ChromaticWidget) :
    def __init__(self, *, viewer: Viewer, wavelength_list : list[int], voxel_size : tuple[int,int,int]):
        super().__init__(viewer=viewer, wavelength_list = wavelength_list, voxel_size = voxel_size)

    def _create_widget(self):
        
        @magicgui(
            reference_wavelength = {"choices" : self.wavelength_list},
            layer_wavelenth = {"choices" :self.wavelength_list}
        )
        def correct_spots(
            Signal : Image,
            reference_wavelength : int,
            layer_wavelenth : int
        ) :

            if not calibration_exist(reference_wavelength, corrected_wavelength=layer_wavelenth) :
                raise UserInputError(f"Not calibration was found for reference wavelength : {reference_wavelength}nm and layer wavelength : {layer_wavelenth}")

            calibration = load_calibration(reference_wavelength=reference_wavelength, corrected_wavelength=layer_wavelenth)

            if Signal.data.ndim == 4 :
                new_signal = np.stack(
                    [apply_polynomial_transform_to_signal(
                image=cast(np.ndarray, fov),
                poly=calibration['polynomial_features'],
                model_x = calibration['x_fit'],            
                model_y = calibration['y_fit'],
                model_z = calibration['z_fit'],       
                voxel_size= np.array(self.voxel_size, dtype=int) 
            ).round().astype(int) for fov in tqdm(cast(np.ndarray, Signal.data), desc="correcting chromatic aberrations", total=len(cast(np.ndarray, Signal.data)))]
                )
            elif Signal.data.ndim == 3 :
                new_signal = apply_polynomial_transform_to_signal(
                    image=cast(np.ndarray, Signal.data),
                    poly=calibration['polynomial_features'],
                    model_x = calibration['x_fit'],            
                    model_y = calibration['y_fit'],
                    model_z = calibration['z_fit'],       
                    voxel_size= np.array(self.voxel_size, dtype=int) 
                ).round().astype(int)

            else :
                raise AssertionError("Unforseen dimension")


            res = LayerDataTuple((
                new_signal,
                {
                    "name" : Signal.name,
                    "opacity" : Signal.opacity,
                    "blending" : Signal.blending,
                    "scale" : Signal.scale,
                    "projection_mode" : Signal.projection_mode,
                    "colormap" : Signal.colormap,
                    "contrast_limits" : Signal.contrast_limits,
                    "gamma" : Signal.gamma,
                    'units' : "nm",
                },
                'Image'
            ))

            return res

        return correct_spots

@register_chromatic_widget
class ChromaticAberrationCalibrator(ChromaticWidget) :
    def __init__(self, voxel_size : tuple, viewer : Viewer, wavelength_list):

        self.model_x = LinearRegression()
        self.model_y = LinearRegression()
        self.model_z = LinearRegression()
        self.polynomial_features = PolynomialFeatures()
        self.polynomial_features_inv = PolynomialFeatures()
        self.inv_model_x = LinearRegression()
        self.inv_model_y = LinearRegression()
        self.inv_model_z = LinearRegression()
        self.calibration_folder = CALIBRATION_FOLDER
        self.voxel_size = (1,1,1)
        self.degree = 2
        self.timestamp = get_datetime()
        self.save_widget = self._create_save_widget()
        
        super().__init__(viewer=viewer, voxel_size=voxel_size, wavelength_list=wavelength_list)

        self.register_widget(self.save_widget)

    def _create_widget(self):
        """
        Perform calibration for chromatic abberration correction and create a layer with corrected signal to evaluate quality of fit.
        """

        @magicgui(
                image_abberation={'label' : 'Image to correct :'},
                spatial_reference={'label' : 'Points reference'},
                spatial_reference_shifted={'label' : 'Points with aberrations'},
                location = {"min" : 0},
                degree={'label' : 'Degree'},
                auto_call=False,
                call_button= "Correct chromatic aberrations",
        )
        def create_corrected_layer(
            image_abberation : Image,
            spatial_reference : Points,
            spatial_reference_shifted : Points,
            location : int,
            degree : int = self.degree,
        ) :
            
            voxel_size = spatial_reference.scale
            if len(voxel_size) == 4 :
                voxel_size = voxel_size[1:]
            self.voxel_size = tuple([int(v) for v in voxel_size]) # save as reference if user save calibration

            #Convert pixel coordinates to nm to account for anisotropy
            coords1 = spatial_reference.data
            coords2 = spatial_reference_shifted.data

            if coords1.shape[1] == 4 :
                coords1 = coords1[coords1[:,0] == location]
                coords1 = coords1[:,1:]
            if coords2.shape[1] == 4 :
                coords2 = coords2[coords2[:,0] == location]
                coords2 = coords2[:,1:]
            
            coords1 = coords1 * voxel_size
            coords2 = coords2 * voxel_size


            beads, dist = match_beads(
                coords1= coords1,
                coords2= coords2,
                max_dist= int(max(voxel_size) * 4)
            )

            self.polynomial_features, self.model_x, self.model_y, self.model_z = fit_polynomial_transform_3d(
                                                beads,
                                                dist, 
                                                degree=degree
                                                )
            
            self.polynomial_features_inv, self.inv_model_x, self.inv_model_y, self.inv_model_z = fit_polynomial_transform_3d(
                                                dist, 
                                                beads,
                                                degree=degree
                                                )
            
            if image_abberation.data.ndim == 4 :
                image_corrected = np.stack(
                    [apply_polynomial_transform_to_signal(
                image=cast(np.ndarray, fov),
                poly=self.polynomial_features,
                model_x=self.model_x,
                model_y=self.model_y,
                model_z=self.model_z,
                voxel_size=voxel_size
            ).round().astype(int) for fov in tqdm(cast(np.ndarray, image_abberation.data), desc="correcting chromatic aberrations", total=len(cast(np.ndarray, image_abberation.data)))]
                )
            elif image_abberation.data.ndim == 3 :
                image_corrected = apply_polynomial_transform_to_signal(
                        image=cast(np.ndarray, image_abberation.data),
                        poly=self.polynomial_features,
                        model_x=self.model_x,
                        model_y=self.model_y,
                        model_z=self.model_z,
                        voxel_size=voxel_size
                    ).round().astype(int)
            else :
                raise AssertionError

            res = LayerDataTuple((
                image_corrected,
                {   "name" : "Interpolation result",
                    "opacity" : image_abberation.opacity,
                    "blending" : image_abberation.blending,
                    "scale" : image_abberation.scale,
                    "projection_mode" : image_abberation.projection_mode,
                    "colormap" : image_abberation.colormap,
                    'units' : "nm",
                    "contrast_limits" : image_abberation.contrast_limits,
                    "gamma" : image_abberation.gamma,},
                "Image"

            ))

            return res

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
        def save_method(
            reference_wavelength : int,
            corrected_wavelength : int,
        ) :
            
            save_fit_model(
                x_fit=self.model_x,
                y_fit=self.model_y,
                z_fit=self.model_z,
                polynomial_features= self.polynomial_features,
                polynomial_features_inv= self.polynomial_features_inv,
                x_inv_fit=self.inv_model_x,
                y_inv_fit=self.inv_model_y,
                z_inv_fit=self.inv_model_z,
                voxel_size=self.voxel_size,
                degree=self.degree,
                timestamp= self.timestamp,
                corrected_wavelength=corrected_wavelength,
                reference_wavelength=reference_wavelength,
            )
        
        return save_method