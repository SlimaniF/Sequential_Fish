import napari, os
import numpy as np
import joblib
from abc import ABC, abstractmethod
from magicgui import magicgui
from magicgui import widgets
from napari.layers import Points, Image, Layer
from typing import Tuple
from bigfish.detection import detect_spots

from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from scipy.ndimage import map_coordinates

from sklearn.neighbors import NearestNeighbors
def match_beads(coords1, coords2, max_dist=5):
    """Match nearest beads between two channels."""
    nn = NearestNeighbors(n_neighbors=1).fit(coords2)
    distances, indices = nn.kneighbors(coords1)
    matches = distances[:, 0] < max_dist
    return coords1[matches], coords2[indices[matches, 0]]

def fit_polynomial_transform_3d(src_points, dst_points, degree=2):
    """Fit 3D polynomial regression mapping coords → dst."""
    poly = PolynomialFeatures(degree)
    X_poly = poly.fit_transform(src_points)
    model_x = LinearRegression().fit(X_poly, dst_points[:, 2])  # x
    model_y = LinearRegression().fit(X_poly, dst_points[:, 1])  # y
    model_z = LinearRegression().fit(X_poly, dst_points[:, 0])  # z
    return poly, model_x, model_y, model_z

def apply_polynomial_transform_3d(
        image : np.array, 
        poly, 
        model_x, 
        model_y, 
        model_z, 
        voxel_size : np.array
        ):
    """Warp 3D image using learned polynomial transform."""
    z, y, x = image.shape
    zz, yy, xx = np.meshgrid(np.arange(z), np.arange(y), np.arange(x), indexing='ij')
    coords = np.stack([zz.ravel(), yy.ravel(), xx.ravel()], axis=1)

    X_poly = poly.transform(coords * voxel_size)
    new_z_nm = model_z.predict(X_poly)
    new_y_nm = model_y.predict(X_poly)
    new_x_nm = model_x.predict(X_poly)

    #convert back to pixel
    if voxel_size.ndim == 1 : voxel_size = np.array([voxel_size])
    new_coords_pixel = np.stack([new_z_nm, new_y_nm, new_x_nm], axis=0) / voxel_size.T

    warped = map_coordinates(image, new_coords_pixel, order=1, mode='reflect').reshape(z, y, x)
    return warped


class NapariWidget(ABC) :
    """
    Common super class for custom widgets added to napari interface during run
    Each sub class as a specific function, but the widget can be acess with attribute .widget
    """
    def __init__(self):
        self.widget = self._create_widget()

    @abstractmethod
    def _create_widget(self) :
        """
        This should return a widget you can add to the napari (QWidget)
        """
        pass

class BeadsDetector(NapariWidget) :
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
    
class ChromaticAberrationCorector(NapariWidget) :
    def __init__(self, save_path = os.getcwd()):
        super().__init__()

        self.model_x = LinearRegression()
        self.model_y = LinearRegression()
        self.model_z = LinearRegression()
        self.inv_model_x = LinearRegression()
        self.inv_model_y = LinearRegression()
        self.inv_model_z = LinearRegression()
        self.save_path = save_path
        self.voxel_size = (1,1,1)

    def _create_widget(self):

        @magicgui(
                image_abberation={'label' : 'Image to correct :'},
                spatial_reference_shifted={'label' : 'Points with aberrations'},
                spatial_reference={'label' : 'Points reference'},
                auto_call=False,
                call_button= "Correct chromatic aberrations",
        )
        def create_corrected_layer(
            image_abberation : Image,
            spatial_reference_shifted : Points,
            spatial_reference : Points,
            degree = 2,
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
            

            image_corrected = apply_polynomial_transform_3d(
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

        return create_corrected_layer
    
    def _create_save_widget(self) :

        @magicgui(
                auto_call=False, 
                call_button= "Save calibration"
                
                )
        def save_fit_model() :
            joblib.dump({
                'x_fit' : self.model_x,
                'y_fit' : self.model_y,
                'z_fit' : self.model_z,
                'x_inv_fit' : self.inv_model_x,
                'y_inv_fit' : self.inv_model_y,
                'z_inv_fit' : self.inv_model_z,
                'voxel_size' : self.voxel_size,
                'degree' : self.degree,
                'reference_wavelength' : self.reference_wavelength,
                'corrected_wavelength' : self.corrected_wavelength,
                'timestamp' : None, #TODO
            },
            filename=self.save_path + "/chromatic"
            )

    

# main_test
from Sequential_Fish.tools.utils import open_image

image_location = "/media/SSD_floricslimani/Fish_seq/Davide/Chromatic abberations/"

image = open_image(image_location + "/Beads-Field_01.tif")
image.shape

channels = [0,1,2]
colors = ['green', 'red', 'blue']

os.environ["QT_QPA_PLATFORM"] = "xcb"
Viewer = napari.Viewer()

for chan, color in zip(channels, colors) :
    Viewer.add_image(image[...,chan], name= f"Channel_{chan}", scale=(200,97,97), blending='additive', colormap = color, interpolation2d= 'cubic')

#Widgets
beads_detector = BeadsDetector()
abberration_corrector = ChromaticAberrationCorector()
right_container = widgets.Container(widgets=[beads_detector.widget, abberration_corrector.widget], labels=False)
Viewer.window.add_dock_widget(right_container, name='Beads detector', area='right')

napari.run()
