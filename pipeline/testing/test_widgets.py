import napari, os
import numpy as np
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

def apply_polynomial_transform_3d(volume, poly, model_x, model_y, model_z, voxel_size):
    """Warp 3D volume using learned polynomial transform."""
    z, y, x = volume.shape
    zz, yy, xx = np.meshgrid(np.arange(z), np.arange(y), np.arange(x), indexing='ij')
    coords = np.stack([zz.ravel(), yy.ravel(), xx.ravel()], axis=1)

    X_poly = poly.transform(coords * voxel_size)
    new_z_nm = model_z.predict(X_poly)
    new_y_nm = model_y.predict(X_poly)
    new_x_nm = model_x.predict(X_poly)

    #convert back to pixel
    voxel_size = np.array([voxel_size])
    new_coords_pixel = np.stack([new_z_nm, new_y_nm, new_x_nm], axis=0) / voxel_size.T

    warped = map_coordinates(volume, new_coords_pixel, order=1, mode='reflect').reshape(z, y, x)
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
    def __init__(self):
        super().__init__()

    def _create_widget(self):

        @magicgui(
            beads_image = {'label' : 'Image'},
            threshold = {'label' : 'Threshold', 'min': 1},
            voxel_size = {'annotation' : Tuple[int,int,int], 'label' : 'Voxel size (zyx)'},
            beads_radius = {'annotation' : Tuple[int,int,int], 'label' : 'Beads radius (zyx)'},
            auto_call=False,
            call_button="Detect beads"
                
        )
        def detect_beads(
            beads_image : Image,
            threshold : int = 490,
            beads_radius =(300,150,150),
            voxel_size = (200,97,97),
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

            return detected_beads
        
        
        return detect_beads
    
class ChromaticAberrationCorector(NapariWidget) :
    def __init__(self):
        super().__init__()

    def _create_widget(self):

        @magicgui(
                image_abberation={'label' : 'Image to correct :'},
                spatial_reference={'label' : 'Points reference'},
                spatial_reference_shifted={'label' : 'Points with aberrations'},
                voxel_size = {'annotation' :Tuple[int,int,int], 'label' : "Voxel size (zyx)"},
                auto_call=False,
                call_button= "Correct chromatic aberrations",
        )
        def create_corrected_layer(
            image_abberation : Image,
            spatial_reference : Points,
            spatial_reference_shifted : Points,
            voxel_size : Tuple[int,int,int],
            degree = 2,
        ) ->  Image :
            
            beads, dist = match_beads(
                coords1= spatial_reference.data,
                coords2= spatial_reference_shifted.data,
                max_dist= int(max(voxel_size) * 4)
            )

            poly, model_x, model_y, model_z = fit_polynomial_transform_3d(
                                                dist, 
                                                beads,
                                                degree=degree
                                                )
            
            image_corrected = apply_polynomial_transform_3d(
                image_abberation.data,
                poly=poly,
                model_x=model_x,
                model_y=model_y,
                model_z=model_z,
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
right_container = widgets.Container(widgets=[beads_detector.widget], labels=False)
left_container = widgets.Container(widgets=[abberration_corrector.widget], labels=False)
Viewer.window.add_dock_widget(right_container, name='Beads detector', area='right')
Viewer.window.add_dock_widget(left_container, name='Aberration corrector', area='right')

napari.run()
