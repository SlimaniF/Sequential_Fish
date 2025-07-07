import numpy as np
from scipy.ndimage import map_coordinates
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures

def apply_polynomial_transform_to_signal(
        image : np.ndarray, 
        poly : PolynomialFeatures, 
        model_x : LinearRegression, 
        model_y : LinearRegression, 
        voxel_size : np.ndarray,
        model_z : LinearRegression = None, 
        ):
    """Warp 3D image using learned polynomial transform."""
    z, y, x = image.shape
    zz, yy, xx = np.meshgrid(np.arange(z), np.arange(y), np.arange(x), indexing='ij')
    coords = np.stack([zz.ravel(), yy.ravel(), xx.ravel()], axis=1)

    if isinstance(voxel_size, tuple) : voxel_size = np.array(voxel_size)

    X_poly = poly.transform(coords * voxel_size)
    if not model_z is None : 
        new_z_nm = model_z.predict(X_poly)
    else :
        print(coords.shape)
        new_z_nm = coords[:,0] #TODO test
    new_y_nm = model_y.predict(X_poly)
    new_x_nm = model_x.predict(X_poly)

    #convert back to pixel
    if voxel_size.ndim == 1 : voxel_size = np.array([voxel_size])
    new_coords_pixel = np.stack([new_z_nm, new_y_nm, new_x_nm], axis=0) / voxel_size.T

    warped = map_coordinates(image, new_coords_pixel, order=1, mode='reflect').reshape(z, y, x)
    return warped

def apply_polynomial_transform_spots(
        coords : np.ndarray,
        poly : PolynomialFeatures, 
        model_x : LinearRegression, 
        model_y : LinearRegression, 
        voxel_size : np.ndarray,
        model_z : LinearRegression = None,
        ) :
    """
    Correct chromatic abberrations for spots using pre-calibrated polynomial interpolation.
    """


    monosomes = poly.transform(coords * voxel_size)
    new_y_nm = model_y.predict(monosomes)
    new_x_nm = model_x.predict(monosomes)
    if not model_z is None :
        new_z_nm = model_z.predict(monosomes)
    else :
        new_z_nm = coords[:,0]
    new_coords_pixel = np.stack([new_z_nm, new_y_nm, new_x_nm], axis=1) / voxel_size

    return new_coords_pixel

def get_polynomial_features(degree : int) :
    poly = PolynomialFeatures(degree)
    return poly