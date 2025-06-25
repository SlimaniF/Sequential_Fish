from numpy import ndarray
from sklearn.linear_model import LinearRegression
from typing import TypedDict

class Calibration(TypedDict) :
    x_fit : LinearRegression
    y_fit : LinearRegression
    z_fit : LinearRegression
    x_inv_fit : LinearRegression
    y_inv_fit : LinearRegression
    z_inv_fit : LinearRegression
    voxel_size : ndarray
    degree : int
    reference_wavelength : int
    corrected_wavelength : int
    timestamp : str