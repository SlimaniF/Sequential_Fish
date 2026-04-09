from numpy import ndarray
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from typing import TypedDict
import pandas as pd

class table_dict_type(TypedDict) :
    Acquisition : pd.DataFrame
    Detection : pd.DataFrame
    Spots : pd.DataFrame
    Clusters : pd.DataFrame
    Drift : pd.DataFrame
    Cell : pd.DataFrame
    Colocalisation : pd.DataFrame
    Gene_map : pd.DataFrame

class Calibration(TypedDict) :
    x_fit : LinearRegression
    y_fit : LinearRegression
    z_fit : LinearRegression
    x_inv_fit : LinearRegression
    y_inv_fit : LinearRegression
    z_inv_fit : LinearRegression
    polynomial_features : PolynomialFeatures
    polynomial_features_inv : PolynomialFeatures
    voxel_size : ndarray
    degree : int
    reference_wavelength : int
    corrected_wavelength : int
    timestamp : str
