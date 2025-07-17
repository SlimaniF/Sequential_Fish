from numpy import ndarray
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from typing import TypedDict
from typing import Tuple, List, Dict, Union
from pydantic import BaseModel #Runtime type check

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

class PipelineParameters(BaseModel) :
    """
    Fixed structure of data : type is inforced through pydantic BaseModel
    """
    
    #Global
    VOXEL_SIZE : Tuple[int,int,int]
    WASHOUT_KEY_WORD : str
    HAS_BEAD_CHANNEL : bool

    #Input
    MAP_FILENAME  : str
    CYCLE_KEY : str
    GENES_NAMES_KEY : List[str]
    WAVELENGTH_LIST : List[int]
    cycle_regex : str
    FOLDER_KEYS : Dict[str,str]

    #Segmentation
    MODEL_DICT : Dict[str,str]
    OBJECT_SIZE_DICT : Dict[str,int]

    #Detection
    detection_MAX_WORKERS : int
    SPOT_SIZE : Tuple[int,int,int]
    ALPHA : float
    BETA : float
    GAMMA : float
    CLUSTER_SIZE : int
    MIN_SPOT_PER_CLUSTER : int
    ARTIFACT_RADIUS : int
    DETECTION_SLICE_TO_REMOVE : List[Union[int, None]]

    #Drift
    BEAD_SIZE : Tuple[int,int,int]
    DRIFT_SLICE_TO_REMOVE : List[Union[int,None]]

    #Quantification
    COLOC_DISTANCE : int
    quantif_MAX_WORKERS : int