from numpy import ndarray
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from typing import TypedDict
from typing import Optional, List, Dict
from pydantic import BaseModel, Field

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
    
    #PATH to run folder, images and location to save
    RUN_PATH : Optional[str] = Field(default="")
    FISH_FOLDER : str = Field(default= "FISH_Z-stacks")
    SAVE_PATH : str = Field(default="/analysis/")
    
    #Microscope parameters
    VOXEL_SIZE : tuple[int,int,int] = Field(default=(200,97,97))
    BEAD_CHANNEl : int | None = Field(default=None)
    DAPI_CHANNEl : int  = Field(default=-1) #Default to last channel
    WAVELENGTH_LIST : list[int] | None = Field(default=None)

    #Folder organization
    WASHOUT_KEY_WORD : str = Field(default="Washout")
    MAP_FILENAME : str = Field(default="cycle_file.xlsx")
    GENES_NAMES_KEY : list[str] = Field(default_factory= lambda : ["Gene0", ])
    CYCLE_KEY : str = Field(default="Cycle n.")
    cycle_regex : str = Field(default=r"img(\d+)_000_000000_0000000000.ome.tif")

    #Segmentation parameters
    MODEL_DICT : dict[str,str] = Field(default_factory=lambda : {'nucleus_model' : 'cpsam', 'cytoplasm_model' : 'cpsam'})
    OBJECT_SIZE_DICT : dict[str,int] = Field(default_factory=lambda : {'nucleus_size' : 60, 'cytoplasm_size' : 80})
    DO_3D_SEGMENTATION : bool = Field(default=False)
    FLOW_THRESHOLD  : float = Field(default=0.4) #Not used in 3D
    CELLPROB_THRESHOLD : float = Field(default=0.0)
    MIN_SIZE : int = Field(default=15)
    PLOT_VISUALS : bool = Field(default=True)
    
    #Detection parameters
    detection_MAX_WORKERS : int = Field(default= 4)
    SPOT_SIZE : tuple[int,int,int] = Field(default= (300,140,140))
    ALPHA : float = Field(default=0.5)
    BETA : float = Field(default=1)
    GAMMA : int = Field(default=2)
    CLUSTER_SIZE : int = Field(default=400)
    MIN_SPOT_PER_CLUSTER : int = Field(default=5)
    ARTIFACT_RADIUS : int = Field(default=1400)
    DETECTION_SLICE_TO_REMOVE : list[int | None] = Field(default_factory=lambda: [5,None])
    
    #Drift parameters
    BEAD_SIZE : tuple[int,int,int] = Field(default= (200,200,200))
    DRIFT_SLICE_TO_REMOVE : list[int | None] = Field(default_factory=lambda:[5,5])
    REFERENCE_CYCLE : int = Field(default=0)
    DO_HIGHPASS_FILTER : bool = Field(default= False)
    
    #Quantification parameters
    COLOC_DISTANCE : int = Field(default=200)
    quantif_MAX_WORKERS : int = Field(default=10)

    @classmethod
    def from_default_parameters(cls) :
        """Create an instance with default parameters"""
        return cls()

    def get_filename(self) :
        return "pipeline_settings.json"

class AnalysisParameters(BaseModel) :
    """
    Fixed structure of data : type is inforced through pydantic BaseModel class.
    """

    #Plots
    frameon : bool = Field(default=True)
    
    #Preprocessing
    FILTER_RNA : List[str] | None = Field(default=None)
    RENAME_RULE : Dict[str,str] | None = Field(default=None)

    #Distributions
    distribution_measures : List[str] | None = Field(default=None)

    #Chromatic abberration
    reference_wavelength : int | None = Field(default=555) #None to ignore chromatic abberations correction

    #Density analysis
    min_diversity : int = Field(default=2)
    min_spots_number : int = Field(default=2)
    cluster_radius : int = Field(default=0)

    #Co-localization analysis
    coloc_distance : int = Field(default=0)
    coloc_significance : float = Field(default=10e-3)

    @classmethod
    def from_default_parameters(cls) :
        """Create an instance with default parameters"""
        return cls()

    def get_filename(self) :
        return "analysis_settings.json"