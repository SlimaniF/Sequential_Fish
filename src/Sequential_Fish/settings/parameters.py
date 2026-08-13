
from pydantic import BaseModel, Field
import warnings

from pydantic.fields import FieldInfo

class ParametersModel(BaseModel) :
    """
    Subclass this class to create Pydantic data model enforcing type when file is loaded and to be used for dynamic QForm GUI allowing user to modify paramters.
    Subclasses need to define each class variable with default value to Field with default argument or default_factory set.
    In Field info pass json_schema_extra with a dict containing 'tab' key to organize parameters in different tabs in GUI.
    """

    @classmethod
    def from_default_parameters(cls) :
        """Create an instance with default parameters"""
        return cls()
    
    @classmethod
    def get_tab_names(cls) -> list[str]:
        """Return all names found in json_schema_extra with key 'tab'."""

        tab_names = set()
        for field_name, field_info in cls.model_fields.items() :
            json_schema_extra = field_info.json_schema_extra
            if json_schema_extra is None :
                warnings.warn(f"{field_name} was not assigned to a tab in Field description. Pass a 'json_schema_extra' with 'tab' to assign it to a tab. Ignore this warning to assign to 'other' tab.")
                tab_names.add("other")
            elif isinstance(json_schema_extra, dict) :
                tab_names.add(json_schema_extra.get("tab","other"))
            else :
                raise AssertionError("json_schema_extra not None nor dict")
        
        return sorted(tab_names)

    @classmethod
    def get_parameters_from_tab(cls, tab_name : str) -> dict[str, FieldInfo]:

        tab_parameters = {}
        for field_name, field_info in cls.model_fields.items() :
            json_schema_extra = field_info.json_schema_extra
            if json_schema_extra is None :
                found_tab_name = "other"
            elif isinstance(json_schema_extra, dict) :
                found_tab_name = json_schema_extra.get("tab","other")
            else :
                raise AssertionError("json_schema_extra not None nor dict")

            if tab_name == found_tab_name :
                tab_parameters[field_name] = field_info

        return tab_parameters


class PipelineParameters(ParametersModel) :
    
    #PATH to run folder, images and location to save
    FISH_FOLDER : str = Field(default= "FISH_Z-stacks", json_schema_extra={"tab" : "general"})
    LOCATION_KEYWORD : str = Field(default="Location-", json_schema_extra={"tab" : "general"})
    
    #Microscope parameters
    VOXEL_SIZE : tuple[int,int,int] = Field(default=(200,97,97), json_schema_extra={"tab" : "general"})
    BEAD_CHANNEl : int | None = Field(default=None, json_schema_extra={"tab" : "general"})
    DAPI_CHANNEl : int  = Field(default=-1, json_schema_extra={"tab" : "general"}) #Default to last channel
    WAVELENGTH_LIST : list[int] | None = Field(default=None, json_schema_extra={"tab" : "general"})

    #Folder organization
    WASHOUT_KEY_WORD : str = Field(default="Washout", json_schema_extra={"tab" : "general"})
    MAP_FILENAME : str = Field(default="cycle_file.xlsx", json_schema_extra={"tab" : "general"})
    GENES_NAMES_KEY : list[str] = Field(default_factory= lambda : ["Gene0", ], json_schema_extra={"tab" : "general"})
    CYCLE_KEY : str = Field(default="Cycle n.", json_schema_extra={"tab" : "general"})
    cycle_regex : str = Field(default=r"img(\d+)_000_000000_0000000000.ome.tif", json_schema_extra={"tab" : "general"})

    #Segmentation parameters
    DO_3D_SEGMENTATION_NUCLEUS : bool = Field(default=False, json_schema_extra={"tab" : "segmentation"})
    DO_3D_SEGMENTATION_CYTOPLASM : bool = Field(default=False, json_schema_extra={"tab" : "segmentation"})
    SEGMENT_ONLY_NUCLEI : bool = Field(default=False, json_schema_extra={"tab" : "segmentation"})
    MODEL_DICT : dict[str,str] = Field(default_factory=lambda : {'nucleus_model' : 'cpsam', 'cytoplasm_model' : 'cpsam'}, json_schema_extra={"tab" : "segmentation"})
    OBJECT_SIZE_DICT : dict[str,int] = Field(default_factory=lambda : {'nucleus_size' : 60, 'cytoplasm_size' : 80}, json_schema_extra={"tab" : "segmentation"})
    FLOW_3D_SMOOTH : dict[str,int] = Field(default_factory=lambda : {'nucleus' : 0, 'cytoplasm' : 0}, json_schema_extra={"tab" : "segmentation"})
    FLOW_THRESHOLD  : dict[str,float] = Field(default_factory=lambda : {'nucleus' : 0.4, 'cytoplasm' : 0.4}, json_schema_extra={"tab" : "segmentation"}) #Not used in 3D
    CELLPROB_THRESHOLD : dict[str,float] = Field(default_factory=lambda : {'nucleus' : 0., 'cytoplasm' : 0.}, json_schema_extra={"tab" : "segmentation"})
    MIN_SIZE : dict[str,int] = Field(default_factory=lambda : {'nucleus' : 15, 'cytoplasm' : 15}, json_schema_extra={"tab" : "segmentation"})
    PLOT_VISUALS : bool = Field(default=True, json_schema_extra={"tab" : "segmentation"})
    
    #Detection parameters
    detection_MAX_WORKERS : int = Field(default= 4, json_schema_extra={"tab" : "detection"})
    SPOT_SIZE : tuple[int,int,int] = Field(default= (300,140,140), json_schema_extra={"tab" : "detection"})
    ALPHA : float = Field(default=0.5, json_schema_extra={"tab" : "detection"})
    BETA : float = Field(default=1, json_schema_extra={"tab" : "detection"})
    GAMMA : int = Field(default=2, json_schema_extra={"tab" : "detection"})
    CLUSTER_SIZE : int = Field(default=400, json_schema_extra={"tab" : "detection"})
    MIN_SPOT_PER_CLUSTER : int = Field(default=5, json_schema_extra={"tab" : "detection"})
    ARTIFACT_RADIUS : int = Field(default=1400, json_schema_extra={"tab" : "detection"})
    DETECTION_SLICE_TO_REMOVE : list[int | None] = Field(default_factory=lambda: [5,None], json_schema_extra={"tab" : "detection"})
    
    #Drift parameters
    BEAD_SIZE : tuple[int,int,int] = Field(default= (200,200,200), json_schema_extra={"tab" : "drift"})
    DRIFT_SLICE_TO_REMOVE : list[int | None] = Field(default_factory=lambda:[5,5], json_schema_extra={"tab" : "drift"})
    REFERENCE_CYCLE : int = Field(default=0, json_schema_extra={"tab" : "drift"})
    DO_HIGHPASS_FILTER : bool = Field(default= False, json_schema_extra={"tab" : "drift"})
    
    #Quantification parameters
    COLOC_DISTANCE : int = Field(default=200, json_schema_extra={"tab" : "quantification"})
    quantif_MAX_WORKERS : int = Field(default=10, json_schema_extra={"tab" : "quantification"})

    @classmethod
    def get_filename(cls) :
        return "pipeline_settings.json"

class AnalysisParameters(ParametersModel) :
    """
    Fixed structure of data : type is inforced through pydantic BaseModel class.
    """

    #Plots
    frameon : bool = Field(default=True, json_schema_extra={"tab" : "General"})
    
    #Preprocessing
    FILTER_RNA : list[str] | None = Field(default=None, json_schema_extra={"tab" : "General"})
    FILTER_CYCLE : dict[str,list[int]] | None = Field(default=None, json_schema_extra={"tab" : "General"})
    RENAME_RULE : dict[str,str] | None = Field(default=None, json_schema_extra={"tab" : "General"})

    #Distributions
    distribution_measures : list[str] | None = Field(default=None, json_schema_extra={"tab" : "Distribution"})

    #Chromatic abberration
    reference_wavelength : int | None = Field(default=555, json_schema_extra={"tab" : "ChromaticAbberations"}) #None to ignore chromatic abberations correction

    #Density analysis
    min_diversity : int = Field(default=2, json_schema_extra={"tab" : "Density"})
    min_spots_number : int = Field(default=2, json_schema_extra={"tab" : "Density"})
    cluster_radius : int = Field(default=0, json_schema_extra={"tab" : "Density"})

    #Co-localization analysis
    coloc_distance : int = Field(default=0, json_schema_extra={"tab" : "Colocalization"})
    coloc_significance : float = Field(default=10e-3, json_schema_extra={"tab" : "Colocalization"})
    foci_rnas : list[str] | None = Field(default=None, json_schema_extra={"tab" : "Colocalization"})

    #Dashboard
    drift_checker : tuple[str,str] = Field(default=("",""), json_schema_extra={"tab" : "General"})
    chroma_checker : tuple[str,str] = Field(default=("",""), json_schema_extra={"tab" : "General"})

    #Multivariate_exploratory
    control_genes : list[str] | None = Field(default=None, json_schema_extra={"tab" : "Multivariate"})

    @classmethod
    def get_filename(cls) :
        return "analysis_settings.json"