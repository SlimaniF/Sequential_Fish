import pandas as pd
from typing import TypedDict, Optional

from pydantic import BaseModel, Field


class table_dict_type(TypedDict) :
    Acquisition : pd.DataFrame
    Detection : pd.DataFrame
    Spots : pd.DataFrame
    Clusters : pd.DataFrame
    Drift : pd.DataFrame
    Cell : pd.DataFrame
    Colocalisation : pd.DataFrame
    Gene_map : pd.DataFrame

class pipeline_parameters(BaseModel) :
    RUN_PATH : Optional[str] = Field(default="")
    VOXEL_SIZE : tuple[int,int,int] = Field(default=(200,97,97))
    WASHOUT_KEY_WORD : str = Field(default="Washout")
    BEAD_CHANNEl : int | None = Field(default=None)
    DAPI_CHANNEl : int  = Field(default=-1) #Default to last channel
    MAP_FILENAME : str = Field(default="cycle_file.xlsx")
    CYCLE_KEY : str = Field(default="Cycle n.")
    GENES_NAMES_KEY : list[str] = Field(default_factory= lambda : ["Gene0", ])
    cycle_regex : str = Field(default=r"img(\d+)_000_000000_0000000000.ome.tif")
    FISH_FOLDER : str = Field(default= "FISH_Z-stacks")
    MODEL_DICT : dict[str,str] = Field(default_factory=lambda : {'nucleus_model' : 'cpsam', 'cytoplasm_model' : 'cpsam'})
    OBJECT_SIZE_DICT : dict[str,int] = Field(default_factory=lambda : {'nucleus_size' : 60, 'cytoplasm_size' : 80})
    DO_3D_SEGMENTATION : bool = Field(default=False)
    PLOT_VISUALS : bool = Field(default=True)
    detection_MAX_WORKERS : int = Field(default= 4)
    SPOT_SIZE : tuple[int,int,int] = Field(default= (300,140,140))
    ALPHA : float = Field(default=0.5)
    BETA : float = Field(default=1)
    GAMMA : int = Field(default=2)
    CLUSTER_SIZE : int = Field(default=400)
    MIN_SPOT_PER_CLUSTER : int = Field(default=5)
    ARTIFACT_RADIUS : int = Field(default=1400)
    DETECTION_SLICE_TO_REMOVE : list[int | None] = Field(default_factory=lambda: [5,None])
    SAVE_PATH : str = Field(default="/analysis/")
    BEAD_SIZE : tuple[int,int,int] = Field(default= (200,200,200))
    DRIFT_SLICE_TO_REMOVE : list[int | None] = Field(default_factory=lambda:[5,5])
    DO_HIGHPASS_FILTER : bool = Field(default= False)
    REFERENCE_CYCLE : int = Field(default=0)
    COLOC_DISTANCE : int = Field(default=200)
    quantif_MAX_WORKERS : int = Field(default=10)

    @classmethod
    def from_default_parameters(cls) :
        """Create an instance with default parameters"""
        return cls()