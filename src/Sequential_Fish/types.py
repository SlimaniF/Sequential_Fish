import pandas as pd
from typing import TypedDict

from pydantic import BaseModel


class table_dict_type(TypedDict) :
    Acquisition : pd.DataFrame
    Detection : pd.DataFrame
    Spots : pd.DataFrame
    Clusters : pd.DataFrame
    Drift : pd.DataFrame
    Cell : pd.DataFrame
    Colocalisation : pd.DataFrame
    Gene_map : pd.DataFrame

class parameters_dict(BaseModel) :
    RUN_PATH : str
    VOXEL_SIZE : tuple[int,int,int]
    WASHOUT_KEY_WORD : str
    HAS_BEAD_CHANNEL : bool
    MAP_FILENAME : str
    CYCLE_KEY : str
    GENES_NAMES_KEY : list[str]
    cycle_regex : str
    FOLDER_KEYS : dict[str,str]
    MODEL_DICT : dict[str,str]
    OBJECT_SIZE_DICT : dict[str,int]
    PLOT_VISUALS : bool
    detection_MAX_WORKERS : int
    SPOT_SIZE : tuple[int,int,int]
    ALPHA : float
    BETA : float
    GAMMA : int
    CLUSTER_SIZE : int
    MIN_SPOT_PER_CLUSTER : int
    ARTIFACT_RADIUS : int
    DETECTION_SLICE_TO_REMOVE : list[int | None]
    SAVE_PATH : str
    BEAD_SIZE : tuple[int,int,int]
    DRIFT_SLICE_TO_REMOVE : list[int | None]
    DO_HIGHPASS_FILTER : bool
    COLOC_DISTANCE : int
    quantif_MAX_WORKERS : int
