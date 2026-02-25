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
