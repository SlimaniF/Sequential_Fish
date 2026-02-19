"""
Submodule for data post processing, eg Filtering...
"""

import pandas as pd
from .analysis_parameters import FILTER_RNA
from ..tools import safe_merge_no_duplicates

def Spots_filtering(
    Spots : pd.DataFrame,
    Detection : pd.DataFrame = None,
    Cell : pd.DataFrame = None,
    filter_washout= True,
    segmentation_filter= True,
    ) :
    
    if filter_washout : 
        Spots = Spots.loc[~Spots['is_washout']]
    
    if segmentation_filter :
        Spots = Spots.loc[Spots['cell_label'] != 0]
        # Spots = Spots.loc[] #Create couple(location, label) and try if spots couple are cell couple.
    
    
    if (not Cell is None) and (not Detection is None) :
        
        if 'location' not in Spots.columns :
            Spots = safe_merge_no_duplicates(
                Spots,
                Detection,
                keys='location',
                on='detection_id'
            )
        
        Spots = pd.merge(
            Spots,
            Cell,
            how='inner',
            left_on= ['location','cell_label','detection_id'],
            right_on= ['location','label','detection_id']
        )
        
    return Spots

def RNA_filtering(df_with_target : pd.DataFrame, rna_to_filter : 'list[str]' = FILTER_RNA) :
    
    if 'target' not in df_with_target : raise KeyError('"target" column was not found in dataframe columns.')
    
    df_with_target = df_with_target.drop(
        df_with_target[df_with_target['target'].isin(rna_to_filter)].index,
        axis=0
    )
    
    return df_with_target