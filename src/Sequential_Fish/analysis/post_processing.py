"""
Submodule for data post processing : cycle/rna filtering, data merging, preparing cached data...
"""

import os
import logging

import pandas as pd
from ..tools import safe_merge_no_duplicates
from ..chromatic_abberrations import correct_Spots_dataframe
from ..settings import AnalysisParameters
from .colocalisation import colocalisation_truth_df


def Spots_post_processing(
        Spots : pd.DataFrame,
        Cell : pd.DataFrame,
        Detection : pd.DataFrame,
        Acquisition : pd.DataFrame,
        Gene_map : pd.DataFrame,
        reference_wavelength : int | None = None
) :
    """
    Filter Spots in washout artifacts and Spots in cells found on segmentation edges.
    """
    logging.info("Starting post-processing operations") 

    Spots = _Spots_filtering(
        Spots,
        filter_washout=True,
        segmentation_filter=True,
        Cell=Cell,
        Detection=Detection,
    )

    Spots = _spots_merge_data(
        Spots=Spots,
        Acquisition=Acquisition,
        Detection=Detection,
        Gene_map=Gene_map,
        Cell=Cell
    )


    if not reference_wavelength is None :
        logging.info(
            f"Correcting chromatic abberations :\n\t -> reference wavelength : {reference_wavelength}\n\t -> found wavelengths : {list(Detection["wavelength"].unique())}"
            )
    
        Spots = correct_Spots_dataframe(
            Detection=Detection,
            Spots=Spots,
            reference_wavelength= reference_wavelength
        ) 

    return Spots

def _Spots_filtering(
    Spots : pd.DataFrame,
    Detection : pd.DataFrame = None,
    Cell : pd.DataFrame = None,
    filter_washout= True,
    segmentation_filter= True,
    ) :
    """
    Filters :
        -> Washout
        -> Spots out of segmentation
        -> Spots in discarded cells
    """
    
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
            Cell.loc[:,['location','label','detection_id']],
            how='inner',
            left_on= ['location','cell_label','detection_id'],
            right_on= ['location','label','detection_id']
        )
        
    return Spots

def _spots_merge_data(
    Spots : pd.DataFrame,
    Acquisition : pd.DataFrame,
    Detection : pd.DataFrame,
    Gene_map : pd.DataFrame,
    Cell : pd.DataFrame,
    ) :
    """
    Merge required information into Spots df
    """
    
    Detection = safe_merge_no_duplicates(
        Detection,
        Acquisition,
        on= ['acquisition_id'],
        keys=['cycle','location', 'fish_reodered_shape']
    )

    Detection = safe_merge_no_duplicates(
        Detection,
        Gene_map,
        on= ['cycle','color_id'],
        keys=['target']
    )

    Spots =safe_merge_no_duplicates(
        Spots,
        Detection,
        on= 'detection_id',
        keys= ['location','target', 'voxel_size', 'fish_reodered_shape']
    )

    Spots = safe_merge_no_duplicates(
        Spots,
        Cell.rename(columns={'label' : 'cell_label'}),
        on=['acquisition_id','detection_id','cell_label'],
        keys=['cell_id']
    )


    return Spots

def apply_user_configuration(
        Gene_map : pd.DataFrame,
        Detection : pd.DataFrame,
        Spots : pd.DataFrame,
        rename_rule,
        filter_rna,
        filter_cycle
) :
    """
    Apply user configuration : renaming targets, rna filter, cycle filter.

    # Returns
        Gene_map, Detection, Spots

    """

    logging.info(f"Renaming targets using rule : {rename_rule}")
    Gene_map = _rename_targets(
        Gene_map,
        rule=rename_rule
        )

    logging.info(f"Removing RNAs from analysis : {filter_rna}")
    Detection, Spots = _remove_rna(
        Gene_map=Gene_map,
        filter_rna=filter_rna,
        Detection=Detection,
        Spots=Spots
    )

    logging.info(f"Removing cycles from analysis using rule : {filter_cycle}")
    Detection, Spots = _remove_cycles(
        Gene_map=Gene_map,
        filter_cycles=filter_cycle,
        Detection=Detection,
        Spots=Spots
    )

    return Gene_map, Detection, Spots


def _rename_targets(
        Gene_map : pd.DataFrame,
        rule : dict | None = None
) :
    if not rule is None :
        Gene_map["target"] = Gene_map['target'].replace(rule)

    return Gene_map

def _remove_rna(
        Gene_map : pd.DataFrame,
        filter_rna : list[str] | None,
        Detection : pd.DataFrame,
        Spots : pd.DataFrame
) :
    if filter_rna is None :
        return Detection, Spots

    loc_map = Gene_map.loc[~Gene_map["target"].isin(filter_rna), ["cycle","color_id"]]
    filtered_detection  = Detection.loc[(Detection["cycle"].isin(loc_map["cycle"])) & (Detection["color_id"].isin(loc_map["color_id"]))]
    filtered_spots = Spots.loc[Spots["detection_id"].isin(filtered_detection["detection_id"])]

    return filtered_detection, filtered_spots

def _remove_cycles(
        Gene_map : pd.DataFrame,
        filter_cycles : dict[str,list[int]] | None,
        Detection : pd.DataFrame,
        Spots : pd.DataFrame,        
) :

    if filter_cycles is None :
        return Detection, Spots
    
    for target, cycles in filter_cycles.items() :
        loc_map = Gene_map.loc[
            (Gene_map["target"] == target) & (Gene_map["cycle"].isin(cycles)) ,
            ["cycle","color_id"]
            ]
        filtered_detection  = Detection.loc[(Detection["cycle"].isin(loc_map["cycle"])) & (Detection["color_id"].isin(loc_map["color_id"])),["detection_id"]]
        filtered_spots = Spots.loc[~Spots["detection_id"].isin(filtered_detection.squeeze())]

    return filtered_detection, filtered_spots


def update_cache(
        run_path : str,
        analysis_parameters : AnalysisParameters
) :

    _cache_colocalization_data(run_path, analysis_parameters)

def _cache_colocalization_data(
        run_path : str, 
        analysis_parameters :  AnalysisParameters
        ) :

    run_coloc_truth_table = False

    cached_data_path = os.path.join(run_path,"result_tables", "coloc_truth_table.feather")
    cache_attr = ["coloc_distance", "FILTER_CYCLE", "RENAME_RULE"] # cached data contains all run RNAs on purpose they are filtered when loading the table and using it for figures.
    if os.path.isfile(cached_data_path) : 
        df = pd.read_feather(cached_data_path, columns=["spot_id"])
        attrs : dict = df.attrs

        for key in cache_attr : 
            assert key in attrs.keys(), f"{key} was not found in cache metadata."
            if attrs[key] != getattr(analysis_parameters,key) : run_coloc_truth_table = True

    else :
        run_coloc_truth_table = True

    if run_coloc_truth_table :
        logging.info("Update in parameters, computing co-localization truth df.")
        logging.disable(logging.INFO)
        df = _create_cache_colocalization_truth_table(
            run_path, 
            colocalization_distance=analysis_parameters.coloc_distance, 
            cycle_filter_rule=analysis_parameters.FILTER_CYCLE,
            rename_rule=analysis_parameters.RENAME_RULE,
            reference_wv=analysis_parameters.reference_wavelength
            )

        for key in cache_attr :
            df.attrs[key] = getattr(analysis_parameters,key)
        df.to_feather(cached_data_path)
        logging.disable(logging.NOTSET)
        logging.info("cache update completed.")
        

def _create_cache_colocalization_truth_table(
        run_path : str, 
        colocalization_distance : int, 
        cycle_filter_rule : dict,
        rename_rule : dict,
        reference_wv : int
) -> pd.DataFrame :

    result_path = os.path.join(run_path,"result_tables")
    Acquisition = pd.read_feather(os.path.join(result_path,"Acquisition.feather"))
    Detection = pd.read_feather(os.path.join(result_path,"Detection.feather"))
    Gene_map = pd.read_feather(os.path.join(result_path,"Gene_map.feather"))
    Spots = pd.read_feather(os.path.join(result_path,"Spots.feather"))
    Cell = pd.read_feather(os.path.join(result_path,"Cell.feather"))

    Gene_map, Detection, Spots = apply_user_configuration(
        Gene_map=Gene_map,
        Detection=Detection,
        Spots=Spots,
        rename_rule=rename_rule,
        filter_cycle=cycle_filter_rule,
        filter_rna=None
    )

    Spots = Spots_post_processing(
        Spots=Spots,
        Cell=Cell,
        Detection=Detection,
        Acquisition=Acquisition,
        Gene_map=Gene_map,
        reference_wavelength=reference_wv
    )

    Spots.info()

    coloc_truth_df= colocalisation_truth_df(
        Spots=Spots,
        colocalisation_distance=colocalization_distance
    )

    coloc_truth_df.info()

    return coloc_truth_df