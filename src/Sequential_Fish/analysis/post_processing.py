"""
Submodule for data post processing : cycle/rna filtering, data merging, preparing cached data...
"""

import os
import logging

import pandas as pd
import numpy as np
from ..tools import safe_merge_no_duplicates
from ..chromatic_abberrations import correct_Spots_dataframe
from ..settings import AnalysisParameters
from .colocalisation import colocalisation_truth_df, compute_z_score_frame, create_coloc_rate_expectancy, _compute_cell_distribution_populations


def Spots_post_processing(
        Spots : pd.DataFrame,
        Cell : pd.DataFrame,
        Detection : pd.DataFrame,
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


    return Spots, Detection

def apply_user_configuration(
        Gene_map : pd.DataFrame,
        Detection : pd.DataFrame,
        Spots : pd.DataFrame,
        Cell : pd.DataFrame,
        Acquisition : pd.DataFrame,
        rename_rule,
        filter_rna,
        filter_cycle,
        foci_rnas,
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

    Spots,Detection = _spots_merge_data(
        Spots=Spots,
        Acquisition=Acquisition,
        Detection=Detection,
        Gene_map=Gene_map,
        Cell=Cell
    )

    if not foci_rnas is None :
        logging.info(f"Adding free and clustered populations for rnas : {foci_rnas}")

        Spots = _add_foci_to_analysis(
            Spots,
            foci_rnas=foci_rnas
        )

    logging.info(f"Removing RNAs from analysis : {filter_rna}")

    if not filter_rna is None :
        Detection = Detection.loc[~Detection['target'].isin(filter_rna)]
        Spots = Spots.loc[~Spots['target'].isin(filter_rna)]

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
    cache_attr = ["coloc_distance", "FILTER_CYCLE", "RENAME_RULE","foci_rnas"] # cached data contains all run RNAs on purpose they are filtered when loading the table and using it for figures.
    if os.path.isfile(cached_data_path) : 
        coloc_truth_df = pd.read_feather(cached_data_path, columns=["spot_id"])
        attrs : dict = coloc_truth_df.attrs

        for key in cache_attr : 
            assert key in attrs.keys(), f"{key} was not found in cache metadata."
            if attrs[key] != getattr(analysis_parameters,key) : run_coloc_truth_table = True

    else :
        run_coloc_truth_table = True

    if run_coloc_truth_table :
        logging.info("Update in parameters, computing co-localization truth coloc_truth_df.")
        logging.disable(logging.INFO)

        result_path = os.path.join(run_path,"result_tables")
        Acquisition = pd.read_feather(os.path.join(result_path,"Acquisition.feather"))
        Detection = pd.read_feather(os.path.join(result_path,"Detection.feather"))
        Gene_map = pd.read_feather(os.path.join(result_path,"Gene_map.feather"))
        Spots = pd.read_feather(os.path.join(result_path,"Spots.feather"))
        Cell = pd.read_feather(os.path.join(result_path,"Cell.feather"))

        Spots = Spots_post_processing(
            Spots=Spots,
            Cell=Cell,
            Detection=Detection,
            reference_wavelength=analysis_parameters.reference_wavelength
        )

        Gene_map, Detection, Spots = apply_user_configuration(
            Gene_map=Gene_map,
            Detection=Detection,
            Spots=Spots,
            Cell=Cell,
            Acquisition=Acquisition,
            rename_rule=analysis_parameters.RENAME_RULE,
            filter_cycle=analysis_parameters.FILTER_CYCLE,
            filter_rna=None,
            foci_rnas=analysis_parameters.foci_rnas
        )

        coloc_truth_df = _create_cache_colocalization(
            Spots=Spots,
            colocalization_distance=analysis_parameters.coloc_distance, 
            )

        zscore_df = _create_zscore_table(
            Spots=Spots,
            Cell=Cell,
            coloc_truth_df=coloc_truth_df,
            voxel_size=Detection["voxel_size"].at[0],
            colocalisation_distance=analysis_parameters.coloc_distance
        )

        for key in cache_attr :
            coloc_truth_df.attrs[key] = getattr(analysis_parameters,key)
        coloc_truth_df.to_feather(cached_data_path)

        zscore_path = os.path.join(run_path,"analysis","data","coloc_zscores.csv")
        os.makedirs(os.path.dirname(zscore_path),exist_ok=True)
        zscore_df.to_csv(zscore_path, sep=";")

        logging.disable(logging.NOTSET)
        logging.info("cache update completed.")
        

def _create_cache_colocalization(
        Spots : pd.DataFrame,
        colocalization_distance : int, 
) -> pd.DataFrame :

    coloc_truth_df= colocalisation_truth_df(
        Spots=Spots,
        colocalisation_distance=colocalization_distance
    )

    coloc_truth_df = safe_merge_no_duplicates(
        coloc_truth_df,
        Spots,
        on="spot_id",
        keys="cell_id"
    )

    return coloc_truth_df


def _create_zscore_table(
        Spots : pd.DataFrame,
        Cell : pd.DataFrame,
        coloc_truth_df : pd.DataFrame,
        voxel_size : tuple,
        colocalisation_distance : int,
) :

    rna_list = Spots["target"].unique().tolist()

    coloc_rates, selfcoloc_rates = create_coloc_rate_expectancy(
        Spots=Spots,
        Cell=Cell,
        voxel_size=voxel_size,
        colocalisation_distance=colocalisation_distance,
        RNA_list=rna_list
    )

    # Model values for expectancy and std
    abundancies = _compute_cell_distribution_populations(Spots)
    cell_ids = coloc_rates.index
    multi_index = pd.MultiIndex.from_product([rna_list, abundancies.index]) #one line per couple (rna,cell_id)
    multi_index.names = ['target','cell_id']

    ## Initialise empty dataframes for expectancy/std
    expected_event_count = pd.DataFrame(columns= pd.Index(rna_list), index=multi_index, dtype=float)
    expected_event_count_std = pd.DataFrame(columns= pd.Index(rna_list), index=multi_index, dtype=float)

    ## Filling
    for rna in rna_list :
        prod = coloc_rates.multiply(abundancies[rna],axis=0)
        product_index = pd.MultiIndex.from_product([[rna], cell_ids])

        ###Expectancy
        expected_event_count.loc[product_index, :] = prod.values # E = n*p
        expected_event_count.loc[rna,[rna]] = (selfcoloc_rates[rna] * abundancies[rna]).values # Correction for selfcoloc

        ###Std
        prod = coloc_rates.multiply(abundancies[rna],axis=0).multiply((1-coloc_rates),axis=0) #std = sqrt(np(1-p))
        prod = prod.apply(np.sqrt)
        product_index = pd.MultiIndex.from_product([[rna], cell_ids])
        expected_event_count_std.loc[product_index, :] = prod.values

    # Coloc measurements
    coloc_truth_df = pd.merge( #Filter spots belonging to RNA distributions removed in user configuration
        coloc_truth_df,
        Spots.loc[:,["spot_id"]],
        on='spot_id',
    )

    measure_coloc_events = coloc_truth_df.groupby(['target','cell_id'])[rna_list].sum()
    assert not measure_coloc_events.isna().any().any(), "Analysis shouldn't yield nan at this point. Make sure that that previously nan values are safe to discard and discard them prior to this point." #Then safe to ignore nan values as nan values comes from 0 abundancies, ie cell without spot which should be discarded when comptuting co-localization statistic.
    coloc_rates = measure_coloc_events.divide(abundancies).replace(np.inf,np.nan).replace(-np.inf,np.nan)
    
    #Zscore computation
    zscore_frame = compute_z_score_frame(
        measured_colocalisation_events=measure_coloc_events,
        expected_colocalisation_events=expected_event_count,
        expected_standard_deviation= expected_event_count_std,
    )

    return zscore_frame


def _add_foci_to_analysis(
    Spots : pd.DataFrame, 
    foci_rnas : list[str]
    ) :
    """
    Separate part of Spots data to analyse spots by populations : clustered and free.
    """

    foci_spots = Spots.loc[(Spots["target"].isin(foci_rnas)) & (Spots["population"] == "clustered")]
    foci_spots.loc[:,["target"]] = foci_spots["target"].str.cat(['_clustered']*len(foci_spots))
    free_spots = Spots.loc[(Spots["target"].isin(foci_rnas)) & (Spots["population"] == "free")]
    free_spots.loc[:,["target"]] = free_spots["target"].str.cat(['_free']*len(free_spots))

    free_spots["spot_id"] = np.arange(1, len(free_spots)+1) + Spots["spot_id"].max()
    foci_spots["spot_id"] = np.arange(1, len(foci_spots)+1) + free_spots["spot_id"].max()

    Spots = pd.concat([
        Spots,
        foci_spots,
        free_spots
    ], ignore_index=True)



    return Spots