"""
Main script to call for analysis pipeline.
"""
import os
import pandas as pd
from typing import cast
import numpy as np

from ..settings import get_settings
from ..settings import AnalysisParameters

from .post_processing import Spots_filtering
from .density import density_analysis
from .distributions import distributions_analysis
from .colocalisation import main as coloc_main
from ..tools import safe_merge_no_duplicates
from ..chromatic_abberrations import correct_Spots_dataframe
from .dashboards import main as launch_dashboards


ANALYSIS_MODULES = ['all','distributions' ,'density', 'pipeline_metrics', 'pair-colocalization', 'colocalization']

def run(run_path,*args) :
    
    if '-h' in args or '--help' in args :    
        print(f"Avalaible modules are {ANALYSIS_MODULES}")
        return True
    
    if run_path is None : quit()
    else : print(run_path)

    analysis_parameters = cast(AnalysisParameters, get_settings(run_path, settings_name="analysis"))
    
    REQUIRED_TABLES = ["Acquisition", "Detection", "Spots", "Clusters", "Drift", "Gene_map", "Cell"]
    truth_table = np.array([os.path.isfile(os.path.join(run_path, "result_tables", table +".feather")) for table in REQUIRED_TABLES], dtype=bool)
    if not all(truth_table) :
        raise FileNotFoundError(f"All data tables could not be found in result directory. Missing {np.array(REQUIRED_TABLES)[~truth_table]}.")

    Acquisition = pd.read_feather(run_path + "/result_tables/Acquisition.feather")
    Detection = pd.read_feather(run_path + "/result_tables/Detection.feather")
    Spots = pd.read_feather(run_path + "/result_tables/Spots.feather")
    Drift = pd.read_feather(run_path + "/result_tables/Drift.feather")
    Gene_map = pd.read_feather(run_path + "/result_tables/Gene_map.feather")
    Cell = pd.read_feather(run_path + "/result_tables/Cell.feather")

    #Post-processing
    if not analysis_parameters.reference_wavelength is None :
        print("Correcting chromatic abberations :")
        print(f"reference wavelength : {analysis_parameters.reference_wavelength}")
        print(f"found wavelengths : {list(Detection["wavelength"].unique())}")
        Spots = correct_Spots_dataframe(#Chromatic abberation correction
            Detection=Detection,
            Spots=Spots,
            reference_wavelength= analysis_parameters.reference_wavelength
        ) 
    unfiltered_Spots = Spots.copy()
    
    #Rename target
    if not analysis_parameters.RENAME_RULE is None :
        Gene_map["target"] = Gene_map['target'].replace(analysis_parameters.RENAME_RULE)

    #Filter RNA
    if not analysis_parameters.FILTER_RNA is None :
        loc_map = Gene_map.loc[Gene_map["target"].isin(analysis_parameters.FILTER_RNA), ["cycle","color_id"]]
        filtered_detection  = Detection.loc[(Detection["cycle"].isin(loc_map["cycle"])) & (Detection["color_id"].isin(loc_map["color_id"])),["detection_id"]]
        Spots = Spots.loc[~Spots["detection_id"].isin(filtered_detection.squeeze())]

    #Filter cycles :
    if not analysis_parameters.FILTER_CYCLE is None :
        for target, cycles in analysis_parameters.FILTER_CYCLE.items() :
            loc_map = Gene_map.loc[
                (Gene_map["target"] == target) & (Gene_map["cycle"].isin(cycles)) ,
                ["cycle","color_id"]
                ]
            filtered_detection  = Detection.loc[(Detection["cycle"].isin(loc_map["cycle"])) & (Detection["color_id"].isin(loc_map["color_id"])),["detection_id"]]
            Spots = Spots.loc[~Spots["detection_id"].isin(filtered_detection.squeeze())]

    
    Spots = Spots_filtering(
        Spots,
        filter_washout=True,
        segmentation_filter=True,
        Cell=Cell,
        Detection=Detection
    )
    
    Spots_with_washout = Spots_filtering(
        unfiltered_Spots,
        filter_washout=False,
        segmentation_filter=True,
        Cell=Cell,
        Detection=Detection
    )

    Spots = _spots_merge_data(
        Spots=Spots,
        Acquisition=Acquisition,
        Detection=Detection,
        Gene_map=Gene_map,
        Cell=Cell
    )

    Spots_with_washout = _spots_merge_data(
        Spots=Spots_with_washout,
        Acquisition=Acquisition,
        Detection=Detection,
        Gene_map=Gene_map,
        Cell=Cell
    )



    if "distributions" in args or "all" in args :
        if not analysis_parameters.distribution_measures is None and len(analysis_parameters.distribution_measures) > 0:
            distribution_sucess = distributions_analysis(
                Acquisition=Acquisition,
                Detection=Detection,
                Cell=Cell,
                Spots=Spots,
                Gene_map=Gene_map,
                run_path=run_path,
                disibutions_measures= analysis_parameters.distribution_measures
            )
            if not distribution_sucess :
                print("Error raised during distribution analysis. Please check log in ~analysis/distribution_analysis folder.")
    
    if "density" in args  or "all" in args:
        
        density_sucess = density_analysis(
            Acquisition=Acquisition,
            Detection=Detection,
            Spots=Spots,
            Gene_map=Gene_map,
            run_path=run_path,
            min_number_spots=analysis_parameters.min_spots_number,
            min_diversity=analysis_parameters.min_diversity,
            cluster_radius=analysis_parameters.cluster_radius
        )
        if not density_sucess :
            print("Error raised during density analysis. Please check log in ~analysis/density_analysis folder.")

#GENERAL DATA QUALITY DASHBOARDS
    any_pipeline_metrics = any((
        "pipeline" in args,
        "pipeline_metrics" in args,
        "pipeline metrics" in args,
    ))
    if any_pipeline_metrics or "all" in args:
        sucess = launch_dashboards(
            run_path,
            Spots=Spots,
            Cell=Cell
        )
        if not sucess :
            print("Error raised during dashboards analysis. Please check log in ~analysis/dashboards/ folder.")


# COLOCALIZATION
    any_paircoloc = any((
        'coloc' in args,
        'colocalisation' in args,
        'colocalization' in args,
        'pair' in args,
        'pair-colocalisation' in args,
        'pair-colocalization' in args,
        
    ))
    if any_paircoloc or "all" in args:

        if not analysis_parameters.foci_rnas is None :
            Spots = add_foci_to_analysis(
                Spots,
                foci_rnas=analysis_parameters.foci_rnas
            )
        
        coloc_sucess = coloc_main(
            filtered_Spots=Spots,
            Cell=Cell,
            Acquisition=Acquisition,
            Detection=Detection,
            Gene_map=Gene_map,
            colocalisation_distance=analysis_parameters.coloc_distance,
            run_path=run_path,
            significance= analysis_parameters.coloc_significance,
            frameon=analysis_parameters.frameon
        )

        if not coloc_sucess :
            print(f"Error raised during coloc analysis. Please check log in ~analysis/co_localization/")


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

def add_foci_to_analysis(
    Spots : pd.DataFrame, 
    foci_rnas : list[str]
    ) :
    """
    Separate part of Spots data to analyse spots by populations : clustered and free.
    """

    save_len = len(Spots.dropna())

    foci_spots = Spots.loc[(Spots["target"].isin(foci_rnas)) & (Spots["population"] == "clustered")]
    foci_spots.loc[:,["target"]] = foci_spots["target"].str.cat(['_clustered']*len(foci_spots))
    Spots.loc[(Spots["target"].isin(foci_rnas)) & (Spots["population"] == "clustered")] = foci_spots

    free_spots = Spots.loc[(Spots["target"].isin(foci_rnas)) & (Spots["population"] == "free")]
    free_spots.loc[:,["target"]] = free_spots["target"].str.cat(['_free']*len(free_spots))
    Spots.loc[(Spots["target"].isin(foci_rnas)) & (Spots["population"] == "free")] = free_spots

    assert save_len == len(Spots.dropna())

    return Spots