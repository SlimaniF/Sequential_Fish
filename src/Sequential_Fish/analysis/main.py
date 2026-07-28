"""
Main script to call for analysis pipeline.
"""
import os
import logging
import pandas as pd
from typing import cast
import numpy as np

from ..settings import get_settings
from ..settings import AnalysisParameters

from .post_processing import Spots_post_processing, apply_user_configuration
from .density import density_analysis
from .distributions import distributions_analysis
from .colocalisation import main as coloc_main
from .dashboards import main as launch_dashboards


ANALYSIS_MODULES = ['all','distributions' ,'density', 'pipeline_metrics', 'pair-colocalization', 'colocalization']

def run(run_path,*args) :
    
    if '-h' in args or '--help' in args :    
        print(f"Avalaible modules are {ANALYSIS_MODULES}")
        return True

    #Init
    analysis_parameters = cast(AnalysisParameters, get_settings(run_path, settings_name="analysis"))
    REQUIRED_TABLES = ["Acquisition", "Detection", "Spots", "Clusters", "Drift", "Gene_map", "Cell"]
    truth_table = np.array([os.path.isfile(os.path.join(run_path, "result_tables", table +".feather")) for table in REQUIRED_TABLES], dtype=bool)
    if not all(truth_table) :
        raise FileNotFoundError(f"All data tables could not be found in result directory. Missing {np.array(REQUIRED_TABLES)[~truth_table]}.")

    log_path = os.path.join(run_path,"analysis","analysis_log.log")
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
    )

    #Loading tables
    Acquisition = pd.read_feather(run_path + "/result_tables/Acquisition.feather")
    Detection = pd.read_feather(run_path + "/result_tables/Detection.feather")
    Spots = pd.read_feather(run_path + "/result_tables/Spots.feather")
    Drift = pd.read_feather(run_path + "/result_tables/Drift.feather")
    Gene_map = pd.read_feather(run_path + "/result_tables/Gene_map.feather")
    Cell = pd.read_feather(run_path + "/result_tables/Cell.feather")

    #User defined filters
    Gene_map, Detection, Spots = apply_user_configuration(
        Gene_map=Gene_map,
        Detection=Detection,
        Spots=Spots,
        rename_rule=analysis_parameters.RENAME_RULE,
        filter_cycle=analysis_parameters.FILTER_CYCLE,
        filter_rna=analysis_parameters.FILTER_RNA
    )

    #Post-processing
    Spots = Spots_post_processing(
        Spots=Spots,
        Cell=Cell,
        Detection=Detection,
        Acquisition=Acquisition,
        Gene_map=Gene_map,
        reference_wavelength = analysis_parameters.reference_wavelength
    )

    # Call to analysis submodules
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
            Cell=Cell,
            Acquisition=Acquisition,
            Gene_map=Gene_map,
            Drift=Drift,
            Detection=Detection,
            coloc_range=analysis_parameters.coloc_distance,
            drift_checker= analysis_parameters.drift_checker,
            chroma_checker=analysis_parameters.chroma_checker
        )
        if not sucess :
            print("Error raised during dashboards analysis. Please check log in ~analysis/dashboards/ folder.")


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