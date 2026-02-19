"""
Main script to call for analysis pipeline.
"""
import pandas as pd

from .post_processing import Spots_filtering
from .density import density_analysis
from .distributions import distributions_analysis
from ..run_saves import select_path

ANALYSIS_MODULES = ['all','distributions' ,'density', 'pipeline_metrics', 'pair-colocalization', 'colocalization']

def run(*args) :
    
    if '-h' in args or '--help' in args :    
        print(f"Avalaible modules are {ANALYSIS_MODULES}")
        return True
    
    run_path = select_path()
    if run_path is None : quit()
    else : print(run_path)
    
    Acquisition = pd.read_feather(run_path + "/result_tables/Acquisition.feather")
    Detection = pd.read_feather(run_path + "/result_tables/Detection.feather")
    Spots = pd.read_feather(run_path + "/result_tables/Spots.feather")
    Drift = pd.read_feather(run_path + "/result_tables/Drift.feather")
    Gene_map = pd.read_feather(run_path + "/result_tables/Gene_map.feather")
    Cell = pd.read_feather(run_path + "/result_tables/Cell.feather")

    #Post-processing
    unfiltered_Spots = Spots.copy()
    
    #Rename target
    from .analysis_parameters import RENAME_RULE, frameon
    Gene_map["target"] = Gene_map['target'].replace(RENAME_RULE)
    
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
    
    if "distributions" in args or "all" in args :
        
        from .analysis_parameters import distribution_measures
        
        distribution_sucess = distributions_analysis(
            Acquisition=Acquisition,
            Detection=Detection,
            Cell=Cell,
            Spots=Spots,
            Gene_map=Gene_map,
            run_path=run_path,
            disibutions_measures= distribution_measures
        )
        if not distribution_sucess :
            print("Error raised during distribution analysis. Please check log in ~analysis/distribution_analysis folder.")
    
    if "density" in args  or "all" in args:
        
        from .analysis_parameters import min_diversity, min_spots_number, cluster_radius
        density_sucess = density_analysis(
            Acquisition=Acquisition,
            Detection=Detection,
            Spots=Spots,
            Gene_map=Gene_map,
            run_path=run_path,
            min_number_spots=min_spots_number,
            min_diversity=min_diversity,
            cluster_radius=cluster_radius
        )
        if not density_sucess :
            print("Error raised during density analysis. Please check log in ~analysis/density_analysis folder.")
        
    any_pipeline_metrics = any((
        "pipeline" in args,
        "pipeline_metrics" in args,
        "pipeline metrics" in args,
    ))
    if any_pipeline_metrics or "all" in args:
        from .pipeline_metrics import pipeline_metrics
        drift_sucess = pipeline_metrics(
            Acquisition=Acquisition,
            Detection=Detection,
            Gene_map= Gene_map,
            Spots_with_washout=Spots_with_washout,
            Unfiltered_spots=unfiltered_Spots,
            Cell=Cell,
            Drift=Drift,
            run_path= run_path
        )
        if not drift_sucess :
            print("Error raised during Drift analysis. Please check log in ~analysis/pipeline_metrics/ folder.")

    any_paircoloc = any((
        'coloc' in args,
        'colocalisation' in args,
        'colocalization' in args,
        'pair' in args,
        'pair-colocalisation' in args,
        'pair-colocalization' in args,
        
    ))
    if any_paircoloc or "all" in args:
        from .colocalisation import main as coloc_main
        from .analysis_parameters import coloc_distance, coloc_significance
        
        coloc_sucess = coloc_main(
            filtered_Spots=Spots,
            Cell=Cell,
            Acquisition=Acquisition,
            Detection=Detection,
            Gene_map=Gene_map,
            colocalisation_distance=coloc_distance,
            run_path=run_path,
            significance= coloc_significance,
            frameon=frameon
        )

        if not coloc_sucess :
            print(f"Error raised during coloc analysis. Please check log in ~analysis/co_localization/")