import os
import logging

import pandas as pd
from matplotlib.pyplot import close

from .cell_quality import cell_dashboard
from .signal_quality import signal_quality_dashboard

def main(
    run_path : str,
    Spots : pd.DataFrame,
    Cell : pd.DataFrame,
    Acquisition : pd.DataFrame,
    Gene_map : pd.DataFrame,
    Drift : pd.DataFrame,
    Detection : pd.DataFrame,
    coloc_range : int,
    drift_checker : tuple[str,str],
    chroma_checker : tuple[str,str]
    ) :


    logging.info("Generating cell data board.")
    cell_figure = cell_dashboard(
        run_path=run_path,
        Cell=Cell,
        Spots=Spots
    )
    save_folder = os.path.join(run_path,"analysis")
    os.makedirs(save_folder,exist_ok=True)
    cell_figure.savefig(os.path.join(save_folder,"cell_dashboard.svg"))
    close(cell_figure)

    logging.info("Generating signal quality data board.")
    signal_figure = signal_quality_dashboard(
        run_path=run_path,
        Acquisition=Acquisition,
        Cell=Cell,
        Gene_map=Gene_map,
        Detection=Detection,
        Drift=Drift,
        coloc_range=coloc_range,
        drift_checker=drift_checker,
        chroma_checker=chroma_checker
    )
    if not signal_figure is None :
        signal_figure.savefig(os.path.join(save_folder,"signal_dashboard.svg"))
    close(signal_figure)

    return True