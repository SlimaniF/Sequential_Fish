import os

import pandas as pd
from matplotlib.pyplot import close

from .cell_quality import cell_dashboard

def main(
    run_path : str,
    Spots : pd.DataFrame,
    Cell : pd.DataFrame
    ) :

    cell_figure = cell_dashboard(
        run_path=run_path,
        Cell=Cell,
        Spots=Spots
    )
    save_folder = os.path.join(run_path,"analysis")
    os.makedirs(save_folder,exist_ok=True)
    cell_figure.savefig(os.path.join(save_folder,"cell_dashboard.svg"))
    close(cell_figure)

    return True