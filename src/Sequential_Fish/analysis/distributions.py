import os
import logging, traceback
import pandas as pd
import matplotlib.pyplot as plt

from ..tools import safe_merge_no_duplicates
from .utils import distribution_super_plot, merge_data
from .post_processing import RNA_filtering
from .utils import get_xlabels

def distributions_analysis(
    Acquisition : pd.DataFrame,
    Detection : pd.DataFrame,
    Cell : pd.DataFrame,
    Spots : pd.DataFrame,
    Gene_map : pd.DataFrame,
    disibutions_measures : 'list[str]',
    run_path :str,
) :
    output_path = run_path + "/analysis/distribution_analysis/"
    os.makedirs(output_path, exist_ok=True)
    
    log_file = output_path + "/distribution_analysis_log.log"
    logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    force= True
)
    
    try :
        print("Starting distribution analysis...")
        logging.info(f"New density analysis")
        logging.info(f"Distribution_measures :\n{disibutions_measures}")
        
        Detection, Cell, Spots = merge_data(
            Acquisition=Acquisition,
            Detection=Detection,
            Cell=Cell,
            Spots=Spots,
            Gene_map=Gene_map
        )
        
        Cell = safe_merge_no_duplicates(
            Cell,
            Detection,
            on='detection_id',
            keys='cycle'
        )
        Cell = Cell.loc[~Cell['target'].str.contains('Washout')]
        Cell = RNA_filtering(Cell)

        for measure in disibutions_measures :
        
            data = Cell.groupby('target')[measure].apply(list)

            fig = plt.figure(figsize=(16,8))
            ax = fig.gca()
            ax = distribution_super_plot(
                data,
                ax,
                ylabel=measure,
                title= f"Distribution of {measure} per cell",
            )

            if 'index' in measure :
                min_x,max_x,min_y,max_y = plt.axis()
                ax.plot([min_x, max_x], [1,1], '--b')
            
            xlabels = get_xlabels(ax)
            cycles = Cell.groupby('target')['cycle'].first()
            cell_numbers = Cell.groupby('target')['cell_id'].count()
            for i, label_data in enumerate(zip(xlabels.copy(), cycles, cell_numbers)) :
                label, cycle, cell_number = label_data
                xlabels[i] = f"{label}\n({cycle})\n[{cell_number}]"
            
            if len(xlabels) > 15 :
                ax.set_xticklabels(xlabels, rotation=30)
            else :
                ax.set_xticklabels(xlabels, rotation=0)

            plt.savefig(output_path + f"/{measure}.svg")
            plt.close()
    
    except Exception as e :
        logging.error(f"analysis failed :\n{traceback.format_exc()}")
        
        
        return False
        
    else :
        logging.info(f"analysis succeed")
        
        return True
        