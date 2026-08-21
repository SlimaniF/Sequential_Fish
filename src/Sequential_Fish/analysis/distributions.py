import os
import logging, traceback, warnings
import pandas as pd
import matplotlib.pyplot as plt

from .utils import distribution_super_plot
from .utils import get_xlabels

def distributions_analysis(
    Detection : pd.DataFrame,
    Cell : pd.DataFrame,
    distributions_measures : 'list[str]',
    run_path :str,
    washout_keyword : str | None = None,
) :
    
    output_path = run_path + "/analysis/distribution_analysis/"
    os.makedirs(output_path, exist_ok=True)
    
    try :
        logging.info(f"Starting distributions plotting for following metrics : {distributions_measures}") 

        if not washout_keyword is None :
            Detection = Detection.loc[~Detection["target"].str.contains(washout_keyword)]

        Cell=pd.merge(
            Cell,
            Detection.loc[:,["detection_id","target","cycle"]],
            on="detection_id",
        )
        
        for measure in distributions_measures :

            if not measure in Cell.columns : 
                logging.error(f"{measure} was not found in Cell columns.")
                warnings.warn(f"{measure} was not found in Cell columns.")
                continue
        
            data = Cell.groupby('target')[measure].apply(list)
            fig = generate_distribution_graph(
                data=data,
                Cell=Cell,
                measure=measure,
            )
            fig.savefig(output_path + f"/{measure}.svg")


    except Exception :
        logging.error(f"Distributions plotting failed :\n{traceback.format_exc()}\n")
        
        
        return False
        
    else :
        logging.info("Distribution plotting success\n")
        
        return True


def generate_distribution_graph(
        data : pd.Series,
        measure : str,
        Cell : pd.DataFrame,
) :

    fig = plt.figure(figsize=(16,8))
    ax = fig.gca()
    ax = distribution_super_plot(
        data,
        ax,
        ylabel=measure,
        title= f"Distribution of {measure} per cell",
    )

    if 'index' in measure :
        min_x,max_x,*_ = plt.axis()
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

    return fig