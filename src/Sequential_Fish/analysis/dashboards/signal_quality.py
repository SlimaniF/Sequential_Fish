"""
This dashboard is dedicated to bleaching and washout and maybe later background removal
"""

from matplotlib.axes import Axes
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from Sequential_Fish.tools import safe_merge_no_duplicates
import matplotlib.gridspec as gridspec

def plot_intensity_bleach_curves(
    ax : Axes,
    Cell : pd.DataFrame,
    color_id : int,
    mean_intensity_per_cell : pd.DataFrame
    ) -> Axes:

    if color_id in Cell["color_id"].unique() :
        target = "fish_mean_mean_signal"
    else :
        target = "nucleus_mean_mean_signal"
        color_id = 0

    view = Cell.loc[(~Cell["target"].str.contains("Washout",case=False)) & (Cell["color_id"] == color_id),["color_id","location","cycle",target]]

    print("color_id : ",color_id)
    print(view)

    sns.lineplot(
        view,
        x="cycle",
        y=target,
        hue="color_id",
        ax=ax
    )

    return ax

def _get_mean_intensity_per_cell(
    Cell,
    ) :
    fish_results = Cell.groupby(["color_id","location","cycle"])["fish_mean_mean_signal"].apply(tuple).rename("mean signal")
    dapi_results : pd.DataFrame = Cell.groupby(["color_id","location","cycle"])["nucleus_mean_mean_signal"].apply(list).rename("mean signal").reset_index(drop=False)
    dapi_results["color_id"] = dapi_results["color_id"].max() +1 
    dapi_results["mean signal"] = dapi_results["mean signal"].apply(tuple)
    dapi_results = dapi_results.drop_duplicates()

    return  pd.concat([fish_results,dapi_results], axis=0, ignore_index=True)



if __name__ == "__main__" :
    run_path = "/media/SSD_floricslimani/Fish_seq/POLR2/2026-03-24 - HeLa_POLR2_Run17_tiff/"

    Acquisition = pd.read_feather(run_path +"/result_tables/Acquisition.feather")
    Cell = pd.read_feather(run_path +"/result_tables/Cell.feather")
    Gene_map = pd.read_feather(run_path +"/result_tables/Gene_map.feather")
    Detection = pd.read_feather(run_path +"/result_tables/Detection.feather")

    Detection = safe_merge_no_duplicates(
        Detection,
        Acquisition,
        on="acquisition_id",
        keys=["cycle"]
    )

    Detection = safe_merge_no_duplicates(
        Detection,
        Gene_map,
        on=["cycle","color_id"],
        keys=["target"]
    )
    Cell = safe_merge_no_duplicates(
        Cell,
        Detection,
        on="detection_id",
        keys=["target","cycle","color_id"]
    )

    mean_intensity_per_cell = _get_mean_intensity_per_cell(Cell)
    color_number = Detection["color_id"].nunique() +1

    fig = plt.figure(figsize=(40,20))
    grid = gridspec.GridSpec(3,color_number)

    for color in range(color_number) :
        ax = fig.add_subplot(grid[0,color])
        ax  = plot_intensity_bleach_curves(
            ax= ax,
            Cell=Cell,
            color_id=color,
            mean_intensity_per_cell=mean_intensity_per_cell
            )
    

    plt.show()