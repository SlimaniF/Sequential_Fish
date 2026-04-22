"""
This dashboard is dedicated to bleaching and washout and maybe later background removal
"""
from math import log10
import os
from typing import Literal, TypedDict, cast

import numpy as np
import pandas as pd

from skimage.metrics import structural_similarity as ssim
from skimage.metrics import mean_squared_error
from scipy.stats import pearsonr
from pebble import ProcessPool, ThreadPool

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap

from Sequential_Fish.tools import safe_merge_no_duplicates

def plot_intensity_bleach_curves(
    ax : Axes,
    Cell : pd.DataFrame,
    color_id : int,
    ) -> Axes:

    if color_id in Cell["color_id"].unique() :
        target = "fish_mean_mean_signal"
    else :
        target = "nucleus_mean_mean_signal"
        color_id = 0

    view = Cell.loc[(~Cell["target"].str.contains("Washout",case=False)) & (Cell["color_id"] == color_id),["color_id","wavelength","location","cycle",target]]
    if target == "fish_mean_mean_signal":
        wavelength = view["wavelength"].iat[0]
        ax.set_title(f"Channel {color_id} - {wavelength} nm")
    else :
        wavelength = 461
        ax.set_title(f"Nucleus {color_id} - {wavelength} nm")

    color_list = [(228,26,28),(255,255,51),(255,127,0),(77,175,74),(166,206,227),(55,126,184),(152,78,163)]
    for a,color in enumerate(color_list) : 
        color_list[a] = tuple(np.array(color)/255)

    cmap = ListedColormap(color_list)

    color = cmap(
        -(wavelength - 750)/400
    )

    sns.lineplot(
        view,
        x="cycle",
        y=target,
        color="black", err_kws={"facecolor": color},
        ax=ax
    )


    ax.set_ylabel("Signal (a.u)",size = 10)
    sns.despine(ax=ax)
    ax.set_xlim(0)

    return ax

def plot_endcycle_coloc_rates(
    ax : Axes,
    drift_checker : float,
    chroma_checker : float,
    ) -> Axes :

    data = pd.DataFrame({
        "type" : ["Drift","Abberations"],
        "colocalization rate" : [drift_checker, chroma_checker]
    })
    data["colocalization rate"] *= 100
    print(data)

    sns.barplot(
        data,
        x="type",
        y="colocalization rate",
        hue="type", palette=["#fb9a99", "#80b1d3"],
        edgecolor = "black",
        ax=ax
    )

    ax.set_ylim(0,100)
    ax.get_yaxis().set_visible(False)
    sns.despine(ax=ax, left=True)
    ax.set_xlabel("")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30)

    drift_txt = ax.text(0, data.iloc[0].at["colocalization rate"] /2, f"{round(data.iloc[0].at["colocalization rate"],1)} %")
    abb_txt = ax.text(1, data.iloc[1].at["colocalization rate"] /2, f"{(round(data.iloc[1].at["colocalization rate"],1))} %")

    plt.setp([drift_txt,abb_txt], size=10, fontweight="bold", color="white", rotation = 15, va="center", ha="center")

    return ax

def _get_coloc_rates_for_checkers(
    run_path : str,
    **checkers : dict[str,tuple[str,str]]
    ) :

    table_path = os.path.join(run_path,"analysis","co_localization","datasheet","coloc_rates_mean.csv")
    print(table_path)
    if not os.path.isfile(table_path) : raise FileNotFoundError("Couldn't find coloc rates table, run first co-localization analysis")
    coloc_table = pd.read_csv(table_path, sep=";").set_index("target")

    coloc_rates = {key : float(coloc_table.at[str(pair[0]),str(pair[1])].replace(",",'.')) for key,pair in checkers.items()}

    return coloc_rates

class SimilarityMetrics(TypedDict) :
    spatial_correlation : float
    structural_similarity : float
    mean_squared_error : float
    peak_snr : float

def compute_similarity_metrics(
    image1 : np.ndarray,
    image2 : np.ndarray
    ) :

    if np.issubdtype(image1.dtype, np.floating) :
        max_intensity = 1
    else :
        max_intensity = np.iinfo(image1.dtype).max

    res = {
        "spatial_correlation" : pearsonr(image1.flatten(), image2.flatten()),
        "structural_similarity" : ssim(image1, image2, full=True),
        "mean_squared_error" : mean_squared_error(image1, image2),
    }
    res["peak_snr"] = 10 * log10(max_intensity**2 / res["mean_squared_error"])

    return SimilarityMetrics(**res)


class DriftMeans(TypedDict) :
    z_mean_drift : pd.Series
    y_mean_drift : pd.Series
    x_mean_drift : pd.Series
    euclidian_mean_drift : pd.Series


def plot_drift_means(
    ax : Axes,
    drift_means : DriftMeans,
    measure : Literal["z_mean_drift", "y_mean_drift","x_mean_drift","euclidian_mean_drift",]
    ) :

    data = pd.DataFrame(drift_means).reset_index(drop=False)

    sns.swarmplot(
        data = data,
        x = measure,
        orient="h",
        edgecolor='black',
        linewidth=1,
        color='white',
        ax=ax
    )

    ax.set_xlabel("")
    ax.set_ylabel(measure.split("_")[0], size=10)
    ax.set_xlim(0)

    return ax


def _compute_drift_means(
    Acquisition : pd.DataFrame,
    Drift : pd.DataFrame
    ) :

    Drift = safe_merge_no_duplicates(
        Drift,
        Acquisition,
        on="acquisition_id",
        keys="location"
    )
    Drift["euclidian_drift"] = (Drift["drift_z"].pow(2) + Drift["drift_y"].pow(2) +Drift["drift_x"].pow(2)).apply(np.sqrt)
    return DriftMeans(
        z_mean_drift = cast(pd.Series, Drift.groupby("location")["drift_z"].mean()),
        y_mean_drift = cast(pd.Series, Drift.groupby("location")["drift_y"].mean()),
        x_mean_drift = cast(pd.Series, Drift.groupby("location")["drift_x"].mean()),
        euclidian_mean_drift = cast(pd.Series, Drift.groupby("location")["euclidian_drift"].mean())
        )


def process_similarity_metrics(
    run_path : str,
    Gene_map : pd.DataFrame,
    Acquisition : pd.DataFrame,
    drift_namepair : tuple[str,str],
    abb_namepair : tuple[str,str]
    ) :
    
    pass #TODO

if __name__ == "__main__" :
    run_path = "/home/floric/Documents/SeqfFish/2026-03-24 - HeLa_POLR2_Run17_tiff/"

    Acquisition = pd.read_feather(run_path +"/result_tables/Acquisition.feather")
    Cell = pd.read_feather(run_path +"/result_tables/Cell.feather")
    Gene_map = pd.read_feather(run_path +"/result_tables/Gene_map.feather")
    Detection = pd.read_feather(run_path +"/result_tables/Detection.feather")
    Drift = pd.read_feather(run_path +"/result_tables/Drift.feather")

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
        keys=["target","cycle","color_id","wavelength"]
    )
    Cell["fish_mean_mean_signal"] = Cell["nucleus_mean_mean_signal"]

    #1.
    color_number = Detection["color_id"].nunique() +1

    #2.
    drift_checker = ('POLR2A','POLR2A_endcycle')
    chroma_checker = ('PRPF8','PRPF8_endcycle')
    coloc_rates = _get_coloc_rates_for_checkers(
        run_path=run_path,
        drift_checker=drift_checker,
        chroma_checker=chroma_checker
    )

    #3. Similarity measurements

    #4. Drift
    drift_means = _compute_drift_means(
        Acquisition=Acquisition,
        Drift=Drift
    )


    fig = plt.figure(figsize=(40,25))
    grid = gridspec.GridSpec(5,color_number +2, 
    width_ratios=[1,1] + [2]*color_number,
    height_ratios=[1]*4 + [2],
    wspace=0.5
    )

    #Plot bleach curves
    for color in range(color_number) :
        ax = fig.add_subplot(grid[:2,color+2])
        ax  = plot_intensity_bleach_curves(
            ax= ax,
            Cell=Cell,
            color_id=color,
            )
        
        if color > 0 :
            ax.set_ylabel("")
            sns.despine(ax=ax,left=True)
    
    coloc_checkers_ax = fig.add_subplot(grid[4,0])
    similarity_ax = fig.add_subplot(grid[4,1])

    #Plot drift swarmplots
    max_drift = max([axis.max() for axis in drift_means.values()])
    for i, measure in enumerate(["z_mean_drift","y_mean_drift","x_mean_drift","euclidian_mean_drift"]) :
        ax = fig.add_subplot(grid[i,:2])
        ax = plot_drift_means(
            ax=ax,
            drift_means=drift_means,
            measure=measure
        )

        if i == 0 :
            ax.set_title("Drift (px)")
        ax.set_xlim(0,max_drift+0.5)

        if i < 3 :
            ax.set_xticks([])

    plot_endcycle_coloc_rates(
        ax=coloc_checkers_ax,
        **coloc_rates
    )


    plt.show()