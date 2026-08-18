"""
This dashboard is dedicated to bleaching and washout and maybe later background removal
"""
from math import log10
import os
from typing import Literal, TypedDict, cast

from cycler import K
import numpy as np
import pandas as pd
import warnings
import logging

from skimage.metrics import mean_squared_error
from scipy.stats import pearsonr
from pebble import ProcessPool, ThreadPool

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap

from Sequential_Fish.tools import safe_merge_no_duplicates, open_cycle
import time as t
from tqdm import tqdm
import psutil

from Sequential_Fish.tools.utils import open_location

class DriftMeans(TypedDict) :
    z_mean_drift : pd.Series
    y_mean_drift : pd.Series
    x_mean_drift : pd.Series
    euclidian_mean_drift : pd.Series

class SimilarityMetrics(TypedDict) :
    spatial_correlation : float
    structural_similarity : float
    mean_squared_error : float
    peak_snr : float


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


    ax.set_ylabel("Signal (a.u)",size = 15)
    sns.despine(ax=ax)
    ax.set_xlim(0)

    return ax

def plot_endcycle_coloc_rates(
    ax : Axes,
    coloc_range : int,
    drift_checker : float,
    chroma_checker : float,
    ) -> Axes :

    data = pd.DataFrame({
        "type" : ["Drift","Abberations"],
        "colocalization rate" : [drift_checker, chroma_checker]
    })
    data["colocalization rate"] *= 100

    sns.barplot(
        data,
        x="type",
        y="colocalization rate",
        hue="type", palette=["#fb9a99", "#80b1d3"],
        edgecolor = "black",
        ax=ax
    )

    ax.set_ylim(0,100)
    ax.set_yticks([])
    sns.despine(ax=ax, left=True)
    ax.set_xlabel("")
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30)

    drift_txt = ax.text(0, data.iloc[0].at["colocalization rate"] /2, f"{round(data.iloc[0].at["colocalization rate"],1)} %")
    abb_txt = ax.text(1, data.iloc[1].at["colocalization rate"] /2, f"{(round(data.iloc[1].at["colocalization rate"],1))} %")

    plt.setp([drift_txt,abb_txt], size=10, fontweight="bold", color="white", rotation = 15, va="center", ha="center")
    ax.set_ylabel(f"{coloc_range} nm")

    return ax

def plot_drift_means(
    ax : Axes,
    drift_means : DriftMeans,
    measure : Literal["z_mean_drift", "y_mean_drift","x_mean_drift","euclidian_mean_drift",]
    ) :

    data = pd.DataFrame(drift_means).reset_index(drop=False)

    sns.stripplot(
        data = data,
        x = measure,
        orient="h",
        edgecolor='black',
        linewidth=1,
        size=10,
        alpha=0.8,
        color='white',
        ax=ax
    )

    ax.set_xlabel("")
    ax.set_ylabel(measure.split("_")[0].capitalize(), size=15)
    ax.set_xlim(0)

    return ax

def plot_similarity_metrics(
    ax : Axes,
    similarity_df : pd.DataFrame,
    ) :

    sns.stripplot(
        data=similarity_df,
        x="type",
        y="mean_squared_error",
        hue="type",
        palette="bright",
        ax=ax,
        edgecolor="black",
        linewidth=1,
        alpha=0.8
    )

    sns.despine(ax=ax, bottom=False, left=True)

    return ax

def plot_drift_heatmap(
    ax : Axes,
    Drift
    ) :

    Drift = Drift.loc[Drift["cycle"] != 0]
    Drift.loc[:,["found"]] = (Drift["drift_z"] != 0) | (Drift["drift_y"] != 0) | (Drift["drift_x"] != 0)
    Drift = Drift.pivot(index="location",columns = "cycle", values="found",)

    sns.heatmap(
        Drift.transpose(), 
        cmap= ListedColormap(["#F89891","#B9FBBA"]),
        linecolor="white",
        linewidths=0.5,
        cbar=False,
        ax=ax
        )

    ax.set_xticks([0,len(Drift)])
    ax.set_xticklabels(["1",str(len(Drift))])
    ax.set_yticks([0,len(Drift.columns)])
    ax.set_yticklabels(["1",str(len(Drift.columns))])

    ax.set_xlabel("Location", size=15)
    ax.set_ylabel("Cycle", size=15)

    sns.despine(ax=ax, top=False, right=False)

    return ax

def _get_coloc_rates_for_checkers(
    run_path : str,
    **checkers : tuple[str,str]
    ) :

    table_path = os.path.join(run_path,"analysis","data","coloc_rates_mean.csv")
    if not os.path.isfile(table_path) : raise FileNotFoundError("Couldn't find coloc rates table, run first co-localization analysis")
    coloc_table = pd.read_csv(table_path, sep=";")
    coloc_table = coloc_table.set_index("target")

    coloc_rates = {key : float(coloc_table.at[pair[0],pair[1]]) for key,pair in checkers.items()}

    return coloc_rates

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


def compute_similarity_metrics(
    image1 : np.ndarray,
    image2 : np.ndarray
    ) :

    if np.issubdtype(image1.dtype, np.floating) :
        max_intensity = 1
    else :
        max_intensity = np.iinfo(image1.dtype).max

    ms_err = mean_squared_error(image1, image2)
    peak_snr = 10 * log10(max_intensity**2 / ms_err)

    res = {
        "spatial_correlation" : pearsonr(image1.ravel(), image2.ravel())[0],
        # "structural_similarity" : ssim(image1, image2, data_range=max_intensity),
        "mean_squared_error" : ms_err,
        "peak_snr" : peak_snr,
    }

    return SimilarityMetrics(**res)

def _process_location(
    Acquisition : pd.DataFrame,
    location : str,
    drift_initcycle : tuple[int,int],
    drift_endcycle : tuple[int,int],
    abb_initcycle : tuple[int,int],
    abb_endcycle : tuple[int,int],
    ) :

    location_stack = open_location(Acquisition, location)
    location_stack = location_stack.max(axis=1)
    drift_init = location_stack[drift_initcycle[0],...,drift_initcycle[1]]
    drift_end = location_stack[drift_endcycle[0],...,drift_endcycle[1]]
    abb_init = location_stack[abb_initcycle[0],...,abb_initcycle[1]]
    abb_end = location_stack[abb_endcycle[0],...,abb_endcycle[1]]
    del location_stack
    
    with ThreadPool() as pool :
        futures = pool.map(
            compute_similarity_metrics,
            [drift_init, abb_init],
            [drift_end,abb_end]
        )

    res = pd.DataFrame(futures.result())
    res["location"] = location
    res["type"] = ["drift", "abberation"]

    return res


def process_similarity_metrics(
    Gene_map : pd.DataFrame,
    Acquisition : pd.DataFrame,
    Drift : pd.DataFrame,
    drift_namepair : tuple[str,str],
    abb_namepair : tuple[str,str],
    run_path : str,
    ) :

    print(Gene_map)
    print(Gene_map['target'])
    Gene_map = Gene_map.set_index("target", verify_integrity=True)
    drift_initcycle = Gene_map.at[drift_namepair[0],"cycle"], Gene_map.at[drift_namepair[0],"color_id"]
    drift_endcycle = Gene_map.at[drift_namepair[1],"cycle"], Gene_map.at[drift_namepair[1],"color_id"]
    abb_initcycle = Gene_map.at[abb_namepair[0],"cycle"], Gene_map.at[abb_namepair[0],"color_id"]
    abb_endcycle = Gene_map.at[abb_namepair[1],"cycle"], Gene_map.at[abb_namepair[1],"color_id"]

    print("start pool")
    with ProcessPool(max_workers=4) as pool:
        locations = Acquisition["location"].unique()
        len_location = len(locations)
        futures = pool.map(
            _process_location,
            [Acquisition]*len_location,
            locations,
            [drift_initcycle]*len_location,
            [drift_endcycle]*len_location,
            [abb_initcycle]*len_location,
            [abb_endcycle]*len_location,
            )
        res = [r for r in list(tqdm(futures.result(), desc="Computing similarity measures", total=len_location))]
    res = pd.concat(res,axis=0)

    print("saving similarity metrics...")
    os.makedirs(os.path.join(run_path,"data"), exist_ok=True)
    res.to_csv(os.path.join(run_path,"data","similarity.csv"), sep=';')
    print("done")
    

    return res

def signal_quality_dashboard(
    run_path : str,
    Acquisition: pd.DataFrame,
    Cell: pd.DataFrame,
    Gene_map: pd.DataFrame,
    Detection: pd.DataFrame,
    Drift: pd.DataFrame,
    coloc_range : int,
    drift_checker : tuple[str,str],
    chroma_checker : tuple[str,str],
    ) :

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

    if any((not checker in Detection["target"].unique().tolist() for checker in drift_checker + chroma_checker)) :
        warnings.warn(f"Couldn't find one of {drift_checker +chroma_checker} in detection. Available columns are {Detection["target"].unique()}. \nEnsure correct spelling or that it was not filtered out for signal quality dashboard.")
        logging.error(f"Couldn't find one of {drift_checker +chroma_checker} in detection. Available columns are {Detection["target"].unique()}. \nEnsure correct spelling or that it was not filtered out for signal quality dashboard.")
        return None

    Cell = pd.merge(
        Cell,
        Detection.loc[:,["detection_id", "target","cycle","color_id","wavelength"]],
        on="detection_id",
    )

    #1.
    color_number = Detection["color_id"].nunique() +1

    #2.
    coloc_rates = _get_coloc_rates_for_checkers(
        run_path=run_path,
        drift_checker=drift_checker,
        chroma_checker=chroma_checker
    )

    #3. Similarity measurements
    similarity_path = os.path.join(run_path,"data","similarity.csv")
    if not os.path.isfile(similarity_path) :
        similarity_df = process_similarity_metrics(
            Gene_map=Gene_map,
            Acquisition=Acquisition,
            Drift=Drift,
            drift_namepair=drift_checker,
            abb_namepair=chroma_checker,
            run_path=run_path
            )
    else :
        similarity_df = pd.read_csv(similarity_path, sep=";")

    #4. Drift
    drift_means = _compute_drift_means(
        Acquisition=Acquisition,
        Drift=Drift
    )


    fig = plt.figure(figsize=(20,13))
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
    drift_heatmap_ax = fig.add_subplot(grid[3:, 2:])

    #Plot drift swarmplots
    max_drift = max([cast(pd.Series, axis).max() for axis in drift_means.values()])
    MEASURES : list[Literal["z_mean_drift","y_mean_drift","x_mean_drift","euclidian_mean_drift"]] = ["z_mean_drift","y_mean_drift","x_mean_drift","euclidian_mean_drift"]
    for i, measure in enumerate(MEASURES) :
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

    coloc_checkers_ax = plot_endcycle_coloc_rates(
        ax=coloc_checkers_ax,
        coloc_range = coloc_range,
        **coloc_rates
    )


    similarity_ax = plot_similarity_metrics(
        similarity_ax,
        similarity_df
    )

    drift_heatmap_ax = plot_drift_heatmap(
        drift_heatmap_ax,
        Drift=Drift
    )

    return fig