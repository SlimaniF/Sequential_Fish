from concurrent import futures
import os
from typing import cast

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.gridspec as gridspec
from matplotlib.axes import Axes
import seaborn as sns
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram

import time as t
from pebble import ThreadPool



def plot_cell_number_per_location(
    ax: Axes, 
    cell_distribution : pd.DataFrame
    ) -> Axes:

    g = sns.barplot(
        data=cell_distribution, 
        x="cell number", y="location", 
        hue="location",palette="bright",
        legend=False,
        alpha=.6, 
        ax=ax,
        )
    sns.despine(left=True)
    xlabel = ax.set_xlabel("Cell number")
    ylabel = ax.set_ylabel("Location")

    yticks = ax.get_yticks()
    yticklabels = ax.get_yticklabels()
    # Only show every 5th tick label
    ax.set_yticks(yticks[::5])
    ax.set_yticklabels(yticklabels[::5])

    #Display mean value
    mean = float(cell_distribution["cell number"].median())
    mean_line = ax.plot(
        [mean,mean],
        [-1,len(cell_distribution)],
        "--",
        label='median',
        color="#231A42"
    )
    ax.set_ylim(len(cell_distribution),-0.5,)
    ax.legend()

    plt.setp([xlabel,ylabel],  size=15)

    return ax




def plot_analyzed_cell_number(
    ax : Axes,
    total_cell_number : int,
    analyzed_cell_number : int
) -> Axes :

    discarded_cell_fraction = (total_cell_number - analyzed_cell_number) / total_cell_number
    analyzed_cell_fraction = analyzed_cell_number / total_cell_number
    
    outer_pie, labels, *_ = ax.pie(
        [discarded_cell_fraction, analyzed_cell_fraction],
        explode=[0.1,0.1],
        radius= 1.3, 
        colors=["red","green"],
        labels = [f"{total_cell_number-analyzed_cell_number}\nFiltered",f"{analyzed_cell_number}\nAnalyzed"],
        labeldistance=0.8,
        )
    central_text = ax.text(0,0,f"{total_cell_number}\ndetected\ncells")

    plt.setp(outer_pie, width=0.5, edgecolor='black',)
    plt.setp(labels, weight="bold", horizontalalignment="center",va="center")
    plt.setp(central_text, size= 30, weight="bold", horizontalalignment="center",va="center")

    return ax

def _create_cell_distribution(Cell : pd.DataFrame) -> pd.DataFrame:
    cell_distribution = pd.DataFrame(Cell.groupby("location")["cell_id"].nunique()).rename(columns={"cell_id" : "cell number"}).reset_index(drop=False)
    cell_distribution["location"] = cell_distribution["location"].str.replace("Location-","")

    return cell_distribution

def _count_total_cell_number(
    run_path : str,
) -> int :

    segmentation_folder = os.path.join(run_path, "segmentation")

    def _count_cell_in_location(path : str) :
        arr = np.load(path)["cytoplasm"]
        return len(np.unique(arr))-1

    path_list = [os.path.join(segmentation_folder,file) for file in os.listdir(segmentation_folder)]
    
    with ThreadPool() as workpool :
        futures = workpool.map(
            _count_cell_in_location,
            path_list
        )

    cell_number = sum(futures.result())


    return cell_number

def _count_analyzed_cell_number(Cell : pd.DataFrame) -> int :
    analyzed_cell_number = Cell.groupby("location")["cell_id"].nunique().sum()
    return analyzed_cell_number
    

def plot_spot_is_found_heatmap(
    ax : Axes,
    spot_is_found_per_cell : pd.DataFrame
    ) ->Axes :

    linkage_matrix = linkage(spot_is_found_per_cell)
    try :
        dn = dendrogram(linkage_matrix, no_plot=True)
    except RecursionError as e :
        if str(e) == "maximum recursion depth exceeded" :
            target_order = np.arange(len(spot_is_found_per_cell))
        else :
            raise e
    else :
        target_order = dn["leaves"]
        del dn
    
    spot_is_found_per_cell = spot_is_found_per_cell.iloc[target_order]
 
    sns.heatmap(
        spot_is_found_per_cell.transpose(), 
        cmap= ListedColormap(["#F89891","#B9FBBA"]),
        cbar=False,
        ax=ax
        )
    
    sns.despine(ax=ax, bottom=False, left=False,top=False,right=False)
    ax.set_ylabel("")

    ax.set_yticks(ax.get_yticks(), cast(list[str],ax.get_yticklabels()), rotation = -15, va="bottom")
    ax.set_xlabel("Population was found in cell",size=15)
    ax.set_xticks([])
    ax.get_xaxis().set_label_position("top")


    return ax




def _create_spot_per_cell_distribution(Spots : pd.DataFrame) :
    spot_per_cell_distribution : pd.DataFrame =  Spots.loc[Spots["cell_id"] != 0].groupby(["location","cell_id","target"])["spot_id"].nunique().rename('spot number')
    spot_per_cell_distribution = spot_per_cell_distribution.reset_index(drop=False).pivot(columns="target",index=["location","cell_id"],values="spot number")

    return spot_per_cell_distribution


def plot_spot_per_cell_distribution(
    ax : Axes,
    spot_per_cell_distribution : pd.DataFrame,
    ) -> Axes :

    data = spot_per_cell_distribution.melt()
    data["target"] = data["target"].str.replace("_","\n")

    sns.boxplot(
        data,
        x="target",
        y="value",
        # hue="target", palette="bright",
        color="white",
        linewidth=2,
        linecolor="black",
        showfliers=False,
        ax=ax
        )

    xlabel = ax.set_xlabel("")
    ylabel = ax.set_ylabel("Single molecule per cell (#)")
    plt.setp([xlabel,ylabel],  size=15)


    return ax

def cell_dashboard(
    run_path : str,
    Cell : pd.DataFrame,
    Spots : pd.DataFrame
    ) :
    #1.
    total_cell_number = (_count_total_cell_number(run_path))
    analyzed_cell_number = (_count_analyzed_cell_number(Cell))

    #2.
    cell_distribution = _create_cell_distribution(Cell)

    #3.
    spot_per_cell_distribution = _create_spot_per_cell_distribution(Spots)
    spot_is_found_per_cell = spot_per_cell_distribution > 0

    #4.


    fig = plt.figure(figsize=(40,20))
    grid = gridspec.GridSpec(2,2,width_ratios=[1,3],height_ratios=[1,2])
    topleft = fig.add_subplot(grid[0,0])
    topright = fig.add_subplot(grid[0,1])
    bottomleft = fig.add_subplot(grid[1,0])
    bottomright = fig.add_subplot(grid[1,1])

    topleft = plot_analyzed_cell_number(
        topleft,
        total_cell_number,
        analyzed_cell_number
    )

    bottomleft = plot_cell_number_per_location(
        ax=bottomleft,
        cell_distribution=cell_distribution
    )

    topright = plot_spot_is_found_heatmap(
        topright,
        spot_is_found_per_cell
    )

    bottomright = plot_spot_per_cell_distribution(
        ax=bottomright,
        spot_per_cell_distribution=spot_per_cell_distribution
    )

    return fig

