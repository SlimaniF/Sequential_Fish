import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as p
import logging, traceback

from matplotlib.figure import Figure

from ..tools import safe_merge_no_duplicates
from .utils import get_colors_list, merge_data

T_SCORE = 1.833

def cell_number(
    Cell : pd.DataFrame,
    frameon=False
) :
    
    data = Cell.loc[~Cell['target'].str.contains('Washout')]
    data = data.groupby(['target'], as_index=False).agg({
        'cell_id' : 'nunique'
        }).sort_values('target', ascending=False)

    colors = get_colors_list(len(data))

    fig = plt.figure(figsize=(15,4), frameon=frameon)
    ax = fig.gca()
    ax.bar(
        x= data['target'],
        height= data['cell_id'],
        color=colors,
        width=0.6
    )
    ax.set_ylabel('cell count')

    return fig

def _test_washout_parity(Spots : pd.DataFrame) :
    parity = Spots.loc[Spots['target'].str.contains('Washout')]['cycle'] % 2
    
    assert len(parity.unique()) == 1, f"Washout detected for pair and unpair cycles.\n{Spots.loc[Spots['target'].str.contains('Washout')].groupby('target')['cycle'].unique()}"
    
    washout_are_uneven_number = (parity == 1).all()
    
    return washout_are_uneven_number
    
def Spots_washout_filtering(
    Spots_with_washout : pd.DataFrame,
    frameon=False,
) :
    
    washout_data = Spots_with_washout.groupby(['cycle','target']).agg({
    'is_washout' : 'sum'
    })
    washout_are_uneven_numbers = _test_washout_parity(Spots_with_washout)
    
    
    fig = plt.figure(figsize=(20,8), frameon=frameon)
    ax = fig.gca()

    if washout_are_uneven_numbers : modulo = 0 #then look for even numbers

    data = washout_data[(washout_data.index.get_level_values(0) %2 == modulo)]

    x = np.arange(len(data))
    y = data['is_washout']
    ax.plot(x,y, color='red')
    ax.legend()

    cycles = data.index.get_level_values(0)
    targets = data.index.get_level_values(1)
    xlabels = [f"{target}\n({cycle})" for target, cycle in zip(targets,cycles)]
    ax.set_xticks(x, xlabels)

    fig.text(0.1, 0.5, 'Filtered spots count (#)', va='center', rotation='vertical')
    ax.set_xlabel('RNA (cycle)')

    return fig
    
def edge_and_segmentation_filtering(
    Unfiltered_spots : pd.DataFrame,
    Unfiltered_cells : pd.DataFrame,
    frameon = False,
) :
    
    spots_count_per_cell = Unfiltered_spots.groupby(['location','cell_label']).agg({
    'spot_id' : 'count',
    'target' : 'nunique',
    })
    
    remaining_cells = Unfiltered_cells.rename(columns={"label" : "cell_label"}).groupby(['location','cell_label'])['cell_id'].nunique()
    spots_count_per_cell['not_on_edge'] = spots_count_per_cell.index.isin(remaining_cells.index)
    spots_count_per_cell = spots_count_per_cell.rename(columns={'spot_id' : 'spot_count'})
    
    total_filtered = spots_count_per_cell[~spots_count_per_cell['not_on_edge']].agg({
    'spot_count' : ['sum', 'count']
    })
    total_filtered.index = pd.Index(['spot_count','cell_count'])
    total_filtered.columns = pd.Index(['filtered'])
    
    total_remaining = spots_count_per_cell[spots_count_per_cell['not_on_edge']].agg({
    'spot_count' : ['sum', 'count']
    })
    total_remaining.index = pd.Index(['spot_count','cell_count'])
    total_remaining.columns = pd.Index(['remaining'])
    data = pd.concat([total_filtered, total_remaining], axis=1)
    
    
    fig = plt.figure(figsize=(20,10),frameon=frameon)
    left,right = fig.subplots(1,2)

    labels = ['filtered','remaining']
    colors = ['red','green']

    x = np.arange(len(labels))

    #Cell
    left.bar(x,data.loc['cell_count',labels], color =colors, label= labels, alpha = 0.7, edgecolor= 'black')
    left.set_ylabel('cell count (#)')
    left.legend()
    left.set_xticks(x,labels)
    left.set_title("Cells")

    #Spots
    right.bar(x,data.loc['spot_count',labels], color =colors, label= labels, alpha = 0.7, edgecolor= 'black')
    right.set_ylabel('spot count (#)')
    right.legend()
    right.set_xticks(x,labels)
    right.set_title("Spots")

    fig.text(x=0.5, y=0.90, ha='center', s="Segmentation + Edge filtering", fontdict={'weight' : 'bold', 'size' : '14'})
    return fig

def cell_area(Cell : pd.DataFrame, frameon=False) :
    area_data = Cell.groupby('cell_id')['cell_area'].first()
    
    fig = plt.figure(figsize=(10,10), frameon=frameon)
    ax = fig.gca()
    
    count,bins,_ = ax.hist(
        area_data,
        bins=30,
        color='green',
        alpha = 0.5,
        edgecolor = 'black',
        )
    mean_area, median_area = area_data.mean(), area_data.median()
    xmin,xmax,ymin,ymax = ax.axis()
    ax.plot([mean_area,mean_area], [0, xmax], '--r', label= f'mean : {round(mean_area)} px')
    ax.plot([median_area,median_area], [0, xmax], '--', label=f'median : {round(median_area)} px')
    
    ax.axis([xmin,xmax,ymin,ymax])
    ax.legend()
    ax.set_xlabel('Cell area (px)')
    ax.set_ylabel('Count (#)')
    
    return fig

def _alignement_success(
    Acquisition : pd.DataFrame,
    Detection : pd.DataFrame,
    Drift : pd.DataFrame
    ) :
    
    Drift = safe_merge_no_duplicates(
        Drift,
        Acquisition,
        on= 'acquisition_id',
        keys='cycle'
    )
    
    dapi_mask = (~Drift['error'].isna()) & (Drift['drift_type'] == 'dapi')
    fish_mask = (~Drift['error'].isna()) & (Drift['drift_type'] == 'fish')
    dapi_sucess = Drift.loc[dapi_mask & ((Drift["drift_z"] != 0) | (Drift["drift_y"] != 0) | (Drift["drift_x"] != 0))]

    dapi_sucess_number = len(dapi_sucess)
    dapi_test_number = len(Drift[dapi_mask])
    dapi_fail_number = dapi_test_number - dapi_sucess_number

    dapi_sucess_percentage = dapi_sucess_number*100/dapi_test_number
    dapi_fail_percentage = dapi_fail_number*100/dapi_test_number
    fish_sucess = Drift.loc[fish_mask & ((Drift["drift_z"] != 0) | (Drift["drift_y"] != 0) | (Drift["drift_x"] != 0))]

    fish_sucess_number = len(fish_sucess)
    fish_test_number = len(Drift[fish_mask])
    fish_fail_number = fish_test_number - fish_sucess_number

    fish_sucess_percentage = fish_sucess_number*100/fish_test_number
    fish_fail_percentage = fish_fail_number*100/fish_test_number
    
    fig, axes = plt.subplots(1,2, figsize=(16,8), sharey=True)
    axes : 'list[plt.Axes]'
    fish, dapi = axes

    fish.bar(['sucess', 'failure'], [fish_sucess_percentage, fish_fail_percentage], color=['green', 'red'], alpha = 0.5, edgecolor= 'black')
    dapi.bar(['sucess', 'failure'], [dapi_sucess_percentage, dapi_fail_percentage], color=['green', 'red'], alpha = 0.5, edgecolor= 'black')

    fish.set_title('Fish')
    dapi.set_title('Dapi')
    fish.set_ylabel("Drift correction resulted in (%)")

    return fig

def _alignement_sucess_per_cycle(
    Drift : pd.DataFrame,
    Acquisition : pd.DataFrame,
    Detection : pd.DataFrame,
) :
    
    Drift = safe_merge_no_duplicates(
        Drift,
        Acquisition,
        keys= 'cycle',
        on= 'acquisition_id'
    )
    
    Drift['sucess'] = (Drift['drift_z'] != 0) | (Drift['drift_y'] != 0) | (Drift['drift_x'] | 0)
    per_cycle_sucess = Drift.loc[~Drift['error'].isna()].groupby(['drift_type','cycle']).agg({
        'sucess' : ['sum','count'],
    })

    per_cycle_sucess.columns = ['sucess', 'total']
    per_cycle_sucess['sucess_percentage'] = per_cycle_sucess['sucess']*100/per_cycle_sucess['total']
    
    fig = plt.figure(figsize=(8,8))
    ax = fig.gca()

    X = per_cycle_sucess.loc['fish'].index
    Y = per_cycle_sucess.loc['fish']['sucess_percentage']
    color = ['blue','orange']
    labels = ['fish', 'washout']

    legend_symb = [p.Rectangle((0,0),0,1, color= 'blue'), p.Rectangle((0,0),0,1, color= 'orange')]
    plt.legend(legend_symb, labels)

    ax.bar(X, Y, color=color, alpha = 0.8, edgecolor='black')
    ax.set_ylabel('Sucess rate')
    ax.set_xlabel('cycle')
    ax.set_title('Alignement sucess rate per cycle')

    return fig

def _compute_euclidian_drift(df: pd.DataFrame, voxel_size) :
    df.loc[:,['drift_z', 'drift_y', 'drift_x']] *= voxel_size
    df_sq = df.apply(np.square).sum(axis=1)
    df['euclidian_drift'] = df_sq.apply(np.sqrt)
    return df

def _drift_distance(
    Drift : pd.DataFrame,
    Acquisition : pd.DataFrame,
    voxel_size
) :
    
    Drift = safe_merge_no_duplicates(
        Drift,
        Acquisition,
        keys= 'cycle',
        on= 'acquisition_id'
    )
    
    fish_drift = Drift[(~Drift['error'].isna()) & (Drift['drift_type'] == 'fish')].loc[:,['cycle','drift_z', 'drift_y', 'drift_x']]
    fish_drift = _compute_euclidian_drift(fish_drift, voxel_size)
    dapi_drift = Drift.loc[(~Drift['error'].isna()) & (Drift['drift_type'] == 'dapi')].loc[:,['cycle','drift_z', 'drift_y', 'drift_x']]
    dapi_drift = _compute_euclidian_drift(dapi_drift, voxel_size)
    
    fig, axes = plt.subplots(1,2, figsize=(16,8), sharey=True)
    fish, dapi = axes

    fish_data = fish_drift.drop(columns='cycle')
    dapi_data = dapi_drift.drop(columns='cycle')

    fish_mean = fish_data.apply(np.abs).mean(axis=0)

    fish_std = fish_data.apply(np.abs).std(axis=0)
    fish_count = fish_data.count(axis=0)
    confidence_interval = (T_SCORE*fish_std) / fish_count.apply(np.sqrt)

    colors = ['gray', 'gray', 'gray', 'red']

    fish.bar(
        fish_mean.index, 
        fish_mean, 
        yerr=confidence_interval, 
        capsize=5,
        color=colors,
        alpha=0.8,
        edgecolor='black'
        )

    xmin,xmax,ymin,ymax = fish.axis()
    fish.plot([xmin,xmax], [0,0],'k')
    fish.set_title('Fish')
    fish.set_ylabel('Absolute drift distance average (nm)')

    dapi_mean = dapi_data.apply(np.abs).mean(axis=0)
    dapi_std = dapi_data.apply(np.abs).std(axis=0)
    dapi_count = dapi_data.count(axis=0)
    confidence_interval = (T_SCORE*dapi_std) / dapi_count.apply(np.sqrt)

    dapi.bar(
        dapi_mean.index, 
        dapi_mean, 
        yerr=confidence_interval, 
        capsize=5,
        color=colors,
        alpha=0.8,
        edgecolor='black'
        )

    xmin,xmax,ymin,ymax = dapi.axis()
    dapi.set_title('Dapi')

    fig.suptitle(f"Voxel size (zyx) : {voxel_size}")

    return fig

def _drift_distance_per_cycle(
    Drift : pd.DataFrame,
    Acquisition : pd.DataFrame,
    voxel_size
) :
    
    Drift = safe_merge_no_duplicates(
        Drift,
        Acquisition,
        keys= 'cycle',
        on= 'acquisition_id'
    )
    
    fish_drift = Drift[(~Drift['error'].isna()) & (Drift['drift_type'] == 'fish')].loc[:,['cycle','drift_z', 'drift_y', 'drift_x']]
    fish_drift = _compute_euclidian_drift(fish_drift, voxel_size)
    fish_data = fish_drift.groupby('cycle').agg({
    'euclidian_drift' : ['mean', 'std', 'count']
    })

    fish_data[('euclidian_drift','confidence_interval')] = (T_SCORE * fish_data.loc[:,("euclidian_drift",'std')]) / fish_data.loc[:,("euclidian_drift",'count')].apply(np.sqrt)
    
    fig = plt.figure(figsize=(8,8))
    ax = fig.gca()

    X = fish_data.index
    Y = fish_data.loc[:,('euclidian_drift','mean')]
    error_bar = fish_data.loc[:,('euclidian_drift','confidence_interval')]
    color = ['green','lightblue']
    labels = ['fish', 'washout']

    legend_symb = [p.Rectangle((0,0),0,1, color= color[0]), p.Rectangle((0,0),0,1, color= color[1])]
    plt.legend(legend_symb,labels)

    ax.bar(X, Y, color=color, yerr = error_bar, capsize=5, alpha= 0.8, edgecolor='black')
    ax.set_ylabel('Euclidian drift distance (nm)')
    ax.set_xlabel('cycle')
    ax.set_title('Drift distance per cycle')

    xmin,xmax,ymin,ymax = plt.axis()

    ax.set_ylim(0,ymax)

    return fig

def drift_statistics(
    Acquisition : pd.DataFrame,
    Detection : pd.DataFrame,
    Drift : pd.DataFrame
    ) :
    #Get voxel_size
    voxel_size = tuple(Detection['voxel_size'].iat[0])
    voxel_size = [int(i) for i in voxel_size]
    
    sucess_rate_fig = _alignement_success(Acquisition, Detection, Drift)
    sucess_rate_per_cycle_fig = _alignement_sucess_per_cycle(Drift, Acquisition, Detection)
    drift_distance_fig = _drift_distance(Drift=Drift, Acquisition=Acquisition, voxel_size=voxel_size)
    distance_per_cycle_fig = _drift_distance_per_cycle(Drift, Acquisition=Acquisition, voxel_size=voxel_size)
    
    figures = {
        "sucess_rate" : sucess_rate_fig,
        "sucess_rate_per_cycle" : sucess_rate_per_cycle_fig,
        "drift_distance" : drift_distance_fig,
        "distance_per_cycle" : distance_per_cycle_fig,
    }
    
    return figures

def pipeline_metrics(
    Acquisition : pd.DataFrame,
    Detection : pd.DataFrame,
    Gene_map : pd.DataFrame,
    Spots_with_washout : pd.DataFrame,
    Unfiltered_spots : pd.DataFrame,
    Cell : pd.DataFrame,
    Drift : pd.DataFrame,
    run_path : str,
    frameon = True,
) :
    
    Detection, Cell, Unfiltered_spots = merge_data(
            Acquisition=Acquisition,
            Detection=Detection,
            Cell=Cell,
            Spots=Unfiltered_spots,
            Gene_map=Gene_map
        )
    
    Spots_with_washout = safe_merge_no_duplicates(
        Spots_with_washout,
        Unfiltered_spots,
        on='spot_id',
        keys='target'
    )
    
    output_path = run_path + "/analysis/pipeline_metrics/"
    os.makedirs(output_path, exist_ok=True)
    
    log_file = output_path + "/pipeline_metrics_log.log"
    logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    force= True
    )
    error_count = 0
    
    analysis_dict = {
        'cell_number' : cell_number, 
        'washout_filtering' : Spots_washout_filtering, 
        'edges_figures_filtering' : edge_and_segmentation_filtering,
        'cell_area' : cell_area,
    }
    smfishtools = {
        'cell_number' : {'Cell' : Cell, 'frameon' : frameon},
        'washout_filtering' : {'Spots_with_washout' : Spots_with_washout, 'frameon' : frameon},
        'edges_figures_filtering' : {'Unfiltered_cells' : Cell, 'Unfiltered_spots' : Unfiltered_spots, 'frameon' : frameon},
        'cell_area' : {'Cell' : Cell, 'frameon' : frameon},
    }
    
    logging.info(f"Start general pipeline metrics analysis")
    for analysis_name in analysis_dict.keys() :
        logging.info(f"Start {analysis_name} analysis")
        try :
            analysis = analysis_dict[analysis_name]
            parameters = smfishtools[analysis_name]
            fig : Figure = analysis(**parameters)
            fig.savefig(f"{output_path}/{analysis_name}.svg")
            plt.close()
            
        except Exception as e :
            logging.error(f"analysis {analysis_name} failed :\n{e}")
            error_count += 1
            continue

        else : 
            logging.info(f"analysis {analysis_name} succeeded")
            continue

    try :
        logging.info(f"Start Drift analysis")
        drift_figures = drift_statistics(
            Acquisition=Acquisition,
            Detection=Detection,
            Drift=Drift
        )

        os.makedirs(f"{output_path}/drift_statistics/",exist_ok=True)
        for key in drift_figures :
            fig : Figure = drift_figures[key]
            fig.savefig(f"{output_path}/drift_statistics/{key}.svg")

    except Exception as e :
        logging.error(f"Drift analysis failed :\n{traceback.format_exc()}")
        error_count +=1
        return False

    else : 
        logging.info(f"Drift analysis succeeded")
        return True

    finally :
        if error_count > 0 :
            print(f"Pipeline metrics finished with {error_count} error(s).")
        plt.close()