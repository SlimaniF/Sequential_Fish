"""
Submodule for co-localization analysis

# Notation

We try co-localization between 2 populations *population1* --> *population2*, meaning we count how many single from population 1 
co-localize with population 2.
Those populations can be selected from 'all', 'free' or 'clustered' spots.  
    
V is volume in pixel by which we mean the number of available position a single molecule can take in a cell
    
# Notes

- **Volume** :  
    In analysis we can't acess volume of a cell with precision as cell are segmented in 2D, computing volume by multiplying 
Area with z-stack number would consider a lot of pixel that single molecules can't acess in reality. To work around this we compute first
the average number of single molecule per plane and work with modelisation in 2D as our model is dimensionly independent. For each cell
number of single molecule per plane and 'volume' of plane (pixels) then we compute expected colocalization using  
**'plane abudancy'** = volume x spot per plane.

- **Statistical test** :  
    With our modelisation we have acces to probabilistic densities of colocalization events for each cell. Though for each cell 
    (i.e particular set of volume and abudancies) we only have one measure of co-colocalization event meaning we can't statisticaly test
    the relevance of this measurement, anyway it would not be interesting as **the normalised co-localisation rate** already
    is a indication of how far we are from random distribution. When then average these events measurement on the whole cell population,
    the statistical test is to know if this average value (also average normalised valued) is statistically significant compared to what
    we would expect from a set of cells with randomly distributed single molecules. To do so we need to compute, in a way,
    the distribution of distributions of co-localization events for a bach of cells (i.e a set of volume and abundancies).
    
    Hopefully, **the distribution of random variables following normal law also follow normal law** and we can compute its mean and std.
    (from propertis of gaussan and linear combinations)
    Mean is the mean of each individual distribution (knowing volume and abudancies) (law of total expectation)
    Std is the sum of squared variance(or std ???) divided by N squared.
    
    From this new normal distribution we can use usual statistical tests.

"""
import matplotlib.colorbar
import pandas as pd
import numpy as np

from typing import Literal

from ..tools import safe_merge_no_duplicates

from sklearn.neighbors import NearestNeighbors
from scipy.ndimage import distance_transform_edt

"""
1. Data co-localization pariwise rates : Compute actual co-colocalization pairwise rates.
"""

def _get_population_index(
    Spots : pd.DataFrame, 
    population_key:Literal['all','free','clustered']
    ) :
    """
    Get index from Spots data frame to a population : 'all', 'clustered', 'free'
    """
    if population_key == 'all' : population_index = Spots.index
    elif population_key == 'clustered' : population_index = Spots.loc[~Spots['cluster_id'].isna()].index
    elif population_key == 'cluster' : population_index = Spots.loc[~Spots['cluster_id'].isna()].index
    elif population_key == 'clusters' : population_index = Spots.loc[~Spots['cluster_id'].isna()].index
    elif population_key == 'free' : population_index = Spots.loc[Spots['cluster_id'].isna()].index
    else : raise AssertionError("{} incorect key for population_key".format(population_key))

    return population_index

def _create_coordinate_df(
    Spots : pd.DataFrame, 
    population_key:Literal['all','free','clustered']
    ) :
    """
    Prepare dataframe with nanometer coordinates
    """
    
    population_index = _get_population_index(Spots, population_key)
    coordinates_df = Spots.loc[population_index]

    #converting pixel coordinates to nanometers
    coordinates_df['voxel_size_z'], coordinates_df['voxel_size_y'], coordinates_df['voxel_size_x'] = list(zip(*coordinates_df['voxel_size']))
    coordinates_df['z'] *= coordinates_df['voxel_size_z'] 
    coordinates_df['y'] *= coordinates_df['voxel_size_y'] 
    coordinates_df['x'] *= coordinates_df['voxel_size_x'] 

    coordinates_df['coordinates'] = list(zip(coordinates_df['z'], coordinates_df['y'], coordinates_df['x']))
    coordinates_df = coordinates_df.groupby(['location','target'])['coordinates'].apply(list)
    return coordinates_df

def _create_neighbor_model_dict(
    spots_coordinates_df : pd.DataFrame, 
    colocalisation_distance : int
    ) :
    """
    Prepare for each single distribution a Nearestneighbor model, population passed to those must be 'population2' in other words the population where co-localization is tried WITH. See submodule description above.
    """
    neighbor_models_dict = dict()
    for idx in spots_coordinates_df.index :
        spot_distribution = spots_coordinates_df.at[idx]
        new_model = NearestNeighbors(radius=colocalisation_distance)
        new_model.fit(spot_distribution)
        neighbor_models_dict[idx] = new_model
        
    return neighbor_models_dict

def _compute_colocalisation_truth_df(
    spots_coordinates_df : pd.DataFrame, 
    Spots : pd.DataFrame, 
    neighbor_models_dict : dict,
    population_1 : Literal['all','clustered','free'],
    self_colocalisation = False,
    ) :
    
    population1_index = _get_population_index(Spots, population_key=population_1)
    RNAs = list(spots_coordinates_df.index.get_level_values(1).unique())
    colocalisation_truth_df = pd.DataFrame(index=population1_index, columns= RNAs, dtype=bool)
    colocalisation_truth_df = colocalisation_truth_df.join(Spots.loc[:,['spot_id','location','target', 'z','y','x','voxel_size']])

    #converting coordinates to nanometers
    colocalisation_truth_df['voxel_size_z'], colocalisation_truth_df['voxel_size_y'], colocalisation_truth_df['voxel_size_x'] = list(zip(*colocalisation_truth_df['voxel_size']))
    colocalisation_truth_df['z'] *= colocalisation_truth_df['voxel_size_z'] 
    colocalisation_truth_df['y'] *= colocalisation_truth_df['voxel_size_y'] 
    colocalisation_truth_df['x'] *= colocalisation_truth_df['voxel_size_x'] 
    colocalisation_truth_df['coordinates'] = list(zip(colocalisation_truth_df['z'], colocalisation_truth_df['y'], colocalisation_truth_df['x']))

    colocalisation_truth_df = colocalisation_truth_df.drop(columns=['z','y','x','voxel_size','voxel_size_z','voxel_size_y','voxel_size_x'])
    
    for location in colocalisation_truth_df['location'].unique() :
        target_idx = colocalisation_truth_df[colocalisation_truth_df['location'] == location].index
        for rna in RNAs :
            model : NearestNeighbors = neighbor_models_dict[(location, rna)]
            coordinates = list(colocalisation_truth_df.loc[target_idx]['coordinates'].apply(np.array,dtype=int))
            coordinates = np.array(coordinates, dtype=int)

            query = model.radius_neighbors(coordinates, return_distance=False)
            if self_colocalisation :
                query = pd.Series(query, index=target_idx).apply(len) > 1
            else :
                query = pd.Series(query, index=target_idx).apply(len).apply(bool) #if count is 0 no colocalisation -> False else True

            colocalisation_truth_df.loc[target_idx,[rna]] = query
    
    return colocalisation_truth_df

def _spots_merge_data(
    Spots : pd.DataFrame,
    Acquisition : pd.DataFrame,
    Detection : pd.DataFrame,
    Gene_map : pd.DataFrame,
    Cell : pd.DataFrame,
    ) :
    """
    Merge required information into Spots df
    """
    
    Detection = safe_merge_no_duplicates(
    Detection,
    Acquisition,
    on= ['acquisition_id'],
    keys=['cycle','location', 'fish_reodered_shape']
)

    Detection = safe_merge_no_duplicates(
        Detection,
        Gene_map,
        on= ['cycle','color_id'],
        keys=['target']
    )

    Spots =safe_merge_no_duplicates(
        Spots,
        Detection,
        on= 'detection_id',
        keys= ['location','target', 'voxel_size', 'fish_reodered_shape']
    )

    Spots = safe_merge_no_duplicates(
        Spots,
        Cell.rename(columns={'label' : 'cell_label'}),
        on=['acquisition_id','detection_id','cell_label'],
        keys=['cell_id']
    )


    return Spots

def colocalisation_truth_df(
    Spots : pd.DataFrame,
    Acquisition : pd.DataFrame,
    Detection : pd.DataFrame,
    Gene_map : pd.DataFrame,
    Cell : pd.DataFrame,
    population_1 : Literal['all','clustered','free'] = 'all',
    population_2 : Literal['all','clustered','free']= 'all',
    colocalisation_distance : int = 400,
    ) :
    """
    
    Create a dataframe where each line corresponds to one spot

    PARAMETERS
    ----------
        Spots must contain voxel_size.
    
    KEYS
    ----
        - `location` : field of fiew reference
        - `target` : from which distribution comes this spot
        - `coordinates` : ...
        - `spot_id` : unique identifier to spot
        
        - `boolean key` : + one key for each different `target` value, representing a boolean value indicating TRUE if this spot co-localize with this distribution
    
    """
    Spots = _spots_merge_data(
        Spots=Spots,
        Detection=Detection,
        Acquisition=Acquisition,
        Gene_map=Gene_map,
        Cell=Cell
    )
    
    population_1_index = _get_population_index(
        Spots, 
        population_key=population_1
        )
    
    real_spots_coordinates_df = _create_coordinate_df(
        Spots, 
        population_key= population_2
        )
    
    neighbor_models_dict = _create_neighbor_model_dict(
        real_spots_coordinates_df, 
        colocalisation_distance=colocalisation_distance
        )
    
    colocalisation_truth_df = _compute_colocalisation_truth_df(
        real_spots_coordinates_df, 
        Spots.loc[population_1_index], 
        neighbor_models_dict,
        population_1=population_1,
        self_colocalisation=False,
        ).set_index('spot_id', drop=False, verify_integrity=True)
    
    self_colocalization = _compute_colocalisation_truth_df(
        real_spots_coordinates_df, 
        Spots.loc[population_1_index], 
        neighbor_models_dict,
        population_1=population_1,
        self_colocalisation=True,
        ).set_index('spot_id', drop=False, verify_integrity=True)
    
    for rna in list(Spots['target'].unique()) :
        target_spots_index = colocalisation_truth_df.loc[colocalisation_truth_df['target'] == rna].index
        colocalisation_truth_df.loc[target_spots_index, [rna]] = self_colocalization.loc[target_spots_index, [rna]]
    
    return colocalisation_truth_df.reset_index(drop=True)

def create_cell_coloc_rates_df(
    Spots : pd.DataFrame, 
    colocalisation_truth_df : pd.DataFrame,
    ) -> pd.DataFrame:
    """
    Create dataframe where each index and columns correspond to `(target, cell_id)` values (i.e one line per distributions, cell_id couple)
    
    Colocalization rate of target i with target j is found at line i, column j (i.e `df.at[i,j]`)
    
    Result is ready for mean calculation or normalisation.
    
    """
    RNA_list = list(Spots['target'].unique())
    RNA_list.sort()
    
    colocalisation_truth_df= safe_merge_no_duplicates(
        colocalisation_truth_df,
        Spots,
        on='spot_id',
        keys='cell_id'
    )

    cell_coloc_rates = colocalisation_truth_df.groupby(['target','cell_id'])[RNA_list].mean() #Normalisation needs to happen after this


    return cell_coloc_rates

def compute_coloc_rates_mean(
    cell_coloc_rates : pd.DataFrame,
    ) :
    
    """
    Return dataframe with one line per distribution with mean value (amongst cell population) of co-localization.
    Colocalization rate of target i with target j is found at line i, column j (i.e `df.at[i,j]`)
    """
    
    coloc_rates = cell_coloc_rates.groupby('target', axis=0, level=0).mean()
    return coloc_rates
    

"""
2. Co-localization normalisation : Computes co-localization scores from modelisation

Aim : normalize cell_coloc_rates from part 1. with modelisation expectancy.
"""

def _compute_corrected_positions_number(
    voxel_size : 'tuple[int]',
    colocalisation_distance : int,
    ) :

    dim = len(voxel_size)

    c = int(
        2*(colocalisation_distance/min(voxel_size))
        )
        
    if c % 2 == 0 : c+=1 #We want c uneven so it is easier to find the center
    center_index = int(
        np.floor(c/2)
    )
    scanned_volume = np.ones((c,)*dim)

    if dim == 3 :
        scanned_volume[center_index,center_index,center_index] = 0
    elif dim == 2 :
        scanned_volume[center_index,center_index,center_index] = 0
    else : raise ValueError(f"Unsupported dimension for voxel size. Must be 2 or 3, it is {dim}")
        
    scanned_volume = distance_transform_edt(scanned_volume, sampling=voxel_size)
    scanned_volume = scanned_volume <= colocalisation_distance
    corrected_positions_number = scanned_volume.sum()

    return corrected_positions_number

def _get_cell_area(
    Cell : pd.DataFrame,
    ) :
    assert (Cell.groupby('cell_id',as_index=True)['cell_area'].unique().apply(len) == 1).all(), "Cell area is not unique for cell_id"
    Cell_area = Cell.groupby('cell_id',as_index=True)['cell_area'].first() #is unique so we take first

    return Cell_area

def _get_spot_per_plane(
    Spots : pd.DataFrame,
    ) :
    Cell_spots_count : pd.DataFrame = Spots.groupby(['cell_id','target','z'], as_index=False)['spot_id'].count()
    Cell_spots_count : pd.DataFrame = Cell_spots_count.groupby(['cell_id','target'], as_index=False)['spot_id'].mean().rename(columns={'spot_id' : 'spot_per_plane'})
    Cell_spots_count = Cell_spots_count.pivot(columns='target',index='cell_id',values='spot_per_plane')

    return Cell_spots_count

def _compute_spot_density(
    spots_per_plane : pd.DataFrame,
    Cell_area : pd.Series,
    RNA_list : 'list[str]',
    ):
    """
    For each distribution (columns) compute density of single molecule per pixel, one line per cell_id which is the co-localization probability.
    Indeed we divide number of single in volume by volume.
    """
    res = spots_per_plane.copy()

    for rna in RNA_list :
        res[rna] = spots_per_plane[rna]/Cell_area
    
    return res

def _compute_selfcoloc_rates(
        Cell_area : pd.DataFrame,
        Cell_abundancies : pd.DataFrame,
        RNA_list : list,
) :
    
    """
    pself colocalization = 1 - (V/k).(1 - (1 - 1/V )**k)
    """

    self_coloc_rates = pd.DataFrame(columns=Cell_abundancies.columns, index= Cell_abundancies.index, dtype=float)
    for rna in RNA_list :
        self_coloc_rates.loc[:,[rna]] = (1-pd.DataFrame((1-1/Cell_area).rename(rna)).pow(Cell_abundancies)).divide(Cell_abundancies)
        self_coloc_rates.loc[:,[rna]] = 1-self_coloc_rates.loc[:,[rna]].multiply(pd.DataFrame(Cell_area.rename(rna)))

    return self_coloc_rates

def create_coloc_rate_expectancy(
    Spots : pd.DataFrame,
    Cell : pd.DataFrame,
    voxel_size : 'tuple[int]',
    colocalisation_distance : int,
    RNA_list : list
) -> 'tuple[pd.DataFrame,pd.DataFrame]' :
    """
    Return dataframes with one line per cell with colocalisation rate of a single molecule with distributions (1 distribution per column).
    Returns (coloc_rates_df, selfcoloc_rates_df)

    """
    
    Spots = Spots.loc[Spots['target'].isin(RNA_list)]

    Cell_area = pd.DataFrame(_get_cell_area(Cell).rename('real_area'))
    corrected_positions_number = _compute_corrected_positions_number(
        voxel_size=voxel_size,
        colocalisation_distance=colocalisation_distance
    )
    Cell_area['corrected_area'] = Cell_area['real_area'] / corrected_positions_number
    
    
    Spots_per_plane = _get_spot_per_plane(Spots)

    coloc_probabilty = _compute_spot_density( # is alreay the co-localization probability
        spots_per_plane=Spots_per_plane,
        Cell_area=Cell_area['corrected_area'],
        RNA_list= RNA_list,
    )

    selfcoloc_probabilty = _compute_selfcoloc_rates(
        Cell_area=Cell_area['corrected_area'],
        Cell_abundancies= Spots_per_plane,
        RNA_list=RNA_list,
    )
    
    return coloc_probabilty, selfcoloc_probabilty

def _compute_cell_distribution_populations(
        Spots : pd.DataFrame,
) :
    Cell_spots_count : pd.DataFrame = Spots.groupby(['cell_id','target'])['spot_id'].count().rename('abundancy')

    Cell_spots_count = Cell_spots_count.reset_index(drop=False).pivot(
        index='cell_id',
        columns='target',
        values='abundancy'
    )

    return Cell_spots_count

def compute_z_score_frame(
        measured_colocalisation_events : pd.DataFrame,
        expected_colocalisation_events : pd.DataFrame,
        expected_standard_deviation : pd.DataFrame, 
        ) :
    if not measured_colocalisation_events.index.equals(expected_colocalisation_events.index) or not measured_colocalisation_events.index.equals(expected_standard_deviation.index) :
        measured_colocalisation_events = measured_colocalisation_events.reindex_like(expected_colocalisation_events)
        
    if not measured_colocalisation_events.columns.equals(expected_colocalisation_events.columns) or not measured_colocalisation_events.columns.equals(expected_standard_deviation.columns) :
        raise ValueError("Cannot compute z score, all dataframe don't share the same columns")
    
    expected_standard_deviation = expected_standard_deviation.replace(0, np.NAN) # Avoid division by 0 error

    z_score_df = measured_colocalisation_events - expected_colocalisation_events
    z_score_df = z_score_df/expected_standard_deviation

    return z_score_df
    
"""
3. Statistical test (p-value) VS Null-model

We already created function for mean computation in a distribution -> `_compute_coloc_rate_expectancy`. But we want the equivalent with standard deviation.

"""

from scipy.stats import shapiro #normality test
from scipy.stats import wilcoxon

def compute_wilcoxon_signed_rank(zscore_distribution)->pd.Series :

    if isinstance(zscore_distribution, (pd.Series, pd.DataFrame)) :
        name = zscore_distribution.name
    else :
        name =  None

    statistic, pvalue = wilcoxon(zscore_distribution)

    if isinstance(pvalue, (int, float)) : pvalue = [pvalue]
    if isinstance(zscore_distribution, pd.DataFrame) :
        res = pd.DataFrame(columns=zscore_distribution.columns, data=pvalue.reshape(1,-1), index=[name])
    else :
        res = pd.Series(pvalue, name=name)
    return res

def compute_pvalue_frame(
        zscore_frame : pd.DataFrame,
) : 
    pvalue_frame = zscore_frame.groupby(axis=0,level=0).apply(lambda x: compute_wilcoxon_signed_rank(x)).droplevel(axis=0, level=1)

    return pvalue_frame    

"""
4. Higher dimension co-localization tests
"""

def get_significative_combination() :
    pass

def get_next_combinations() :
    pass

def get_combinations_abundancies() :
    pass

"""
5. Plots
"""

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, AsinhNorm, LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.cm import get_cmap
import matplotlib

def score_color_scale(vmax=2) :
    
    if vmax >= 2 :
        end_gray = 0.5 + (1/vmax)
    elif vmax <= 0 :
        raise ValueError("vmax cannot be <= 0")
    else :
        end_gray = 0.6
    
    orange_pos  = end_gray + (1 -end_gray) / 2

    colors = [
        (0,'#7d80fc'),
        (0.5,'gray'),
        (end_gray,'yellow'),
        (orange_pos, 'orange'),
        (1,'red'),
        ]
    
    
    colormap = LinearSegmentedColormap.from_list("custom_cmap", colors)

    return colormap

def _plot_heatmap(
        data : pd.DataFrame, 
        vmin : float, 
        vmax : float, 
        cmap = None, 
        norm = None,
        ax : plt.Axes = None, 
        log : bool=False,
        ) :
    
    if ax is None :
        fig = plt.figure(figsize=(12,10))
        ax = fig.gca()
    ax.invert_yaxis()
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
    
    
    if cmap is None : 
        cmap = get_cmap('inferno')
        cmap.set_bad('gray')
    elif isinstance(cmap, str) :
        cmap = get_cmap(cmap)
        cmap.set_bad('gray')
    
    if norm is None :
        if log : 
            norm = LogNorm(vmin=vmin, vmax=vmax)
    
    colormesh = ax.pcolormesh(
        data,
        norm=norm,
        cmap=cmap,
        edgecolor = 'black',
        )

    x_pos = np.arange(len(data.columns)) + 0.5
    y_pos = np.arange(len(data.index)) + 0.5
    ax.set_xticks(x_pos, data.columns, rotation = 30)
    ax.set_yticks(y_pos, data.index)

    cbar : matplotlib.colorbar.Colorbar = plt.colorbar(colormesh, ax=ax, norm=norm)
    min_value = norm.inverse(0)
    mid_value = norm.inverse(0.5)
    int_value = norm.inverse(0.75)
    max_value = norm.inverse(1)

    values = [min_value, mid_value,int_value, max_value]
    if 1 not in values :
        values.append(1)

    cbar.set_ticks(values, labels= values)


    if ax is None :
        return fig
    else :
        return ax
    
def create_pair_colocalisation_figure(
        colocalization_rates : pd.DataFrame,
        zscore_frame : pd.DataFrame,
        pvalue_mask : pd.DataFrame,
        frameon = True,
) :
    
    if len(colocalization_rates.index) != len(colocalization_rates.columns) : raise ValueError("Colocalization rate df is not squared. Did you group by target first ?\n{}".format(colocalization_rates))
    if len(zscore_frame.index) != len(zscore_frame.columns) : raise ValueError("Colocalization rate df is not squared. Did you group by target first ?\n{}".format(zscore_frame))

    cmap = score_color_scale(vmax=20)
    cmap.set_bad('white')

    score_norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=20,)

    #Unsgnificative pvalue set to NaN
    zscore_frame *= pvalue_mask.replace({True : 1, False : np.NaN})
    colocalization_rates *= pvalue_mask.replace({True : 1, False : np.NaN})

    fig = plt.figure(figsize=(24,10), frameon=frameon)
    left,right = fig.subplots(1,2)
    
    right.set_title("Z-score (standardised normalistion)")
    right = _plot_heatmap(zscore_frame, vmin=-1, vmax=20, ax=right, cmap=cmap, norm=score_norm)
    #Hatch nan values
    nan_mask = np.isnan(zscore_frame.to_numpy())
    nan_mask = np.ma.masked_where(~nan_mask, np.ones_like(zscore_frame.to_numpy()))
    right.pcolor(
    nan_mask,
    hatch='///',
    alpha=0.0,  # transparent so you only see the hatching
    edgecolor='none',
    linewidth=0,
    cmap='Greys'
)
    left.set_title("Co-localization rate")
    left = _plot_heatmap(colocalization_rates, vmin=1e-4,log=True, vmax=1, ax=left, cmap=cmap)
    left.pcolor(
    nan_mask,
    hatch='///',
    alpha=0.0,  # transparent so you only see the hatching
    edgecolor='none',
    linewidth=0,
    cmap='Greys'
)
 
    left.set_ylabel('Wich fraction of ..')
    left.set_ylabel('Co-localizes with ..')

    
    return fig


"""

6. Main functions

"""
import os, logging, traceback

def main(
        filtered_Spots : pd.DataFrame,
        Cell : pd.DataFrame,
        Acquisition : pd.DataFrame,
        Detection : pd.DataFrame,
        Gene_map : pd.DataFrame,
        colocalisation_distance : int,
        run_path : str,
        significance : float = 1e-4,
        frameon = True, 
) :
    output_path = run_path + "/analysis/co_localization/"
    os.makedirs(output_path, exist_ok=True)
    
    log_file = output_path + "/error_log.log"
    logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    force= True
    )
    error_count = 0

    try :
        logging.info(f"Pairwise co-localization analysis start")
        sucess = pairwise_colocalization_analysis(
            filtered_Spots=filtered_Spots,
            Cell=Cell,
            Acquisition=Acquisition,
            Detection=Detection,
            Gene_map=Gene_map,
            colocalisation_distance=colocalisation_distance,
            output_path=output_path,
            significance=significance,
            frameon=frameon
        )

    except Exception as e :
        logging.error(f"Pairwise co-localization analysis failed :\n{traceback.format_exc()}")
        error_count += 1
    else :
        logging.info("Pairwise co-localization sucess")

    if error_count == 0 : 
        return True
    else :
        return False

def pairwise_colocalization_analysis(
        filtered_Spots : pd.DataFrame,
        Cell : pd.DataFrame,
        Acquisition : pd.DataFrame,
        Detection : pd.DataFrame,
        Gene_map : pd.DataFrame,
        colocalisation_distance : int,
        output_path : str,
        significance : float = 1e-4,
        frameon = True,
) :
    
    filtered_Spots = _spots_merge_data(
        Spots=filtered_Spots,
        Acquisition=Acquisition,
        Detection=Detection,
        Gene_map=Gene_map,
        Cell=Cell
    )

    voxel_size = Detection['voxel_size'].at[0]

    RNA_list = list(filtered_Spots['target'].unique())
    RNA_list.sort()

    #Coloc rates from models
    coloc_rates, selfcoloc_rates = create_coloc_rate_expectancy(
        Spots=filtered_Spots,
        Cell=Cell,
        voxel_size=voxel_size,
        colocalisation_distance=colocalisation_distance,
        RNA_list=RNA_list
    )

    # Model values for expectancy and std
    abundancies = _compute_cell_distribution_populations(filtered_Spots)
    cell_ids = coloc_rates.index
    multi_index = pd.MultiIndex.from_product([RNA_list, abundancies.index]) #one line per couple (rna,cell_id)
    multi_index.names = ['target','cell_id']

    ## Initialise empty dataframes for expectancy/std
    expected_event_count = pd.DataFrame(columns= RNA_list, index=multi_index, dtype=float)
    expected_event_count_std = pd.DataFrame(columns= RNA_list, index=multi_index, dtype=float)

    ## Filling
    for rna in RNA_list :
        product = coloc_rates.multiply(abundancies[rna],axis=0)
        product_index = pd.MultiIndex.from_product([[rna], cell_ids])

        ###Expectancy
        expected_event_count.loc[product_index, :] = product.values # E = n*p
        expected_event_count.loc[rna,[rna]] = (selfcoloc_rates[rna] * abundancies[rna]).values # Correction for selfcoloc

        ###Std
        product = coloc_rates.multiply(abundancies[rna],axis=0).multiply((1-coloc_rates),axis=0) #std = sqrt(np(1-p))
        product = product.apply(np.sqrt)
        product_index = pd.MultiIndex.from_product([[rna], cell_ids])
        expected_event_count_std.loc[product_index, :] = product.values
    
    # Coloc measurements
    colocalisation_truth = colocalisation_truth_df(
        Spots=filtered_Spots,
        Acquisition=Acquisition,
        Detection=Detection,
        Gene_map=Gene_map,
        Cell=Cell,
        colocalisation_distance=colocalisation_distance
    )

    colocalisation_truth = safe_merge_no_duplicates(
        colocalisation_truth,
        filtered_Spots,
        on='spot_id',
        keys='cell_id'
    )

    measure_coloc_events = colocalisation_truth.groupby(['target','cell_id'])[RNA_list].sum()
    coloc_rates = measure_coloc_events.divide(abundancies)
    
    #Zscore computation
    zscore_frame = compute_z_score_frame(
        measured_colocalisation_events=measure_coloc_events,
        expected_colocalisation_events=expected_event_count,
        expected_standard_deviation= expected_event_count_std,
    )

    #Save datasheet
    os.makedirs(output_path + "/datasheet/",exist_ok=True)
    mean_coloc_rates = coloc_rates.groupby('target',axis=0,level=0).mean()
    mean_coloc_rates.to_excel(output_path + "/datasheet/coloc_rates_mean.xlsx")
    median_zscore = zscore_frame.groupby('target',axis=0,level=0).median()
    median_zscore.to_excel(output_path + "/datasheet/zscore.xlsx")
    
    #p-values computation
    pvalue_frame = compute_pvalue_frame(
        zscore_frame=zscore_frame
    )
    pvalue_mask = pvalue_frame <= significance
    
    #Create graph
    pairwise_coloc_fig = create_pair_colocalisation_figure(
        colocalization_rates= mean_coloc_rates,
        zscore_frame= median_zscore,
        pvalue_mask=pvalue_mask,
        frameon=frameon
    )

    #Save graph
    pairwise_coloc_fig.savefig(output_path + "/pairwise_colocalisation_heatmap.svg")
    plt.close()

    return True