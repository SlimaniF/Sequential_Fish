import pandas as pd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.container as container
import functools
import itertools
from itertools import zip_longest, chain
from math import ceil
from ..tools import safe_merge_no_duplicates

def get_xlabels(ax : plt.Axes) :
    xtickslabels = ax.get_xticklabels()
    xlabels = [text.get_text() for text in xtickslabels]
    
    return xlabels

def get_ylabels(ax : plt.Axes) :
    ytickslabels = ax.get_yticklabels()
    ylabels = [text.get_text() for text in ytickslabels]
    
    return ylabels
    

def _get_min_cluster_radius(voxel_size) :
    return max(voxel_size)
    
    
def _get_voxel_size(Detection : pd.DataFrame) :
    voxel_size = tuple(Detection['voxel_size'].iat[0])
    voxel_size = [int(i) for i in voxel_size]
    
    return voxel_size

def merge_data(
    Acquisition : pd.DataFrame,
    Detection : pd.DataFrame,
    Cell : pd.DataFrame,
    Spots : pd.DataFrame,
    Gene_map : pd.DataFrame,
) :
    """
    Returns : Detection, Cell, Spots
    """
    Detection = safe_merge_no_duplicates(
        Detection,
        Acquisition,
        on= 'acquisition_id',
        keys= ['cycle']
    )
    
    Detection = safe_merge_no_duplicates(
        Detection,
        Gene_map,
        on= ['cycle', 'color_id'],
        keys= 'target'
    )
    
    
    Spots = safe_merge_no_duplicates(
        Spots,
        Detection,
        on='detection_id',
        keys= ['target', 'location']
    )


    Cell = safe_merge_no_duplicates(
        Cell,
        Detection,
        on='detection_id',
        keys='target'
    )
    
    return Detection, Cell, Spots




## PLOTS


def distribution_super_plot(data, ax, title=None, xlabel=None, ylabel=None, showextrema=False, **kwargs) :
    """
    Distribution plot for pd.Series data.
    Series is expected to have a multi-index from 1 to 3 dimensions and Series values are expected to be lists.

    Super plot is returned as plt.Axes

    """
    if type(data) != pd.Series : raise TypeError('data argument passed is of type {0} should be "pd.Series".'.format(type(data)))
    if type(ax) != plt.Axes : raise TypeError('ax argument passed is of type {0} should be "plt.Axes".'.format(type(ax)))

    level = len(data.index.names)

    if level == 1 :
        ax, legend = _distribution_lvl1(data, ax, showextrema=showextrema, **kwargs)
    
    elif level == 2 :
        ax, legend = _distribution_lvl2(data, ax, showextrema=showextrema, **kwargs)

    elif level == 3 :
        ax, legend = _distribution_lvl3(data, ax, showextrema=showextrema, **kwargs)
    
    else : raise LevelError("Unsupported number of dimension in index (should be between 1 and 3).")
    
    legend_size = len(legend[0])
    if legend_size > 5 :
        ncols= (legend_size // 5) + 1
        legend = ax.legend(*legend, ncols=ncols)
    else : legend = ax.legend(*legend)
    
    if type(title) != type(None) : ax.set_title(title)
    if type(xlabel) != type(None) : ax.set_xlabel(xlabel)
    if type(ylabel) != type(None) : ax.set_ylabel(ylabel)

    return ax


class LevelError(IndexError) :
    pass

def _distribution_lvl1(data: pd.Series, ax: plt.Axes, showextrema=False, **kwargs) :

    multi_index = data.index.names
    if len(multi_index) != 1 : raise LevelError("_distribution_lvl1 was called but multi-index dimension does not match.")
    if type(data.iat[0]) != list : raise ValueError("for distribution plot a list is expected for each distributions to plot. {0} was found.".format(type(data.iat[0])))

    colors = get_colors_list(len(data))
    labels = list(data.index)

    ax = violin_plot(
        ax=ax,
        distributions=data,
        labels= labels,
        colors= colors,
        showmeans= True,
        mean_size= 90,
        showextrema=showextrema,
    )

    legend_prep =(
        [plt.scatter(0,0, c= 'white', edgecolors= 'black')],
        ['distribution means']
        )

    return ax, legend_prep

def _distribution_lvl2(data: pd.Series, ax : plt.Axes, alpha= 0.6, showextrema=False, show_distribution_size=True, **kwargs) :

    # print(data)

    # if kwargs.get('sort') :
    if 'sort_parameters' in kwargs.keys() : sort_parameters = kwargs['sort_parameters']
    else : sort_parameters = {}
    #     data = data.sort_index(**sort_parameters)
    # print(data)
    
    data = data.sort_index(**sort_parameters)
    multi_index = data.index.names
    
    if len(multi_index) != 2 : raise LevelError("_distribution_lvl2 was called but multi-index dimension does not match.")
    
    measure = data.name
    distributions: pd.Series = data.reset_index(drop=False).groupby(multi_index[:1], sort=False)[measure].agg(list)
    
    labels_lvl2 = list(distributions.index)
    if show_distribution_size :
        for index, label in enumerate(labels_lvl2) : 
            labels_lvl2[index] = label + '\n' + str([len(dis) for dis in distributions.iat[index]])
    labels_lvl1 = list(data.index.get_level_values(1).unique())
    
    custom_colors = kwargs.get('colors')
    colors_df = make_color_frame(labels= labels_lvl1, custom_colors=custom_colors)

    if 'ascending' in sort_parameters.keys() : ascending = sort_parameters['ascending']
    else : ascending = True
    colors = list(pd.merge(data.reset_index(), colors_df, left_on= multi_index[1], right_on='labels').sort_values(multi_index, ascending=ascending)['colors'])

    ax = violin_plot(
        ax=ax,
        distributions=distributions,
        labels= labels_lvl2,
        colors=colors,
        alpha=alpha,
        linewith=1,
        showmeans= True,
        mean_size= 90,
        showextrema=showextrema,
        multi_violin_plot= True
    )

    legend_prep = (
        [plt.scatter(0,0, c= 'white', edgecolors= 'black')] + [plt.bar(0,0, color= color) for color in colors_df['colors']],
        ['distribution means'] + list(colors_df.index)        
    )


    return ax, legend_prep


def _extract_color_from_legend(legend_prep, remove_first=False) :
    if remove_first :
        legend_prep = (legend_prep[0][1:], legend_prep[1][1:])
    
    artists: list[container.BarContainer] = legend_prep[0]
    colors = [(artist.patches[0].get_facecolor(),) for artist in artists]

    color_df = pd.DataFrame(data= colors, index= legend_prep[1], columns= ['colors'])
    return color_df

def _distribution_lvl3(data: pd.Series, ax: plt.Axes, showextrema, **kwargs) :
    """

    yaxis = measure_value

    level1 = x-axis (axis label)
    level2 = distribution colors (legend)
    level3 = smaller scatter points within distribution whith different shapes (legend)

    """
    multi_index = data.index.names
    if len(multi_index) != 3 : raise LevelError("_distribution_lvl3 was called but multi-index dimension does not match.")
    data = data.sort_index()
    measure = data.name

    #Levevel 1&2
    list_flattener = functools.partial(sum, start=[])
    data_distribution_lvl2: pd.Series = data.reset_index(drop=False).groupby(multi_index[:2])[measure].apply(list).apply(list_flattener)
    ax, legend_lvl2 = _distribution_lvl2(
        data_distribution_lvl2, 
        ax,
        alpha=0.4,
        showextrema=showextrema,
        **kwargs
        )

    #Level 3
    index_lvl3 = data.index.get_level_values(2).unique()

    gen = get_markers_generator()
    next(gen) #skipping 'o' that is used for means of distributions
    marker_list = [next(gen) for i in range(len(index_lvl3))]
    marker_frame = pd.Series(data= marker_list, index= index_lvl3, name= 'markers')
    color_df = _extract_color_from_legend(legend_lvl2, remove_first= True) #First artist is always distributions mean.

    markers = pd.merge(data.reset_index(), marker_frame, left_on= multi_index[2], right_on=multi_index[2]).sort_values(multi_index)['markers']
    colors = pd.merge(data.reset_index(), color_df, left_on= multi_index[1], right_index=True ).sort_values(multi_index)['colors']

    positions = [[i+1]*len(data.xs(idx,axis=0,level=1)) for i,idx in enumerate(data.index.get_level_values(1).unique())]
    print(positions)
    positions = sum(positions,[])

    #random span for points
    seed = 0
    random_gen = np.random.default_rng(seed=seed)
    step_between_positions = 0.20

    assert len(markers) == len(colors) == len(positions) == len(data), 'all length should be identical :\nmarkers length : {0}\ncolors length : {1}\npositions length {2}\ndata length {3}'.format(len(markers), len(colors), len(positions), len(data))

    for position, distribution_lvl3, marker, color in zip(positions, data, markers, colors) :
        x = np.array([position] * len(distribution_lvl3), dtype=float)
        x += (random_gen.random(size=len(distribution_lvl3), dtype=float) * step_between_positions * random_gen.choice([1,-1], size=len(distribution_lvl3)))
        ax.scatter(
            x= x,
            y= distribution_lvl3,
            s = 25,
            color = color,
            marker = marker,
            alpha= 0.1,
            zorder =-100
        )

    legend_prep = (
        [plt.scatter(0,0, c= 'white', edgecolors= 'black')]  + [plt.bar(0,0, color= color) for color in color_df['colors']],
        ['distributions means'] + list(color_df.index)
    )

    return ax, legend_prep

def get_colors_list(size:int = 100, remove_black= False, remove_grey = False, remove_brown= False) -> list:
    """
    Get a list of color from matplotlib.colors of length 'size'.
    100 different shade in the library
    """
    if not isinstance(size, int) : raise TypeError("size should be an int, it is a {0}".format(type(size)))
    if size < 1 : raise ValueError("Size should be >= 1")

    red = _get_red_colors()
    yellow = _get_yellow_colors()
    green = _get_green_colors()
    blue = _get_blue_colors()
    purple = _get_purple_colors()
    brown = _get_brown_colors()
    pink = _get_pink_colors()
    orange = _get_orange_colors()
    black = _get_black_colors()
    grey = _get_grey_colors()

    color_list = list(sum([*zip_longest(red, green, blue, black, orange, purple, grey, yellow, brown, pink)],()))
    length = len(color_list)
    while None in color_list : color_list.remove(None)
    if size > length :
        iteration = ceil(size / length)
        return (color_list * iteration)[:size]

    if remove_black : 
        for color in _get_black_colors() : color_list.remove(color)
    if remove_grey :
        for color in _get_grey_colors() : color_list.remove(color)
    if remove_brown :
        for color in _get_brown_colors() : color_list.remove(color)

    return (color_list)[:size]

def get_color_generator(remove_black=False, remove_grey = False, remove_brown= False) :
    """
    Cycled color generator (100 shades)
    """
    color_list = get_colors_list(remove_black=remove_black, remove_grey=remove_grey, remove_brown=remove_brown)

    color_gen = itertools.cycle(color_list)

    return color_gen


def _get_red_colors() :
    return ["#D0312D", "#990F02", "#60100B", "#7E2811", "#4E0707", "#BC544B", "#680C07"]

def _get_orange_colors():
    return ["#ED7014", "#FCAE1E", "#B56727", "#BE5504", "#D67229", "#E34A27", "#FF8C00"]

def _get_yellow_colors():
    return ["#D6B85A", "#DFC98A", "#C8A951", "#E7C27D", "#BDA55D", "#E4D00A", "#FFEF00"]

def _get_green_colors():
    return ["#3CB043", "#3A5311", "#728C69", "#AEF359", "#5DBB63", "#028A0F", "#234F1E", "#568203", "#4CBB17", "#487800"]

def _get_blue_colors():
    return ["#3944BC", "#63C5DA", "#0A1172", "#281E5D", "#1338BE", "#48AAAD", "#016064", "#2832C2", "#1F456E", "#4682B4"]

def _get_purple_colors():
    return ["#A32CC4", "#7A4988", "#601A35", "#A1045A", "#663046", "#311432", "#9867C5", "#880085"]

def _get_pink_colors():
    return ["#FC94AF", "#F25278", "#FE7D6A", "#FD5DA8", "#E3256B", "#FF6EC7"]

def _get_brown_colors():
    return ["#4B371C", "#231709", "#795C34", "#CC7722", "#65350F", "#652A0E", "#8A3324"]

def _get_black_colors():
    return ["#000000"]

def _get_grey_colors():
    return ["#808080", "#373737", "#594D5B", "#3E3D53", "#9897A9", "#63645E"]

def violin_plot(
        ax: plt.Axes, distributions, labels=None, sub_labels=None, colors=None, xlabel=None, ylabel=None, title=None, y_axis= None,
        linewith = 2, line_color = 'black', alpha= 0.6,
        vertical_plot=True, showmedians= True, showextrema = False, showmeans=False, mean_size = 45, multi_violin_plot= False
                ) :

    """
    Basis for violin plots.
    """

    if type(labels) == type(None) : pass
    elif len(distributions) != len(labels) : raise ValueError("Length of labels and distributions must match.")

    if type(colors) == type(None): pass
    elif len(distributions) != len(colors) and not multi_violin_plot: raise ValueError("Length of colors and distributions must match.")
    
    if type(y_axis) == type(None) : pass
    elif len(y_axis) != 2 : raise ValueError("y_axis is expected to be of length 2 : (ymin, ymax).")
    
    if multi_violin_plot :
        if type(colors) != type(None) :
            if len(colors) == len(distributions) :
                new_colors = []
                for color,distrib in zip(colors, distributions) :
                    new_colors.extend([color]*len(distrib))
                colors = new_colors

        max_individual_violin_number = max([len(distrib) for distrib in distributions]) + 1 #Is the maximum number of violin plotted for one set.
        positions, ticks_positions = multi_violin_plot_positions(distributions)
        distributions = list(chain(*distributions))
        assert len(distributions) == len(positions), "AssertionError : multi_distributions wrongly flattened : positions : {0}, distributions {1}".format(len(positions), len(distributions))

        #colors
        if type(colors) != type(None) :
            if len(colors) == len(distributions) : pass
            else : raise ValueError("Length of colors must either match length of distributions or the number of element in distributions")
    
    else :
        positions = np.arange(1, len(distributions) + 1)
        ticks_positions = np.arange(1, len(distributions) + 1)
        max_individual_violin_number = 1
    #Plot
    violin_plot = ax.violinplot(
        distributions,
        positions=positions, 
        vert= vertical_plot, 
        showmedians=showmedians,
        showextrema= showextrema,
        )

    if type(labels) == type(None) :
        labels = np.arange(1, len(distributions) + 1)
    
    xticks = ax.set_xticks(ticks_positions, labels=labels)


    ax.set_xlim(0.25, len(labels) * max_individual_violin_number + 0.75)
    if type(colors) == type(None) : colors = get_colors_list(len(violin_plot['bodies']))
    for violin, color in zip(violin_plot['bodies'], colors) :
        violin.set_facecolor(color)
        violin.set_alpha(alpha)

    for collection_name in ['cbars', 'cmins', 'cmaxes', 'cmedians'] :
        collection = violin_plot.get(collection_name)
        if type(collection) != type(None) :    
            collection.set_color(line_color)
            collection.set_linewidth(linewith)
    
    if type(y_axis) != type(None) :
        axis = list(ax.axis())
        if y_axis[0] != None : axis[2] = y_axis[0]
        if y_axis[1] != None : axis[3] = y_axis[1]
        ax.axis(axis)
    else : axis = ax.axis()

    if type(sub_labels) != type(None) :
        if len(sub_labels) != len(distributions) : raise ValueError("Length of sub_labels must match number of violins to plot.")
        for sub_label, x_position in zip(sub_labels, positions) :
            ax.text(x = x_position, y= axis[2], s= sub_label, ha= 'center')
    
    if showmeans :
        means = [np.mean(distrib) for distrib in distributions]
        ax.scatter(positions, means, c= colors, s= mean_size, linewidths=0.5, edgecolors='black')

    if type(xlabel) != type(None) : ax.set_xlabel(xlabel)
    if type(ylabel) != type(None) : ax.set_ylabel(ylabel)
    if type(title) != type(None) : ax.set_title(title)

    return ax

def multi_violin_plot_positions(distributions) :

    max_individual_violin_number = max([len(distrib) for distrib in distributions]) + 1#Is the maximum number of violin plotted for one set.

    positions = []
    ticks_positions = []
    for distrib_number, distrib in enumerate(distributions) :
        positions.extend(list(
            np.arange(1, len(distrib) + 1) + (distrib_number * max_individual_violin_number) if len(distrib) > 1 
            else [distrib_number * max_individual_violin_number + (max_individual_violin_number-1)/2 + 1]
        ))

        ticks_positions.append(
            distrib_number * max_individual_violin_number + (len(distrib)-1)/2 + 1 if len(distrib) > 1
            else distrib_number * max_individual_violin_number + (max_individual_violin_number-1)/2 + 1
        )

    return positions, ticks_positions

def get_markers_generator() :
    markers = (marker for marker in ['o','v','D','x','<','>','s','8','p','*','h','P','^'])
    gen = itertools.cycle(markers)
    
def make_color_frame(labels : 'list[str]', custom_colors=[]) :
    """
    Return a color df where labels are passed in index and unique columns 'colors' contains color values. 
    """

    if isinstance(custom_colors, (list,tuple) ) :
        if len(custom_colors) == len(labels) :
            colors = custom_colors
        else : colors = get_colors_list(len(labels))
    else : colors = get_colors_list(len(labels))
    
    
    color_df = pd.DataFrame({
        'labels' : labels,
        'colors' : colors,
    })
    
    return color_df.set_index('labels')