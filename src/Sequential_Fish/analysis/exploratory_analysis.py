"""
Submodule handling multivariate exploratory analysis.  
2 studies are prepared one for co-localization exploration using MCA (Multiple component analysis) and one for RNAs expression using PCA (Principle component analysis).  
Refer to the documentation to read about purpose and interpretations.
"""

import os
import warnings
import traceback
import logging
from itertools import cycle
from itertools import combinations
from itertools import permutations
from typing import Literal

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib import gridspec
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.patches import Patch
import prince
import pingouin as pg
import networkx as nx
import igraph as ig
import leidenalg as la
from scipy.spatial.distance import pdist, squareform
from scipy.stats import chi2
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale



from .colocalisation import colocalisation_truth_df

def main(
        run_path : str,
        Spots : pd.DataFrame,
        colocalisation_distance : int,
        control_genes : list[str] | None,
        thresold_coloc_value : float,
        threshold_coloc_zscore : float,
) :

    exploration_path = os.path.join(run_path,"analysis","multivariate","exploration")
    os.makedirs(exploration_path,exist_ok=True)

    try :
        logging.info("Starting colocalization exploratory analysis.")
        zscore_path = os.path.join(run_path,"analysis","co_localization","datasheet","zscore.csv")
        if not os.path.isfile(zscore_path) : 
            raise FileNotFoundError("Couldn't find colocalization zscore table. Run first co-localization analysis.")

        zscores = pd.read_csv(zscore_path,sep=";")
        zscores = zscores.set_index(zscores.columns[0])


        colocalization_figure = colocalization_exploration(
            Spots, 
            colocalisation_distance, 
            run_path,
            zscores=zscores,
            threshold_coloc_rate=thresold_coloc_value,
            threshold_zscore=threshold_coloc_zscore
            )
        colocalization_figure.savefig(os.path.join(exploration_path,"colocalization.svg"))

    except Exception :
        warnings.warn("Exception raised during multivariate colocalization analysis. Check log")
        logging.error(msg=traceback.format_exc())
        return False

    try :
        logging.info("Starting expression exploratory analysis.")
        expression_exploration_figures_dict = expression_exploration(Spots, run_path, control_genes=control_genes)
        for fig_name, fig in expression_exploration_figures_dict.items() :
            fig.savefig(os.path.join(exploration_path,fig_name + ".svg"))

    except Exception :
        warnings.warn("Exception raised during multivariate expression analysis. Check log")
        logging.error(msg=traceback.format_exc())
        return False

    
        

    return True


def colocalization_exploration(
        Spots : pd.DataFrame,
        colocalisation_distance : int,
        run_path : str,
        zscores : pd.DataFrame,
        threshold_coloc_rate : float,
        threshold_zscore : float,
        na_policy : Literal["fill","drop"] = "fill",
) :
    """
    Perform exploratory analysis to identify structures in co-localization data.

    # Input  
        **Spots** : pd.Dataframe, Spots table after post processing (filter and merge data)  
        **colocalization_distance** : int  

    # Returns  
        **fig** : plt.Figure  

    ## figure  
      
    * Multiple component analysis (MCA) : Scree plot + MCA intertia composition  
    * Hierarchical clustering : Dendogram  
    * Network representation with Louvain communities  
    """

    if "target" not in Spots.columns : 
        raise KeyError("'target' not found in Spots columns. Did you run Spots_post_processing ?")

    data_path = os.path.join(run_path,"analysis","coloc_truth_table.csv")
    if os.path.isfile(data_path) :
        print(f"Using cached data at {data_path}. If you modified user parameters delete cached data.")
        logging.info(f"Using cached data at  : {data_path}")
        data = pd.read_csv(data_path, sep=";").iloc[:,1:]
    else :
        data = _make_data_table_colocalization(Spots, colocalisation_distance)
        data.to_csv(data_path,sep=";")

    if na_policy == "fill" : # Na values when a distribution was not detected in cell
        data = data.fillna(False) 
    elif na_policy == "drop" :
        ref_len = len(data)
        data = data.dropna()
        logging.info(f"{ref_len - len(data)} observation dropped according to 'drop' na_policy.")
    else :
        raise ValueError(f"na_policy should be either 'drop' or 'fill' got {na_policy}")

    # Init figure
    fig = plt.figure(figsize=(18,18))
    grid = gridspec.GridSpec(5,5)
    scree_ax = fig.add_subplot(grid[0,:])
    mca_intertia_ax = fig.add_subplot(grid[1,:3])
    network_graph_ax = fig.add_subplot(grid[2:,:3])
    dendogram_ax = fig.add_subplot(grid[1:,3:])

    #1.Multiple Compnent analysis
    scree_ax, mca_intertia_ax = _multiple_component_analsis(
        data=data,
        scree_ax=scree_ax,
        mca_intertia_ax=mca_intertia_ax
    )

    #2. Network_plots
    network_graph_ax = _network_analysis(
        data, 
        zscores=zscores,
        network_graph_ax=network_graph_ax,
        threshold_value=threshold_coloc_rate,
        threshold_zscore=threshold_zscore
        )


    #3. Hierarchical clustering
    dendogram_ax = _hierarchical_clustering_analysis(data,dendogram_ax)

    return fig

def _make_data_table_colocalization(
        Spots : pd.DataFrame,
        colocalisation_distance : int,
) :
    data = colocalisation_truth_df(
            Spots=Spots,
            colocalisation_distance=colocalisation_distance
        )

    data = data.drop(columns=data.columns[data.columns.str.contains("endcycle")])
    rnas_col = data.columns.tolist()[:-4]

    for rna in rnas_col :
        data.loc[data["target"] == rna, rna] = True # Enforce self co-localization to always be True, this is critical to identify structure during exploratory analysis. Self co-localization or tendency to form monodistribution cluster must be studied elsewhere.

    return data.iloc[:,:-4]

#1.
def _multiple_component_analsis(
        data : pd.DataFrame,
        scree_ax : Axes,
        mca_intertia_ax : Axes
) :

    if not (data.dtypes == bool).all() : 
        raise TypeError("All columns must be booleans did you index columns with interest RNAS ?")
    mca = _run_mca(data)
    inertia_contribution = _compute_intertia_contribution(mca)

    scree_ax = _make_scree_plot(inertia_contribution, scree_ax)
    mca_intertia_ax = _make_contribution_plot(mca, inertia_contribution, mca_intertia_ax)

    return scree_ax, mca_intertia_ax

def _run_mca(
        data : pd.DataFrame
) :

    # Decomposition
    mca = prince.MCA(n_components=len(data.columns)-1)
    mca = mca.fit(data)

    return mca

def _compute_intertia_contribution(mca : prince.MCA) :
    inertia_contribution = pd.DataFrame(
        data= [mca.eigenvalues_ / np.sum(mca.eigenvalues_)],
        columns= mca.column_contributions_.columns,
        ).T.reset_index(drop=False)

    return inertia_contribution

def _find_elbow(intertia_contribution : pd.DataFrame) :

    first_derivative = np.diff(intertia_contribution[0])
    second_derivative = np.diff(first_derivative)

    # Find the elbow point (where the second derivative is largest)
    elbow_point = np.argmax(second_derivative) + 1  # +1 because diff reduces the array size

    return elbow_point

def _make_scree_plot(
        inertia_contribution : pd.DataFrame,
        ax : Axes
        ) :

    ax.plot(
    inertia_contribution["index"],
    inertia_contribution[0],
    "-",
    c="black",
    )

    ax.plot(
        inertia_contribution["index"],
        inertia_contribution[0],
        "o",
        c="lightblue",
        mec="black",
        ms = 10
    )

    ax.set_xticks(inertia_contribution["index"],  inertia_contribution["index"])
    ax.set_ylabel("Contribution to intertia (fraction)")
    ax.set_xlabel("Dimension")
    ax.set_title("Multiple component analysis")
    sns.despine(ax=ax)

    return ax

def _make_contribution_plot(
        mca : prince.mca,
        inertia_contribution : pd.DataFrame,
        ax : Axes
) :

    elbow_point = _find_elbow(inertia_contribution)
    contributions = mca.column_contributions_
    contributions = contributions.iloc[:,:elbow_point+1]
    contributions.index = pd.MultiIndex.from_arrays(zip(*contributions.index.str.split("__").to_list()))
    contributions = contributions.reset_index(drop=False).sort_values(["level_1","level_0"])

    #Colormap
    colors = ["#999999", "#ffffbf", "#1a9641"]
    values = [0, 0.5, 1]  # Adjust as needed for smooth transitions
    cmap = LinearSegmentedColormap.from_list("custom_log_cmap", list(zip(values, colors)))
    norm = LogNorm(vmin=1e-2, vmax=1)

    ax = sns.heatmap(
        contributions.iloc[:,2:].T, 
        annot=False, 
        vmax=1, 
        vmin=0,
        ax=ax,
        linewidths=.5,
        square=True,
        cbar=True,
        norm=norm,
        cmap=cmap,
        cbar_kws={
                "orientation" : "horizontal",
                "pad" : .01
            }
        )

    #Ticks and labels
    ax.xaxis.set_ticks_position("top")
    RNA = contributions["level_0"]
    colors = contributions["level_1"]
    colors = ["#2ca25f" if boolean =="True" else "#e6550d" for boolean in colors]
    ax.set_xticks(ax.get_xticks(), RNA, rotation=30, ha="left")
    for label, color in zip(ax.get_xticklabels(), colors) :
        label.set_color(color)
    
    return ax

#2.
def _network_analysis(
        data : pd.DataFrame,
        zscores : pd.DataFrame,
        network_graph_ax : Axes,
        threshold_zscore : float = 1,
        threshold_value : float =.1,
) :

    G, edges = _construct_network(data, zscores, threshold_zscore=threshold_zscore, threshold_value=threshold_value)
    partition = _find_communities(G)
    network_graph_ax = _make_network_plot(
        G,
        edges,
        network_graph_ax,
        data,
        partition=partition,
        )

    return network_graph_ax

def _construct_network(
        data : pd.DataFrame, 
        zscores:pd.DataFrame,
        threshold_zscore : float = 1,
        threshold_value: float = .1 ,
        ) :
    G = nx.DiGraph()
    edges = {
        "zscores" : {"pairs" : [], "width" : []},
        "value" : {"pairs" : [], "width" : [], "coloc_rate" : []},
    }

    # Add nodes (RNAs)
    for rna in data.columns:
        G.add_node(rna)

    for rna1, rna2 in permutations(data.columns.to_list(),r=2) :

        if rna1 == rna2 : continue
        zscore = zscores.at[rna1,rna2]

        if abs(zscore)>=threshold_zscore  : 
            edges['zscores']["pairs"].append((rna1,rna2))
            edges['zscores']["width"].append(zscore / 50 *20)

        coloc_number = np.sum(([data[rna1] & data[rna2]]))
        population_size = np.sum(data[rna1])

        coloc_rate = coloc_number/population_size
        if coloc_rate >= threshold_value : 
            edges["value"]["pairs"].append((rna1,rna2))
            edges["value"]["width"].append(max(.5,coloc_rate*10))
            edges["value"]["coloc_rate"].append(str(round(coloc_rate*100, 2)) + " %" if coloc_rate >0.01 else "")

        if coloc_rate >= threshold_value or abs(zscore)>=threshold_zscore :
            G.add_edge(rna1,rna2, weight=coloc_rate)

    return G, edges

def _make_network_plot(
        G : nx.Graph,
        edges : dict[Literal["zscores", "value"], tuple],
        ax : Axes,
        data : pd.DataFrame,
        partition : dict | None = None,
) :

    # Assign colors to communities
    if not partition is None :
        colors_com = [partition[node] / len(G.nodes()) for node in G.nodes()]
    else : 
        colors_com = None


    pos = nx.spring_layout(G, seed=1)
    colors = cycle(['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f'])
    colors = [next(colors) for i in range(len(G.nodes))]

    cmap = LinearSegmentedColormap.from_list(
        "network_qualitative_map", 
        list(zip(np.linspace(0,1,len(colors)),colors))
        )

    nx.draw_networkx_edges(
        G,pos=pos, 
        edgelist=edges["zscores"]["pairs"], 
        width= edges["zscores"]["width"],
        edge_color='red',
        )
    nx.draw_networkx_edge_labels(
        G,pos,
        edge_labels= dict(zip(
            edges["value"]["pairs"],
            edges["value"]["coloc_rate"]
            )),

        verticalalignment="bottom"
    )

    
    nx.draw_networkx_edges(
            G,pos=pos, 
            edgelist=edges["value"]["pairs"], 
            width= edges["value"]["width"],
            edge_color=["black" if val > .5 else "gray" for val in edges["value"]["width"]],
            alpha=.3,
            style=':'
            )


    max_pop = data.sum(axis=0).max()

    nx.draw_networkx_nodes(
        G, pos, 
        node_color=colors if partition is None else colors_com, cmap=None if partition is None else cmap,
        node_size= [np.sum(data[node]) / max_pop *2000 for node in G.nodes], 
        edgecolors="black", alpha=.5,
        )
    nx.draw_networkx_labels(
        G, pos,
        )

    return ax


def _find_communities(G : nx.DiGraph) :
    edges = [(u, v, d.get("weight", 1.0)) for u, v, d in G.edges(data=True)]
    igG = ig.Graph.TupleList(edges, directed=True, edge_attrs=["weight"])

    part = la.find_partition(
        igG,
        la.ModularityVertexPartition,        # or RBConfigurationVertexPartition
        weights="weight",
        n_iterations=-1,
    )
    partition = {igG.vs[i]["name"]: part.membership[i] for i in range(igG.vcount())}

    return partition

#3.
def _hierarchical_clustering_analysis(
        data : pd.DataFrame,
        dendogram_ax : Axes
) :

    distance_matrix = _compute_jaccard_matrix(data)
    Z = linkage(distance_matrix, method='complete')
    dendogram_ax = _make_dendogram(Z,data,dendogram_ax)
    
    return dendogram_ax

def _compute_jaccard_matrix(
        data : pd.DataFrame
) :

    # Compute Jaccard distance matrix
    jaccard_dist = pdist(data.T, metric='jaccard')  # .T to compute distances between RNAs (columns)
    jaccard_matrix = squareform(jaccard_dist)  # Convert to square matrix

    return pd.DataFrame(jaccard_matrix, index=data.columns, columns=data.columns)

def _make_dendogram(
        Z : np.ndarray,
        data : pd.DataFrame,
        ax : Axes
        ) :

    dendrogram(Z, labels=data.columns, leaf_rotation=0, orientation="left", ax=ax)
    sns.despine(ax=ax, left=True, right=False)
    plt.title("Hierarchical Clustering")
    plt.xlabel("Jaccard distance")

    return ax

def expression_exploration(
        Spots : pd.DataFrame,
        run_path : str,
        control_genes : list[str] | None,
        na_policy : Literal["drop","fill"] = 'fill',
) -> dict['str', Figure]:
    """
    Perform exploratory analysis to identify structures in RNAs expression data.
    """

    data_path = os.path.join(run_path,"result_tables","expression.csv")
    if os.path.isfile(data_path) :
        print(f"Using cached data at {data_path}. If you modified user parameters delete cached data.")
        logging.info(f"Using cached data at  : {data_path}")
        data = pd.read_csv(data_path, sep=";").iloc[:,1:]
    else :
        data = _make_data_table_expression(Spots)
        data.to_csv(data_path,sep=";")

    if na_policy == "fill" : # Na values when a distribution was not detected in cell
        data = data.fillna(0) 
    elif na_policy == "drop" :
        ref_len = len(data)
        data = data.dropna()
        logging.info(f"{ref_len - len(data)} observation dropped according to 'drop' na_policy.")
    else :
        raise ValueError(f"na_policy should be either 'drop' or 'fill' got {na_policy}")

    res = {}

    #1.Correlation matrix.
    res["correlation_figure"] = expression_correlation_analysis(data)

    #2. Principal component analysis (No normalisation)
    res['principal_component_analysis_figure'] = principal_component_analysis(data, covariates=data.columns.to_list(), covariates_are_control_genes=False)

    #3. Principal compnent analysis normalised with control genes
    if not control_genes is None and not len(control_genes) == 0 :
        res['normalised_principal_component_analysis_figure'] = principal_component_analysis(data, covariates=control_genes, covariates_are_control_genes=True)



    return res



def _make_data_table_expression(
        Spots : pd.DataFrame,
) :
    COLUMNS = [
        "location",
        "spot_id",
        "cluster_id",
        "target",
        "z","y","x",
        "in_nucleus",
        "cell_id",
    ]

    data = Spots.loc[:,COLUMNS]
    data = data.groupby(["location","cell_id","target"])["spot_id"].count().rename("count").reset_index(drop=False)
    data["count"] = data["count"].astype(pd.Int32Dtype())
    data = data.pivot(
    index="cell_id",
    columns="target",
    values="count"
    )

    return data


#1.
def expression_correlation_analysis(data : pd.DataFrame) :

    correlation_matrix = data.corr()
    Z = linkage(correlation_matrix)
    dendogram_dict = dendrogram(Z, no_plot=True)
    hierarchical_order = dendogram_dict["leaves"]

    correlation_matrix = pd.DataFrame(correlation_matrix, columns= data.columns, index=data.columns)
    correlation_matrix = correlation_matrix.iloc[hierarchical_order,hierarchical_order]
    correlation_fig = _make_correlation_plot(correlation_matrix)

    return correlation_fig

def _make_correlation_plot(correlation_matrix : pd.DataFrame) :

    fig_size = max(10, len(correlation_matrix)*0.5)
    fig = plt.figure(figsize=(fig_size+2,fig_size)) #+2 for colorbar space
    ax = fig.gca()

    colors = [
    (0, "#46a7e7"),    # blue at -1
    (0.1, '#ffffff'),  # White at -0.7
    (0.3, "#87cf87"),    # Green at 0.0
    (0.5, "#BEC0BE"),    # Green at 0.0
    (0.7, "#87cf87"),    # Green at 0.0
    (0.9, '#ffffff'),  # White at 0.7
    (1, "#D30707")
    ]    # Red at 1.0

    # Create the colormap
    cmap = LinearSegmentedColormap.from_list('green_white_red', colors)

    sns.heatmap(
        correlation_matrix, 
        annot=True, 
        fmt=".2f", 
        cmap=cmap, 
        vmax=1, vmin=-1,
        linecolor="white",
        linewidths=.2,
    )

    ax.set_xticklabels(ax.get_xticklabels(), rotation=30)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=30)
    ax.set_ylabel('')
    ax.set_xlabel('')

    ax.set_title("Correlation Matrix of RNAs expression")

    return fig


#2.
def principal_component_analysis(
        data : pd.DataFrame,
        covariates : list[str],
        covariates_are_control_genes = False
        ) :

    if any(cov not in data.columns for cov in covariates) :
        raise KeyError(f"At least one covariate from {covariates} was not found in data.columns.")
    if not isinstance(covariates, list) : 
        raise TypeError(f"covariates arguments must be list, {type(covariates)} was passed.")

    #Init plot
    fig = plt.figure(figsize=(22,18))
    grid = gridspec.GridSpec(nrows=9,ncols=8, figure=fig, hspace=.8, wspace=.000001)
    det_metric_ax = fig.add_subplot(grid[0,0])
    shpericity_metric_ax = fig.add_subplot(grid[1,0])
    kmo_metric_ax = fig.add_subplot(grid[2,0])
    scree_ax = fig.add_subplot(grid[:3,2:6])
    pie_ax = fig.add_subplot(grid[:3,6:])
    loadings_ax = fig.add_subplot(grid[4:,:2])
    partial_correlation_ax = fig.add_subplot(grid[4:,2:])

    #partial correlation
    partial_correlation_matrix = _compute_partial_correlation(
        data=data,
        covariates=covariates,
        remove_covariates_from_interest_genes = covariates_are_control_genes
    )
    partial_correlation_ax = _make_partial_correlation_plot(partial_correlation_matrix, covariates, partial_correlation_ax)
    if covariates_are_control_genes :
        data = _normalise_data_with_covariates(data, covariates)

    #sanity metrics
    sanity_metrics = _sanity_metrics(data)
    det_metric_ax = _create_indicator_bar_plot(
        ax=det_metric_ax,
        metric_value=sanity_metrics["corr_det"],
        minvalue=1e-5, maxvalue=1,
        threshold_bad=0.01 , threshold_good=0.1,
        logscale=True
    )
    det_metric_ax.set_xlabel("Correlation matrix determinant")

    shpericity_metric_ax = _create_indicator_bar_plot(
        ax=shpericity_metric_ax,
        metric_value=sanity_metrics["corr_det"],
        minvalue=1e-5, maxvalue=1,
        threshold_bad=0.05, threshold_good=0.01,
        logscale=True
    )
    shpericity_metric_ax.set_xlabel("Bartlett sphericity test")

    kmo_metric_ax = _create_indicator_bar_plot(
        ax=kmo_metric_ax,
        metric_value=sanity_metrics["kmo"],
        minvalue=0, maxvalue=1,
        threshold_bad= .7, threshold_good=.8,
        logscale=False
    )
    kmo_metric_ax.set_xlabel("KMO test")

    



    #Principal component analysis core
    pca, standardised_data = _run_pca(data)
    pie_ax = _make_pca_pie_plot(pca, pie_ax)
    scree_ax = _make_eigenvalues_plot(pca,scree_ax)
    loadings_ax = _make_loadings_plot(pca,standardised_data,loadings_ax)

    return fig

def _normalise_data_with_covariates(
        data : pd.DataFrame,
        covariates : list[str]
) :
    control_median = data.loc[:,covariates].median(axis=1)

    interest_column = data.columns[~data.columns.isin(covariates)]
    normalised_data = data[interest_column].div(control_median,axis=0)

    return normalised_data
    

def _sanity_metrics(data : pd.DataFrame) :
    """
    **Keys** : 'kmo', 'corr_det', 'sphericity'
    """

    _,overall_kmo = calculate_kmo_manual(data)
    sphericity_test_pvalue = bartlett_sphericity_test(data)
    det = np.linalg.det(data.corr())

    return {
        "kmo" : overall_kmo,
        "corr_det" : det,
        "sphericity" : sphericity_test_pvalue
    }

def _compute_partial_correlation(
        data : pd.DataFrame,
        covariates : list[str],
        remove_covariates_from_interest_genes = False
) :

    if remove_covariates_from_interest_genes and data.columns.to_list() == covariates :
        raise ValueError("Cannot remove covariates from interest genes as they are identical")

    covariates_ser = pd.Series(covariates.copy())
    if remove_covariates_from_interest_genes :
        interest_columns = data.columns[~data.columns.isin(covariates_ser)].to_list()
        if len(interest_columns) < 3 : raise ValueError("Incorrect covariates, removal from interest genes yielded less than 3 interest genes.")
    else :
        interest_columns = data.columns.to_list()
    
    genes_combinations = list(combinations(interest_columns,r=2))
    partial_corr_matrix = pd.concat([pg.partial_corr(
        data,
        x=x,
        y=y,
        covar= list(covariates_ser[~(covariates_ser.str.contains(x) | covariates_ser.str.contains(y))])
    ) for x,y in genes_combinations], 
    ignore_index=True)

    partial_corr_matrix["x"] = list(zip(*genes_combinations))[0]
    partial_corr_matrix["y"] = list(zip(*genes_combinations))[1]
    partial_corr_matrix = partial_corr_matrix.pivot(columns="x",index="y",values="r")

    return partial_corr_matrix


def _run_pca(data : pd.DataFrame) :

    standardised_data = scale(data)
    standardised_data = pd.DataFrame(standardised_data, columns=data.columns, index = data.index)

    pca = PCA().fit(standardised_data)

    return pca, standardised_data

def _make_pca_pie_plot(
        pca : PCA, 
        ax : Axes
        ) :

    df = pd.DataFrame({
        "component" : ["PC"+str(i) for i in range(1, len(pca.explained_variance_ratio_)+1)],
        "var_ratio" : pca.explained_variance_ratio_,
        "eigenvalues" : pca.explained_variance_
    })

    df = df.loc[df["eigenvalues"] > 1] #Kaiser criterion
    df = pd.concat([
        df,
        pd.DataFrame({"component" : ["Noise"], "var_ratio" : [1- df["var_ratio"].sum()]})
    ])

    colors = cycle(['#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f','#bf5b17'])
    colors = [next(colors) for i in range(len(df))]
    colors[len(df)-1] = "#ACACAC"

    ax.pie(
        df["var_ratio"],
        colors=colors,
        explode=[0]*(len(df) -1) + [0.1],
        labels=df["var_ratio"].apply(lambda x : str(round(x*100)) + " %" ),
        hatch=[""]*(len(df) -1) + ["/"],
        labeldistance=0.3,
        textprops={"fontweight" : "bold", "size": 14}
        )

    ax.set_xlabel("Component contribution to variance", fontdict={"size" : 12})
    ax.legend(handles=(Patch(facecolor=color) for color in colors), labels=df["component"].to_list())

    return ax

def _make_eigenvalues_plot(
        pca : PCA,
        ax : Axes
) :

    eighenvalues = pca.explained_variance_
    eighenvalues = pd.DataFrame({
        0 : eighenvalues
    }).reset_index(drop=False)

    elbow_point = _find_elbow(eighenvalues)
    ax = _make_scree_plot(eighenvalues, ax)

    xmin,xmax = ax.get_xlim()
    ymin,ymax = ax.get_ylim()
    ax.plot([xmin,xmax],[1,1], "--", c="#fdbb84", label="kaiser criterion")
    ax.scatter(x=eighenvalues.iloc[elbow_point].at["index"], y=eighenvalues.iloc[elbow_point].at[0], 
               s=250, edgecolors="black",c="#fdbb84", label="elbow_point",
               zorder=10
               )

    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    ax.legend()
    ax.set_xlabel("Component")
    ax.set_ylabel("Eigenvalues")
    ax.set_title("Principal component analysis (PCA)", fontdict={"weight" : "bold", "size" : 15})
    
    return ax

def _make_loadings_plot(
        pca : PCA,
        standardised_data : pd.DataFrame,
        ax : Axes,
) :

    # Get loadings (eigenvectors)
    loadings = pd.DataFrame(
        pca.components_.T,
        index=standardised_data.columns,
        columns=[f"PC_{i+1}" for i in range(pca.n_components_)]
    )

    #eigenvalues
    eigenvalues = pca.explained_variance_

    #KAISER CRITERION
    loadings = loadings.loc[:,loadings.columns[eigenvalues >=1]]

    #Plot
    sns.heatmap(loadings, 
                annot=True, fmt=".2f", 
                cmap="coolwarm", vmax=1, vmin=-1, 
                linecolor="black", linewidths=0.5, square=True,
                ax=ax,
                cbar=False
                )

    ax.set_xticklabels(ax.get_xticklabels(), rotation=30)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=30)
    ax.set_ylabel('rna loadings')
    ax.set_xlabel('')

    return ax

def _make_partial_correlation_plot(
        partial_corr_matrix : pd.DataFrame,
        covariates : list[str],
        ax : Axes,
) :

    # Partial corr heatmap
    sns.heatmap(
        partial_corr_matrix, 
        annot=True, fmt=".2f", 
        cmap="coolwarm", vmax=1, vmin=-1, center=0, 
        linewidths=0.2, linecolor="white",
        ax=ax
        )

    ax = plt.gca()
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=30)
    ax.set_ylabel('')
    ax.set_xlabel('')

    if set(covariates) != set(partial_corr_matrix.index.to_list() + partial_corr_matrix.columns.to_list()) :
        ax.set_title("Partial correlation matrix\n" + f"Covariates {covariates}")
    else :
        ax.set_title("Partial correlation matrix\n")

    return ax

def calculate_kmo_manual(data : pd.DataFrame):
    #TODO This function can be optimized
    """
    Compute KMO (Kaiser-Meyer-Olkin) measure of sampling adequacy manually.
    Args:
        data: Pandas DataFrame (samples x variables).
    Returns:
        kmo_values: KMO for each variable (1D numpy array).
        overall_kmo: Overall KMO score (float).
    """
    n = len(data.columns)  # Number of variables
    corr_matrix = data.corr()
    kmo_values = np.zeros(n)
    ALL_GENES = corr_matrix.columns  # Define once outside the loop

    for i in range(n):
        test_gene = ALL_GENES[i]

        # Sum of squared correlations for variable i (excluding self-correlation)
        sum_r_squared = np.sum(corr_matrix.loc[test_gene, :].pow(2)) - np.pow(corr_matrix.loc[test_gene, test_gene], 2)

        # Sum of squared partial correlations for variable i
        sum_a_squared = 0
        for paired_gene in list(ALL_GENES[:i]) + list(ALL_GENES[i+1:]):
            covar = [g for g in ALL_GENES if g not in [test_gene, paired_gene]]
            result = pg.partial_corr(
                data=data,
                x=test_gene,
                y=paired_gene,
                covar=covar
            )
            a_ij = result['r'].values[0]
            sum_a_squared += a_ij ** 2

        # Compute KMO for variable i
        kmo_i = sum_r_squared / (sum_r_squared + sum_a_squared)
        kmo_values[i] = kmo_i

    overall_kmo = np.mean(kmo_values)
    return kmo_values, overall_kmo


def bartlett_sphericity_test(data):
    """
    Perform Bartlett's test of sphericity for PCA suitability.
    Args:
        data: pd.DataFrame.
    Returns:
        p_value: p-value for the test
    """

    n = len(data)  # Number of samples
    p = len(data.columns)  # Number of variables

    # Compute the correlation matrix
    corr_matrix = data.corr()

    # Compute the determinant of the correlation matrix
    det_r = np.linalg.det(corr_matrix)

    # Compute Bartlett's test statistic
    chi_square = - (n - 1 - (2 * p + 5) / 6) * np.log(det_r)

    # Degrees of freedom
    df = p * (p - 1) / 2

    # Compute p-value from chi-square distribution
    p_value = 1 - chi2.cdf(chi_square, df)

    return p_value

def _create_indicator_bar_plot(
    ax : Axes,
    metric_value : int|float,
    minvalue : int|float,
    maxvalue : int|float,
    threshold_good : int|float,
    threshold_bad : int|float,
    logscale = False
) :

    if maxvalue > threshold_good  >= threshold_bad > minvalue :
        ascending = True
    elif minvalue < threshold_good <= threshold_bad < maxvalue :
        ascending = False
    else :
        raise ValueError("Arguments must follow one of following rules : " \
        "\nmaxvalue > threshold_good > threshold_medium > threshold_bad > minvalue\n"
        "minvalue < threshold_good < threshold_medium < threshold_bad < maxvalue")


    if metric_value <= minvalue : 
        if metric_value == 0 :
            metric_value = 0.1
        else :
            metric_value = minvalue*1.2

    if metric_value >= maxvalue :
        metric_value = maxvalue

    ax.barh(
        left=threshold_good if ascending else minvalue,
        y=0.5,
        width=(maxvalue-threshold_good) if ascending else (threshold_good-minvalue),
        height=1,
        facecolor="#31a354",
        edgecolor="black",
        alpha=.7,
        align="center",
    )
    ax.barh(
        left=threshold_bad if ascending else threshold_good,
        y=0.5,
        width=(threshold_good-threshold_bad) if ascending else (threshold_bad - threshold_good),
        height=1,
        facecolor="#fefb4c",
        edgecolor="black",
        alpha=.7,
        align="center",
    )

    ax.barh(
        left=minvalue if ascending else threshold_bad,
        y=0.5,
        width=(threshold_bad-minvalue) if ascending else (maxvalue - threshold_bad),
        height=1,
        facecolor="#de2d26",
        edgecolor="black",
        alpha=1,
        align="center",
    )

    ax.scatter(
        x=metric_value,
        y=0.5,
        s=100,
        facecolor="white",
        edgecolor="black",
        alpha=1,
    )

    xlim=[minvalue,maxvalue+0.2]
    ylim = [-1,3]

    if logscale: 
        ax.set_xscale("log")
        if minvalue <= 0 : raise ValueError("Cannot set minvalue <= 0 with logscale = True.")

    xticks = [minvalue, threshold_bad, threshold_good, maxvalue] if ascending else [minvalue, threshold_good, threshold_bad, maxvalue]

    ax.set_xlim(*xlim)
    ax.set_xticks(xticks, labels=xticks, rotation=30, ha="right")
    ax.set_ylim(*ylim)
    ax.axes.set_yticks([])

    sns.despine(ax=ax, bottom=True)
    ax.axes.yaxis.set_visible(False)

    return ax