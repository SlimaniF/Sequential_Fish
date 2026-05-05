import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram

from matplotlib.colors import LinearSegmentedColormap
import numpy as np

# Define the colors and their positions
colors = [
    (0.0, '#4575b4'),    # Blue at -1
    (0.15, '#91bfdb'),   # Light blue at -0.7
    (0.35, '#e0f3f8'),   # Gray at -0.3
    (0.5, '#ffffff'),    # White at 0
    (0.65, '#fee090'),   # Yellow at 0.3
    (0.85, '#fc8d59'),   # Orange at 0.7
    (1.0, '#d73027')     # Red at 1
]

cmap = LinearSegmentedColormap.from_list('divergent_palette', colors, N=256)

if __name__ == "__main__" :

    run_path = "/media/SSD_floricslimani/Fish_seq/POLR2/2026-03-24 - HeLa_POLR2_Run17_tiff/result_tables/"
    data = pd.read_feather(run_path + "signal_correlations.feather")

    linkage_matrix = linkage(data, optimal_ordering=True)
    dn = dendrogram(linkage_matrix.astype(float), no_plot=True)
    order = np.array(dn["leaves"],dtype=int)
    data = data.iloc[order,order]
    data.loc[:,:] = np.triu(data)

    fig = plt.figure()
    ax = fig.gca()

    sns.heatmap(
        data.T,
        vmin=-1,
        vmax=1,
        cmap=cmap,
        ax=ax,
        square=True,
        linecolor="white",
        linewidths=1
    )

    sns.despine(ax=ax)

    ax.set_xticklabels(ax.get_xticklabels(),rotation=30, ha="right")
    ax.set_title("Signal spatial correlation", size=15, weight="bold")

    plt.show()