import os
import tifffile

import numpy as np
from skimage.metrics import structural_similarity as ssim
from skimage.metrics import mean_squared_error
from scipy.stats import pearsonr
from math import log10

import pandas as pd
import dask.delayed
import dask.array as da
from memory_profiler import profile
from pebble import ProcessPool

from dask.distributed import Client
from dask.diagnostics import ProgressBar

from itertools import cycle

from Sequential_Fish.analysis.dashboards.signal_quality import compute_similarity_metrics 
from Sequential_Fish.tools import reorder_image_stack
from tqdm import tqdm

import time as t

def open_location(
        Acquisition : pd.DataFrame,
        location : str,
) :
    """
    Open all cycles of a location and reorder stacks in order (cycle,z,y,x,channel)
    """

    if not ('location' in Acquisition.index.names and 'cycle' in Acquisition.index.names) :
        Acquisition = Acquisition.set_index(['location','cycle'], verify_integrity=True)
    
    fish_path = Acquisition.at[(location,0), 'full_path']
    
    with tifffile.TiffFile(fish_path) as tif :
        location_stack = tif.asarray()

    stack_map = Acquisition.at[(location,0), 'fish_map']
    location_stack = reorder_image_stack(location_stack, channel_map=stack_map)

    return location_stack

@dask.delayed
def lazy_load_tiff(tiff_path, stack_map) :
    with tifffile.TiffFile(tiff_path) as tif :
        location_stack = tif.asarray()
    location_stack = reorder_image_stack(location_stack, channel_map=stack_map)

    return location_stack

def lazy_open_location(
        Acquisition : pd.DataFrame,
        location : str,
) :
    """
    Open all cycles of a location and reorder stacks in order (cycle,z,y,x,channel)
    """

    if not ('location' in Acquisition.index.names and 'cycle' in Acquisition.index.names) :
        Acquisition = Acquisition.set_index(['location','cycle'], verify_integrity=True)
    
    fish_path = Acquisition.at[(location,0), 'full_path']
    
    stack_map = Acquisition.at[(location,0), 'fish_map']
    location_stack= lazy_load_tiff(fish_path, stack_map)

    return location_stack


def compute_similarity_metrics_dask(
    image1 : np.ndarray,
    image2 : np.ndarray
    ) :

    if np.issubdtype(image1.dtype, np.floating) :
        max_intensity = 1
    else :
        max_intensity = np.iinfo(image1.dtype).max

    res = {
        "spatial_correlation" : pearsonr(image1.flatten(), image2.flatten())[0],
        "structural_similarity" : ssim(image1,image2,full=True, data_range= max_intensity)[0],
        "mean_squared_error" : mean_squared_error(image1, image2),
    }
    res["peak_snr"] = 10 * log10(max_intensity**2 / res["mean_squared_error"])

    return res


@profile
def main() :
    RUN_PATH = "/media/SSD_floricslimani/Fish_seq/POLR2/2026-03-24 - HeLa_POLR2_Run17_tiff/"

    Acquisition = pd.read_feather(RUN_PATH + "/result_tables/Acquisition.feather")
    Acquisition["full_path"] = Acquisition["full_path"].str.replace("/media/SSD_floricslimani/Fish_seq/POLR2/2026-03-24 - HeLa_POLR2_Run17_tiff/",RUN_PATH)
    shape = tuple(Acquisition["fish_reodered_shape"].at[0])

    init_clock = t.time()
    list_dir = os.listdir(RUN_PATH + '/FISH_Z-stacks/')
    measures = []

    for location in tqdm(list_dir[:5]) :
        loc1_stack = open_location(Acquisition,location)
        # loc1_stack = da.from_delayed(loc1_stack, shape=(16,)+ shape, dtype=int)

        image1 = loc1_stack[0,...,0]
        image2 = loc1_stack[-4,...,0]
        del loc1_stack


        measures.append(compute_similarity_metrics_dask(
            image1 = image1,
            image2 = image2
            ))


    print(measures)
    print(f"done; Execution time : {t.time()-init_clock}")


# def process_acquisition(
#     Acquisition,
#     im_shape,
#     location,
# ) : 
#     print(f"computing : {location}")
#     loc1_stack = lazy_open_location(Acquisition,location)
#     loc1_stack = da.from_delayed(loc1_stack, shape=(16,)+ im_shape, dtype=int)

#     image1 = loc1_stack[0,...,0]
#     image2 = loc1_stack[-4,...,0]


#     res = compute_similarity_metrics_dask(
#         image1 = image1,
#         image2 = image2
#         )

#     print(f"returning : {location}")
#     return res

# @profile
# def main() :
#     RUN_PATH = "/media/SSD_floricslimani/Fish_seq/POLR2/2026-03-24 - HeLa_POLR2_Run17_tiff/"

#     Acquisition = pd.read_feather(RUN_PATH + "/result_tables/Acquisition.feather")
#     Acquisition["full_path"] = Acquisition["full_path"].str.replace("/media/SSD_floricslimani/Fish_seq/POLR2/2026-03-24 - HeLa_POLR2_Run17_tiff/",RUN_PATH)
#     shape = tuple(Acquisition["fish_reodered_shape"].at[0])

#     init_clock = t.time()
#     list_dir = os.listdir(RUN_PATH + '/FISH_Z-stacks/')[:5]

#     client = Client(n_workers=4)
#     futures = client.map(
#         process_acquisition,
#         [Acquisition]*len(list_dir),
#         [shape]*len(list_dir),
#         list_dir
#     )

#     measures = client.gather(futures)
#     client.close()
#     print(measures)

#     print(f"done; Execution time : {t.time()-init_clock}")


if __name__ == "__main__" :
    main()