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

from Sequential_Fish.analysis.dashboards.signal_quality import compute_similarity_metrics 
from Sequential_Fish.tools import reorder_image_stack

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
        "spatial_correlation" : da.map_blocks(
            lambda a, b: pearsonr(a.flatten(), b.flatten())[0],
            image1, 
            image2,
            dtype=float,
            ),
        "structural_similarity" : da.map_blocks(
            lambda a,b : ssim(a,b,full=True)[0],
            image1, 
            image2, 
            dtype=float
            ),
        "mean_squared_error" : mean_squared_error(image1, image2),
    }
    res["peak_snr"] = 10 * log10(max_intensity**2 / res["mean_squared_error"])

    return res

@profile
def main() :
    RUN_PATH = "/home/floric/Documents/SeqfFish/2026-03-24 - HeLa_POLR2_Run17_tiff/"

    Acquisition = pd.read_feather(RUN_PATH + "/result_tables/Acquisition.feather")
    Acquisition["full_path"] = Acquisition["full_path"].str.replace("/media/SSD_floricslimani/Fish_seq/POLR2/2026-03-24 - HeLa_POLR2_Run17_tiff/",RUN_PATH)
    shape = tuple(Acquisition["fish_reodered_shape"].at[0])


    print("opening location")
    loc1_stack = lazy_open_location(Acquisition,"Location-01")
    loc1_stack = da.from_delayed(loc1_stack, shape=(16,)+ shape, dtype=int)

    print("selecting images")
    print(loc1_stack.shape)
    image1 = loc1_stack[0,...,0]
    image2 = loc1_stack[-4,...,0]

    print("image1 : ",image1.shape)
    print("image2 : ",image2.shape)

    print("computing similarity metrics")
    res = compute_similarity_metrics_dask(
        image1 = image1,
        image2 = image2
        )

    clock = t.time()
    print(res['spatial_correlation'].compute())
    print(f"spatial_correlation : {t.time()- clock}")
    clock = t.time()
    print(res['structural_similarity'].compute())
    print(f"structural_similarity : {t.time()- clock}")
    clock = t.time()
    print(res['peak_snr'])
    print(f"peak_snr : {t.time()- clock}")

    print("done.")
if __name__ == "__main__" :
    main()