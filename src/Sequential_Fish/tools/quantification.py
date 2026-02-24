import warnings
import numpy as np
import pandas as pd
from skimage.measure import regionprops_table
import bigfish.multistack as multistack
import bigfish.classification as classification
import bigfish.stack as stack

def cell_quantification(
        acquisition_id,
        detection_id,
        spots,
        clusters,
        voxel_size,
        cell_label,
        nucleus_label,
        fish_signal,
        dapi_signal_3D,
        ndim=3,
        talk=False,
        
) :
    
    if len(spots) == 0 :
        spots = np.empty(shape=(0,3), dtype=int)
    if len(clusters) == 0 :
        clusters = np.empty(shape=(0,5),dtype=int)

    if fish_signal.ndim == 3 :
        fish_signal = np.max(fish_signal, axis=0)

    cell_extraction = multistack.extract_cell(
            cell_label=cell_label,
            others_image= {"smfish" : fish_signal},
            ndim=ndim,
            nuc_label=nucleus_label,
            rna_coord= spots,
            others_coord= {'foci' : clusters}
        )
    
    cell_df = pd.DataFrame()
    for cell in cell_extraction :
        new_cell = _main_cell_quant(
            cell=cell,
            voxel_size=voxel_size,
            dapi_stack=dapi_signal_3D,
        )

        new_cell['acquisition_id'] = acquisition_id
        new_cell['detection_id'] = detection_id
        new_cell['cluster_number'] = len(cell['foci'])
        cell_spots = cell['rna_coord']
        if len(cell_spots) > 0 : cell_spots = _keep_spots_in_mask(cell_spots, cell['cell_mask'])
        new_cell['rna_number'] = len(cell_spots)

        cell_df = pd.concat([
            cell_df,
            new_cell
        ],axis=0)

    if talk : print("detection {0} done.".format(detection_id))
    return cell_df


def _main_cell_quant(cell, voxel_size, dapi_stack=None, compute_centrosome= False, centrosome_coords= None) :
    """
    Basic function for cell quantification using bigfish AND COMPUTING DAPI MEASUREMENTS if dapi_stack != None.

    Parameters
    ----------
        cell : dict  
            computed from `bigfish.multistack.extract_cell`
        voxel_size : tuple(z,y,x)
        dapi_stack : np.ndarray
            3D uncropped dapi stack.

    Returns
    -------
        Cell : pd.DataFrame  
    """
    
    #Extracting bigfish cell information
    if isinstance(voxel_size, (list, np.ndarray)) : voxel_size = tuple(voxel_size)
    voxel_size_yx = float(voxel_size[1])
    cell_mask: np.ndarray = cell["cell_mask"]
    nuc_mask = cell["nuc_mask"] 
    rna_coord = cell["rna_coord"]
    smfish = cell["smfish"]
    min_y, min_x, max_y, max_x = cell["bbox"]
    label = cell["cell_id"] # is the label of this cell in cell_label

    if "transcription_site" in cell.keys() :
        ts_coord = cell["transcription_site"]
    else : ts_coord = np.empty(shape=(0,0), dtype= np.int64)

    if "foci" in cell.keys() :
        foci_coord = cell["foci"]
    else : foci_coord = np.empty(shape=(0,0), dtype= np.int64)

    if not isinstance(centrosome_coords, (np.ndarray, list)) and compute_centrosome:
        print("Warning : compute_centrosome is set to True but centrosome_coords parameter is not np.ndarray nor list. \nCompute centrosome has been set to False; centrosome_coords type : {0}".format(type(centrosome_coords)))
        compute_centrosome = False
    if type(centrosome_coords) == list : centrosome_coords = np.array(centrosome_coords, dtype=int)
    
    warnings.filterwarnings("ignore")

    features, features_names = classification.compute_features(
        cell_mask= cell_mask, 
        nuc_mask= nuc_mask, 
        ndim= 3, 
        rna_coord= rna_coord, 
        smfish= smfish, 
        foci_coord= foci_coord, 
        voxel_size_yx= voxel_size_yx, 
        centrosome_coord= centrosome_coords,
        compute_centrosome=compute_centrosome,
        compute_distance=True,
        compute_intranuclear=True,
        compute_protrusion=True,
        compute_dispersion=True,
        compute_topography=True,
        compute_foci=True,
        compute_area=True,
        return_names=True
    )
    warnings.filterwarnings("default")                                                      

    #Custom features
    cell_props_table = regionprops_table(cell_mask.astype(int), properties= ["centroid"])
    cell_coordinates = (float(cell_props_table["centroid-0"] + min_y), float(cell_props_table["centroid-1"] + min_x))
    del cell_props_table
    nucleus_area_px = compute_mask_area(nuc_mask, unit= 'px', voxel_size= voxel_size)
    nucleus_area_nm = compute_mask_area(nuc_mask, unit= 'nm', voxel_size= voxel_size)
 
    #Adding custom features to DataFrame
    features = list(features)
    features.extend([cell_coordinates, label, cell["bbox"], nucleus_area_px, nucleus_area_nm,])
     
    features_names += [ "cell_coordinates", "label", "bbox", "nucleus_area_px", "nucleus_area_nm"]
    
    #signal features
    if type(dapi_stack) != type(None) :
        nucleus_mip_signal_metrics = nucleus_signal_metrics(cell, channel= dapi_stack, projtype= 'mip')
        nucleus_mean_signal_metrics = nucleus_signal_metrics(cell, channel= dapi_stack, projtype= 'mean')
        features.extend([nucleus_mip_signal_metrics["mean"], nucleus_mip_signal_metrics["max"], nucleus_mip_signal_metrics["min"], nucleus_mip_signal_metrics["median"],
                         nucleus_mean_signal_metrics["mean"], nucleus_mean_signal_metrics["max"], nucleus_mean_signal_metrics["min"], nucleus_mean_signal_metrics["median"]])
     
        features_names += ["nucleus_mip_mean_signal","nucleus_mip_max_signal","nucleus_mip_min_signal","nucleus_mip_median_signal",
                           "nucleus_mean_mean_signal","nucleus_mean_max_signal","nucleus_mean_min_signal","nucleus_mean_median_signal"]



    data = features
    header = features_names

    new_cell = pd.DataFrame(data= [data], columns= header)
    return new_cell

def _keep_spots_in_mask(spots,mask:np.ndarray) :
    spot_zip = zip(*spots)
    spots_indexer = tuple(list(spot_zip))
    if mask.ndim == 2 and len(spots[0]) == 3 :
        spots_indexer = spots_indexer[1:]

    spots = spots[mask[spots_indexer]]

    return spots


def compute_mask_area(mask: np.ndarray, unit: str = 'px', voxel_size: tuple= None)-> float:
    """
    Return the area of pbody within cell. 
    
    Parameters
    ----------
        mask: np.ndarray
            mask computed for current cell.
        unit : str
            Unit parameter can be set to either 'px' or 'nm'. If nm is selected voxel_size (z,y,x) or (y,x) has to be given as well.
        voxel_size : tuple
    """
    #Mask must be boolean
    if mask.dtype != bool : raise TypeError("'mask' parameter should be a ndarray with boolean dtype.")

    #Set pixel or nanometers
    if unit.upper() in ["PX", "PIXEL"] : return_pixel = True
    elif unit.upper() in ["NM", "NANOMETER"] and type(voxel_size) in [tuple, list] : return_pixel = False
    elif unit.upper() in ["NM", "NANOMETER"] and voxel_size == None : raise TypeError("When scale is set to nanometer, voxel_size has to be given as a tuple or a list.")
    else : raise ValueError("unit parameter incorrect should be either 'px' or 'nm'. {0} was given".format(unit))
    
    if mask.ndim != 2 : raise ValueError("Only 2D masks are supported")

    if not return_pixel :
        if len(voxel_size) == 2 :
            y_dim = voxel_size[0]
            x_dim = voxel_size[1]
        elif len(voxel_size) == 3 :
            y_dim = voxel_size[1]
            x_dim = voxel_size[2]
        else : raise ValueError("Inapropriate voxel_size length. Should be either 2 or 3, it is {0}".format(len(voxel_size)))

    pixel_number = mask.sum()

    if return_pixel : res = pixel_number
    else : res = pixel_number * y_dim * x_dim

    return res


def nucleus_signal_metrics(cell, channel, projtype = 'mip', use_cell_mask= False) :
    """
    Returns dict containing signal related measures : 'min', 'max', '1 percentile', '9 percentile', 'mean' and 'median'.
      Computed from channel signal in cell's nucleus mask. Signal measures are computed from 2D cell, so channel is projected z-wise according to projtype (provided channel is 3D).
    
        Parameters
        ----------
            cell : dict
                Dictionary computed from bigFish.multistack.extract_cell
            channel : np.ndarray
                Channel from which intensity is computed
            projtype : str
                can either be 'mip' or 'mean'.

        Returns
        -------
            mean_sig : float
        
    """
    min_y, min_x, max_y, max_x = cell["bbox"]
    channel_cropped = channel[..., min_y:max_y, min_x:max_x]
 

    if channel.ndim == 3 :
        if projtype == 'mip' : 
            channel_cropped = stack.maximum_projection(channel_cropped)
        elif projtype == 'mean' :
            channel_cropped = stack.mean_projection(channel_cropped)
    
    if use_cell_mask : nucleus_mask = cell["cell_mask"]
    else : nucleus_mask = cell["nuc_mask"]

    metrics = compute_signalmetrics(channel_cropped, nucleus_mask)
    return metrics

def compute_signalmetrics(signal:np.ndarray, mask: np.ndarray) :
    """Compute 'min', 'max', '1 percentile', '9 percentile', 'mean', 'std' and 'median' value from signal ignoring pixels not in mask.
        
        Parameters
        ----------
            signal : np.ndarray
            mask : np.ndarray

        Returns
        -------
            signalmetrics : dict
               'min', 'max', '1 percentile', '99 percentile', 'mean' and 'median'
    """
    if mask.dtype != bool : raise TypeError("'Mask' parameter should be a ndarray with dtype = bool")

    flat:np.ndarray = signal[mask]
    signalmetrics = {
        "min" : flat.min(),
        "max" : flat.max(),
        "1 percentile" : np.percentile(flat, 1),
        "99 percentile" : np.percentile(flat, 99),
        "mean" : flat.mean(),
        "std" : flat.std(),
        "median" : np.median(flat),
        "sum" : np.sum(flat) 
    }
    return signalmetrics