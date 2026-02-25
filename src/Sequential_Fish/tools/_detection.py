"""
Integrate custom functions to Sequential Fish for multi threaded spot detection using bigfish.
"""

import signal
import warnings
import numpy as np
import pandas as pd
import bigfish.detection as detection
import bigfish.stack as stack
import bigfish.plot as plot

from numpy import nan

from .utils import get_centroids_array, compute_anisotropy_coef, get_centroids_list, _compute_critical_spot_number

def multi_thread_full_detection(
        image,
        voxel_size,
        threshold,
        spot_radius,
        alpha,
        beta,
        gamma,
        artifact_size,
        cluster_radius,
        min_spot_per_cluster,
        filename:str,
        detection_id,
) :

    warning_msgs = []
    try :
        with warnings.catch_warnings(record=True) as wlist :
            res = full_detection(image,voxel_size,threshold,spot_radius,alpha,beta,gamma,artifact_size,cluster_radius,min_spot_per_cluster,filename,detection_id,)
            warning_msgs = [str(warning.message) for warning in wlist]

    except Exception as error :
        warning_msgs.append("Thread {0} : Error occured\n{1}".format(detection_id, str(error)))
        keys = ['spots','spots_post_decomp','clustered_spots_dataframe','clusters_dataframe','clusters','clustered_spots','free_spots','threshold','voxel_size','spot_radius','alpha','beta','gamma','artifact_size','cluster_radius','min_spot_per_cluster',]
        res = dict.fromkeys(keys, nan)
        res['detection_id'] = detection_id
        warning_msgs.append("Thread {0} : Returning nan values".format(detection_id))

    finally :
        return res, warning_msgs


def full_detection (
        image,
        voxel_size,
        threshold,
        spot_radius,
        alpha,
        beta,
        gamma,
        artifact_size,
        cluster_radius,
        min_spot_per_cluster,
        filename:str,
        detection_id,
) :
    """
    Performs full detection process : thread-compatible.

    1. detect spot (bigfish.detection)
    2. cluster_deconvolution (smfishtools.detection)
    3. remove_artificats (smfishtools.detection) (skipped if artificat_size is None)
    4. cluster_detection
    5. plot detection visual (skipped if filename is None)

    RETURN
    ------
    res : dict
        {
        'id'
        'spots'
        'spots_post_decomp'
        'clustered_spots_dataframe',
        'clusters_dataframe'
        'clusters'
        'clustered_spots'
        'free_spots'
        }


    """
    spots, threshold = detection.detect_spots(
            images= image,
            threshold=threshold,
            return_threshold = True,
            voxel_size=voxel_size,
            spot_radius=spot_radius,
        )

    spots_post_decomp = cluster_deconvolution(
        image=image,
        spots=spots,
        spot_radius=spot_radius,
        voxel_size=voxel_size,
        alpha=alpha,
        beta= beta,
        sigma=gamma,
    )

    if type(artifact_size) != type(None) :
        spots_post_decomp = remove_artifact(
            spots_post_decomp, 
            artifact_radius=artifact_size, 
            voxel_size=voxel_size
            )

    spots_cluster_results= cluster_detection(
        spots_post_decomp,
        voxel_size= voxel_size, 
        nb_min_spots=min_spot_per_cluster,
        radius=cluster_radius, 
        keys_to_compute=['clustered_spots_dataframe', 'clusters_dataframe', 'clusters']
        )
        
    clustered_spots_dataframe, clusters_dataframe, clusters = spots_cluster_results['clustered_spots_dataframe'], spots_cluster_results['clusters_dataframe'], spots_cluster_results['clusters']
    clustered_spots = get_centroids_array(clustered_spots_dataframe.query("not cluster_id.isna()"))
    free_spots = get_centroids_array(clustered_spots_dataframe.query("cluster_id.isna()"))
    if type(filename) != type(None) :
        plot.output_spot_tiffvisual(np.max(image, axis=0), [spots, spots_post_decomp, clustered_spots, free_spots], filename, dot_size= 2)
        filename = filename.split('/')[-1]

    clustered_spots_dataframe['intensity'] = image[
        list(clustered_spots_dataframe['z']),
        list(clustered_spots_dataframe['y']),
        list(clustered_spots_dataframe['x']),
    ]

    res = {
        'detection_id' : detection_id,
        'spots' : spots,
        'spots_post_decomp' : spots_post_decomp,
        'clustered_spots_dataframe' : clustered_spots_dataframe,
        'clusters_dataframe' : clusters_dataframe,
        'clusters' : clusters,
        'clustered_spots' : clustered_spots,
        'free_spots' : free_spots,
    }
    
    return res

def cluster_deconvolution(image, spots, spot_radius, voxel_size, alpha, beta, sigma=5, timer= 0) :
    """
    Wrapper handling time out during deconvolution and preprocessing with gaussian background removal --> sigma defines the kernel size if 0 then no denoising is performed.
    --> `smfishtools.detection.spot_decomposition_nobckgrndrmv`
    """

    if len(spots) == 0 : return spots

    if timer > 0 : signal.signal(signal.SIGALRM, detectiontimeout_handler) #Initiating timeout handling
    try :
        if sigma > 0 : im = remove_mean_gaussian_background(image, sigma=sigma, voxel_size=voxel_size)
        else : im = image
        if timer > 0 : signal.alarm(timer)
        spots_postdecomp = spot_decomposition_nobckgrndrmv(im, spots, spot_radius, voxel_size_nm=voxel_size, alpha= alpha, beta= beta)
    except DetectionTimeOutError :
        warnings.warn(" \033[91mCluster deconvolution timeout...\033[0m")
        spots_postdecomp = np.empty((0,0), dtype=int)
    except NoSpotError :
        warnings.warn(" No dense regions to deconvolute.")
        spots_postdecomp = spots
    except ValueError as e :
        if 'x0' in str(e) :
            warnings.warn('x0 is infeasible error raised during cluster deconvolution. (Gaussian fit error)')
            spots_postdecomp = spots
        else :
            raise(e)
    except RuntimeError as e:
        warnings.warn("Run time error {0}".format(e))
        spots_postdecomp = spots
    except Exception as error :
        raise error
    finally :
        if timer > 0 : signal.alarm(0)
    return spots_postdecomp


def cluster_detection(spots, voxel_size, radius = 350, nb_min_spots = 4, keys_to_compute = ["clustered_spots", "clusters"]) :
    """
    Performs `bigfish.detection.cluster_detection()` to detect clusters.
    Then offers possibility to get results sorted in pandas dataframe.

    Parameters
    ----------
        spots : np.ndarray
            Coordinates of the detected spots with shape (nb_spots, 3) or (nb_spots, 2).
        voxel_size : int, float, Tuple(int, float) or List(int, float)
            Size of a voxel, in nanometer. One value per spatial dimension (zyx or yx dimensions). If it's a scalar, the same value is applied to every dimensions.
        radius : int
            The maximum distance between two samples for one to be considered as in the neighborhood of the other. Radius expressed in nanometer.
        nb_min_spots : int
            The number of spots in a neighborhood for a point to be considered as a core point (from which a cluster is expanded). This includes the point itself.
        keys_to_compute : list[str], str
            keys from (clustered_spots, clusters, clustered_spots_dataframe, clusters_dataframe)
                --> clustered_spots : np.ndarray
                    Coordinates of the detected spots with shape (nb_spots, 4) or (nb_spots, 3). One coordinate per dimension (zyx or yx coordinates) plus the index of the cluster assigned to the spot. If no cluster was assigned, value is -1.
                --> clusters : np.ndarray
                    Array with shape (nb_clusters, 5) or (nb_clusters, 4). One coordinate per dimension for the clusters centroid (zyx or yx coordinates), the number of spots detected in the clusters and its index.
                --> clustered_spots_dataframe
                --> clusters_dataframe
    
    Returns
    -------
        res : dict
            keys : keys from `keys_to_compute` argument : (clustered_spots, clusters, clustered_spots_dataframe, clusters_dataframe)    
    """

    if isinstance(keys_to_compute, str) : keys_to_compute = [keys_to_compute]
    elif isinstance(keys_to_compute, list) : pass
    else : raise TypeError("Wrong type for keys_to_compute. Should be list[str] or str. It is {0}".format(type(keys_to_compute)))
    if len(spots) == 0 :
        res = {'clustered_spots' : [], 'clusters' : [], 'clustered_spots_dataframe' : pd.DataFrame(columns= ["id", "cluster_id", "z", "y", "x"]), 'clusters_dataframe' : pd.DataFrame(columns= ["id", "z", "y", "x", "spot_number"])}
        return {key : res[key] for key in keys_to_compute}
    else : res = {}
    clustered_spots, clusters = detection.detect_clusters(spots, voxel_size= voxel_size, radius= radius, nb_min_spots= nb_min_spots)


    if 'clustered_spots' in keys_to_compute :
        res['clustered_spots'] = clustered_spots
        
    if 'clusters' in keys_to_compute : 
        res['clusters'] = clusters

    if 'clustered_spots_dataframe' in keys_to_compute :
        res['clustered_spots_dataframe'] = _compute_clustered_spots_dataframe(clustered_spots)
    
    if 'clusters_dataframe' in keys_to_compute :
        res['clusters_dataframe'] = _compute_cluster_dataframe(clusters)

    return res


def spot_decomposition_nobckgrndrmv(image, spots, spot_radius, voxel_size_nm, alpha= 0.5, beta= 1):
    """ Basically same function as bigfish.detection.decompose_dense but without the remove background gaussian.
    
    Detect dense and bright regions with potential clustered spots and
    simulate a more realistic number of spots in these regions.

    #. We build a reference spot by aggregating predetected spots.
    #. We fit gaussian parameters on the reference spots.
    #. We detect dense regions to decompose.
    #. We simulate as many gaussians as possible in the candidate regions.

    Parameters
    ----------
    image : np.ndarray
        Image with shape (z, y, x) or (y, x).
    spots : np.ndarray
        Coordinate of the spots with shape (nb_spots, 3) or (nb_spots, 2)
        for 3-d or 2-d images respectively.
    voxel_size : int, float, Tuple(int, float) or List(int, float)
        Size of a voxel, in nanometer. One value per spatial dimension (zyx or
        yx dimensions). If it's a scalar, the same value is applied to every
        dimensions.
    spot_radius : int, float, Tuple(int, float) or List(int, float)
        Radius of the spot, in nanometer. One value per spatial dimension (zyx
        or yx dimensions). If it's a scalar, the same radius is applied to
        every dimensions.
    kernel_size : int, float, Tuple(float, int), List(float, int) or None
        Standard deviation used for the gaussian kernel (one for each
        dimension), in pixel. If it's a scalar, the same standard deviation is
        applied to every dimensions. If None, we estimate the kernel size from
        'spot_radius', 'voxel_size' and 'gamma'
    alpha : int or float
        Intensity percentile used to compute the reference spot, between 0
        and 1. The higher, the brighter are the spots simulated in the dense
        regions. Consequently, a high intensity score reduces the number of
        spots added. Default is 0.5, meaning the reference spot considered is
        the median spot.
    beta : int or float
        Multiplicative factor for the intensity threshold of a dense region.
        Default is 1. Threshold is computed with the formula:

        .. math::
            \\mbox{threshold} = \\beta * \\mbox{max(median spot)}

        With :math:`\\mbox{median spot}` the median value of all detected spot
        signals.
    gamma : int or float
        Multiplicative factor use to compute the gaussian kernel size:

        .. math::
            \\mbox{kernel size} = \\frac{\\gamma * \\mbox{spot radius}}{\\mbox{
            voxel size}}

        We perform a large gaussian filter with such scale to estimate image
        background and remove it from original image. A large gamma increases
        the scale of the gaussian filter and smooth the estimated background.
        To decompose very large bright areas, a larger gamma should be set.

    Notes
    -----
    If ``gamma = 0`` and ``kernel_size = None``, image is not denoised.

    Returns
    -------
    spots_postdecom : np.ndarray
        Coordinate of the spots detected, with shape (nb_spots, 3) or
        (nb_spots, 2). One coordinate per dimension (zyx or yx coordinates).

    """
    ndim = image.ndim


    #Spot decomposition    
    if len(spots) == 0 : return spots
    elif len(spots[0]) == 0 : return spots
    # reference spot
    reference_spot = detection.build_reference_spot(
        image=image,
        spots=spots,
        voxel_size=voxel_size_nm, 
        spot_radius= spot_radius,
        alpha=alpha)

    # fit a gaussian function on the reference spot
    if ndim == 3 :
        sigma_z, sigma_yx, amplitude, background = detection.modelize_spot(
            reference_spot=reference_spot, 
            voxel_size= voxel_size_nm, 
            spot_radius= spot_radius)
    else :
        sigma_yx, amplitude, background = detection.modelize_spot(
            reference_spot=reference_spot, 
            voxel_size= voxel_size_nm, 
            spot_radius= spot_radius)

    # detect dense regions
    regions_to_decompose, spots_out_regions, region_size = detection.get_dense_region(
        image= image, 
        spots= spots,
        voxel_size=voxel_size_nm,
        spot_radius= spot_radius,
        beta= beta)

    # precompute gaussian function values
    max_grid = max(200, region_size + 1)
    precomputed_gaussian = detection.precompute_erf(
        ndim= ndim,
        voxel_size=voxel_size_nm,
        sigma=(sigma_z, sigma_yx, sigma_yx),
        max_grid=max_grid)

    # simulate gaussian mixtures
    try :
        spots_in_regions, _ = detection.simulate_gaussian_mixture(
            image= image,
            candidate_regions=regions_to_decompose,
            voxel_size= voxel_size_nm,
            sigma=(sigma_z, sigma_yx, sigma_yx),
            amplitude=amplitude,
            background=background,
            precomputed_gaussian=precomputed_gaussian)
    
    except ValueError as error :
        if "need at least one array to concatenate" in str(error) :
            raise NoSpotError("No dense regions have been found for deconvolution.")
        else : raise error
    except Exception  as error :
        raise error

    spots_postdecomp = np.concatenate((spots_out_regions, spots_in_regions[:, :3]), axis=0)

    return spots_postdecomp



def remove_mean_gaussian_background(image, sigma = 5, voxel_size = (1,1,1)):
    """Removes background of the image using gaussian filter. If the input image is 3D the mean gaussian background is computed from all slices and then substract to the stack.
    
    Parameters
    ----------
        image : np.ndarray
    
        sigma : scalar
            Defines the gaussian filter intensity.
    Returns
    -------
        image_no_background : np.ndarray
            Image with background substracted.
    """

    if not len(voxel_size) == image.ndim : raise ValueError("Inconsistency between voxel_size length and image dimension.")
    anisotropy_coef = compute_anisotropy_coef(voxel_size=voxel_size)
    corrected_sigma = [sigma / anisotropy_i for anisotropy_i in anisotropy_coef]

    image_no_background = stack.remove_background_gaussian(image, corrected_sigma)

    return image_no_background


def remove_artifact(deconvoluted_spots, artifact_radius, voxel_size , spot_density = 2) :

    """
    Artifact are detected as spherical clusters of radius 'artifact_size' and with an average density 'spot_density' of spot within the cluster.
    All spots within the artifact are then removed from deconvoluted_spos.
    
    Critical number of spot is computed as :
    >>> (total_pixel_approximation) * spot_density /100
    >>> with total_pixel_approximation = 4/3*pi*(artifact_radius_xy)²*artifact_radius_z ~rounded to unity

    Parameters
    ----------
        deconvoluted_spots : np.ndarray(z,y,x)
            A dense region decomposition is highly recommended prior to this function
        artifact_radius : int
            in nm
        voxel_size : tuple (z,y,x)
        spot_density : float 
            in range ]0,100]
    """
    
    if spot_density <= 0 or spot_density > 100 : raise ValueError("Spot density must be in range ]0,100]. Current value is {0}".format(spot_density))
    if len(deconvoluted_spots) == 0 : return deconvoluted_spots

    critical_spot_number = _compute_critical_spot_number(radius_nm= artifact_radius, voxel_size=voxel_size, density=spot_density)
    artifacts_df:pd.DataFrame = cluster_detection(deconvoluted_spots, voxel_size=voxel_size, radius= artifact_radius, nb_min_spots=critical_spot_number, keys_to_compute= ['clustered_spots_dataframe'])['clustered_spots_dataframe']
    drop_index = artifacts_df[~artifacts_df["cluster_id"].isna()].index
    artifacts_df = artifacts_df.drop(drop_index, axis= 0)

    clean_spots = get_centroids_list(artifacts_df)
    return np.array(clean_spots, dtype= int)

### Dataframes
def build_Spots_and_Cluster_df(detection_result : dict) :
    """
    Made to build Spots pandas dataframe from result of multi_threaded call to `multi_thread_full_detection`.

    Parameter
    ---------
        detection_result
    """

    SPOTS_COLUMNS = [
        'detection_id',
        'spots',
        'spots_post_decomp',
        'clustered_spots_dataframe',
        'clusters_dataframe',
        'clusters',
        'clustered_spots',
        'free_spots',   
    ]
    
    detection_res = pd.DataFrame(columns=SPOTS_COLUMNS, data=detection_result)
    detection_res = detection_res.set_index('detection_id', verify_integrity=True, drop=False)

    Spots = pd.DataFrame()
    Clusters = pd.DataFrame()

    for detection_id in detection_res.index :
        spots: pd.DataFrame = detection_res.at[detection_id, 'clustered_spots_dataframe']
        clusters: pd.DataFrame = detection_res.at[detection_id, 'clusters_dataframe']

        if isinstance(spots, (int,float)) : 
            spots = pd.DataFrame(columns=['id','cluster_id','population', 'detection_id']) 
        if isinstance(clusters, (int,float)) : 
            clusters = pd.DataFrame(columns=['id','population', 'detection_id']) 
        
        #names
        spots = spots.rename(columns={'id' : 'spot_id'})
        clusters = clusters.rename(columns={'id' : 'cluster_id'})
        
        #is_clustered
        spots['cluster_id'] = spots['cluster_id'].replace(-1,np.nan)
        clusters['cluster_id'] = clusters['cluster_id'].replace(-1,np.nan)
        spots['population'] = spots['cluster_id'].isna()
        spots['population'] = spots['population'].replace({True : 'free', False : 'clustered'})

        #Detection_id
        spots['detection_id'] = detection_id
        clusters['detection_id'] = detection_id

        if not spots.empty :
            Spots = pd.concat([
                Spots,
                spots,
            ],axis=0)

        if not clusters.empty :
            Clusters = pd.concat([
                Clusters,
                clusters,
            ],axis=0)

    Spots = Spots.reset_index(drop=True)
    Clusters = Clusters.reset_index(drop=True)

    return Spots, Clusters

def _compute_clustered_spots_dataframe(clustered_spots) :
    if len(clustered_spots) == 0 : return pd.DataFrame(columns= ["id", "cluster_id", "z", "y", "x"])
    z, y ,x, cluster_index = list(zip(*clustered_spots))
    ids = np.arange(len(clustered_spots))

    df = pd.DataFrame({
        "id" : ids
        ,"cluster_id" : cluster_index
        ,"z" : z
        ,"y" : y
        ,"x" : x
    })

    null_idx = df[df['cluster_id'] == -1].index
    df.loc[null_idx, 'cluster_id'] = np.nan

    return df

def _compute_cluster_dataframe(clusters) :
    if len(clusters) == 0 : return pd.DataFrame(columns= ["id", "z", "y", "x", "spot_number"])
    z, y, x, spots_number, cluster_index = list(zip(*clusters))

    df = pd.DataFrame({
        "id" : cluster_index
        ,"z" : z
        ,"y" : y
        ,"x" : x
        ,"spot_number" : spots_number
    })

    return df


### CUSTOM ERROR CLASS
def detectiontimeout_handler(signum, frame):
    raise DetectionTimeOutError('Acquisition processing timeout.')


class DetectionError(Exception) :
    """Exception class raised during spot detection."""
    pass

class DetectionTimeOutError(DetectionError):
    """'DetectionError' Subclass raised when detection takes too much time."""
    pass

class NoSpotError(DetectionError):
    """'DetectionError' Subclass raised when no spots or too few spots have been detected."""
    pass

class TooManySpotsError(DetectionError):
    """'DetectionError' Subclass raised when too much spots have been detected."""
    pass
