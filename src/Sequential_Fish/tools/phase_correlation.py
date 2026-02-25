import numpy as np

import smfishtools.preprocessing 
from scipy.fft import fftn, fftshift, ifftn, ifftshift
from skimage.registration import phase_cross_correlation

def fft_phase_correlation_drift(
        reference_image : np.ndarray,
        drifted_image : np.ndarray,
        bead_size : tuple = None,
        voxel_size : tuple = None,
        upsample_factor = 1,
        disambiguate = False,
        highpass_filter = False,
) :
    """
    Wrapper around skimage.registration.phase_cross_correlation.
    Before runing phase cross correlation a high pass filter is used to align images with information from object of size `bead_size`.
    Uses correlation of phase in Fourier space to align images computation efficient.
    Outputs dict compatible to transform to pandas DataFrame.

    Parameters
    ----------
        * reference_image : image used as reference to align to.

        * drifted_image : image to find drift correction for.
        
        * upsample_factor : Upsampling factor. Images will be registered to within 1 / upsample_factor of a pixel. For example upsample_factor == 20 means the images will be registered within 1/20th of a pixel. Default is 1 (no upsampling).
        
        * disambiguate : The shift returned by this function is only accurate modulo the image shape, due to the periodic nature of the Fourier transform. If this parameter is set to True, the real space cross-correlation is computed for each possible shift, and the shift with the highest cross-correlation within the overlapping area is returned.

        * highpass_filter : Perform highpass filter in fourier space to keep beads information, requires bead size and voxel size.
        
    Returns
    -------
        Results : dict
            keys :
            * 'drift_z' : int
            * 'drift_y' : int
            * 'drift_x' : int
            * 'error' : float; Translation invariant normalized RMS error between reference_image and drifted_image.
            * 'phasediff' : float; Global phase difference between the two images.

    """
    
    shape = reference_image.shape

    if not bead_size is None :
        bead_spatial_frequencies = np.array([1/object_s for object_s in bead_size])

    #Creating bead shaped frequency filter
    if highpass_filter and not bead_size is None and not voxel_size is None: 

        space = 'fourier'

        min_frequency = min(bead_spatial_frequencies)
        bead_filter = Fourier_frequency_filter(
            shape=shape,
            voxel_size=voxel_size,
            min_frequency=min_frequency,
            max_frequency=None
            )

        #Computing fft and centering frequencies at the image center
        reference_image_fft = fftn(reference_image)
        reference_image_fft = fftshift(reference_image_fft)
        drifted_image_fft = fftn(drifted_image)
        drifted_image_fft = fftshift(drifted_image_fft)

        #Applying filter
        reference_image_fft *= bead_filter
        drifted_image_fft *= bead_filter

    else :
        reference_image_fft = reference_image
        drifted_image_fft = drifted_image
        space = 'real'

    shift, error, phasediff = phase_cross_correlation(
        reference_image=reference_image_fft,
        moving_image=drifted_image_fft,
        space= space,
        upsample_factor=upsample_factor,
    )
    if len(shift) == 3 :
        drift_z, drift_y, drift_x = shift
    elif len(shift) == 2 :
        drift_y, drift_x = shift
        drift_z = np.nan

    else :
        raise ValueError("Incorrect number of dimensions in shift : {0}. Should be 2 or 3".format(len(shift)))

    results = {
        "drift_z" : [round(drift_z)],
        "drift_y" : [round(drift_y)],
        "drift_x" : [round(drift_x)],
        "error" : [error],
        "phasediff" : [phasediff],
    }

    return results

def Fourier_frequency_filter(shape, voxel_size, min_frequency=None, max_frequency=None) :
    frequency_map = _frequency_map(shape, voxel_size)

    if type(min_frequency) == type(None) : min_frequency = 0
    if type(max_frequency) == type(max_frequency) : max_frequency = np.max(frequency_map)

    frequency_filter = (frequency_map >= min_frequency) & (frequency_map <= max_frequency)
    return frequency_filter

def _frequency_map(shape, voxel_size) :
    """
    Euclidian distance frequency map
    """
    voxel_frequency = tuple([1/(voxel_len*axis_len/2)] for axis_len, voxel_len in zip(shape,voxel_size))
    center_coordinates = (np.array(shape) -1) / 2 # float center 
    frequency_map = np.indices(shape, dtype = float)
    for axis_index, axis_frq in enumerate(voxel_frequency) :
        frequency_map[axis_index] = -frequency_map[axis_index] + center_coordinates[axis_index]
        frequency_map[axis_index] *= axis_frq
    frequency_map = np.square(frequency_map)
    frequency_map = np.sum(frequency_map, axis=0)
    frequency_map = np.sqrt(frequency_map)

    return frequency_map