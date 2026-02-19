from .types import parameters_dict 

################

# 0. Global

################
#fullpath to main folder given by experimentalist
RUN_PATH = "/media/SSD_floricslimani/Fish_seq/Davide/FakeRun/"
VOXEL_SIZE = (200,97,97) # size of a pixel in nanometer (z,y,x)
WASHOUT_KEY_WORD = 'Washout' #key for washout in gene map excel
HAS_BEAD_CHANNEL = False

################

# 1. Input

################
# MAP_FILENAME = "HeLa-POLR2_Run9.xlsx"  #filename of required map file giving cycles names
MAP_FILENAME = "cycle_file.xlsx"  #filename of required map file giving cycles names
# MAP_FILENAME = "HeLa-Puro-POLR2_Run10.xlsx"  #filename of required map file giving cycles names
CYCLE_KEY = "Cycle n."
GENES_NAMES_KEY = ["Gene0", ] # Ordered by channel


cycle_regex = r"img(\d+)_000_000000_0000000000.ome.tif" #regex to catch cycle number from tif filename.
FOLDER_KEYS = { # folder names where nucleus and fish images can be found (required keys : 'nucleus', 'fish')
    'nucleus_folder' : "DAPI_Z-stacks",
    'fish_folder' : "FISH_Z-stacks",
}



################

# 2. Segmentation

################

MODEL_DICT = {#cellpose model names, required keys : 'nucleus' and 'cytoplasm'
    'nucleus_model' : 'nuclei',
    'cytoplasm_model' : 'cyto3'
}

OBJECT_SIZE_DICT = {#size in px given to cellpose for prediction
    'nucleus_size' : 140,
    'cytoplasm_size' : 200
}

PLOT_VISUALS = True

################

# 3. Detection

##############
detection_MAX_WORKERS = 4 #Number of threads that can work simultaneously for detection; intensive process so recommended to not excess CPU cores number (currently : 4)

SPOT_SIZE = (300,140,140) #expected size of single molecules in nanometers (not less than voxel size)

#Big fish parameter for spot deconvolution;
ALPHA = 0.5  # aims at building reference spot as the 'alpha-percentile' of spot distribution. ex alpha= 0.5 --> reference spot = distribution median
BETA = 1 # aims at choosing regions to perform deconvolution on. Algorithm looks for regions with bright pixel connected. beta defines what bright pixel are with test : pixel_intensity >= MEDIAN_spot_intensity * beta (independantly of alpha)
GAMMA = 2 # size of kernel for gaussian blur performed before deconvolution.

CLUSTER_SIZE = 400 #size of cluster in nanometer (radius for DBScan algorithm)
MIN_SPOT_PER_CLUSTER = 5 #

ARTIFACT_RADIUS = 1400 # Radius of spots artifact to remove in nanometers.
DETECTION_SLICE_TO_REMOVE = [5,None] # number of slice you want to remove bottom/top If you don't want to remove use None. ie : 1,None removes one slice from the bottom and none from the top


################

# 4. Drift

##############
SAVE_PATH = RUN_PATH + '/visuals/'
BEAD_SIZE = (200, 200, 200) #size of fluorescent beads in nanometers used for fov aligment
DRIFT_SLICE_TO_REMOVE = [5,5] # Number of slice to remove to avoid detecting noise
DO_HIGHPASS_FILTER = False






################

# 4. Quantification

##############
COLOC_DISTANCE = 200 #distance to consider for colocalization events in nanometers
quantif_MAX_WORKERS = 10 #Number of threads to use while computing cells features (small_process, currently good performance with 10)
COLOC_POPULATION = ('all', 'free', 'clustered') # population to consider when computing colocalisation


#### Integrity check ####

if 0 in DETECTION_SLICE_TO_REMOVE : raise ValueError("Error with 'DETECTION_SLICE_TO_REMOVE' parameter. 0 value is not allowed, if you don't want to remove slice use None instead.")
if 0 in DRIFT_SLICE_TO_REMOVE : raise ValueError("Error with 'DRIFT_SLICE_TO_REMOVE' parameter. 0 value is not allowed, if you don't want to remove slice use None instead.")


def get_default_settings() :
    default_settings = {
        "RUN_PATH" : RUN_PATH,
        "VOXEL_SIZE" : VOXEL_SIZE,
        "WASHOUT_KEY_WORD" : WASHOUT_KEY_WORD,
        "HAS_BEAD_CHANNEL" : HAS_BEAD_CHANNEL,
        "MAP_FILENAME" : MAP_FILENAME,
        "CYCLE_KEY" : CYCLE_KEY,
        "GENES_NAMES_KEY" : GENES_NAMES_KEY,
        "cycle_regex" : cycle_regex, 
        "FOLDER_KEYS" : FOLDER_KEYS,
        "MODEL_DICT" : MODEL_DICT,
        "OBJECT_SIZE_DICT" : OBJECT_SIZE_DICT,
        "PLOT_VISUALS" : PLOT_VISUALS, 
        "detection_MAX_WORKERS" : detection_MAX_WORKERS, 
        "SPOT_SIZE" : SPOT_SIZE,
        "ALPHA" : ALPHA, 
        "BETA" : BETA, 
        "GAMMA" : GAMMA, 
        "CLUSTER_SIZE" : CLUSTER_SIZE, 
        "MIN_SPOT_PER_CLUSTER" : MIN_SPOT_PER_CLUSTER, 
        "ARTIFACT_RADIUS" : ARTIFACT_RADIUS, 
        "DETECTION_SLICE_TO_REMOVE" : DETECTION_SLICE_TO_REMOVE,
        "SAVE_PATH" : SAVE_PATH, 
        "BEAD_SIZE" : BEAD_SIZE,
        "DRIFT_SLICE_TO_REMOVE" : DRIFT_SLICE_TO_REMOVE,
        "DO_HIGHPASS_FILTER" : DO_HIGHPASS_FILTER, 
        "COLOC_DISTANCE" : COLOC_DISTANCE, 
        "quantif_MAX_WORKERS" : quantif_MAX_WORKERS, 
    }

    default_settings = parameters_dict(**default_settings)
    return default_settings