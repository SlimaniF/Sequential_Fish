# Keys of DataFrames used during pipeline

## Acquisition

        *acquisition_id	
        *location	
        *cycle	
        *full_path	
        *dapi_full_path	
        Gene1 (Ch1 - Cy5)	
        Gene2 (Ch2 - Cy3)	
        Barcode 1	
        Barcode 2
        threshold_i # for i in colors

## Detection
        *'detection_id
        *'acquisition_id' : sub_data['acquisition_id']
        'visual_name' : [visual_path] * image_number
        'filename' : list(sub_data['full_path'])
        'voxel_size' : [tuple(VOXEL_SIZE)] * image_number
        'spot_size' : [tuple(SPOT_SIZE)] * image_number
        'alpha' : [ALPHA] * image_number
        'beta' : [BETA] * image_number
        'gamma' : [GAMMA] * image_number
        'artifact_radius' : [ARTIFACT_RADIUS] * image_number
        'cluster_size' : [CLUSTER_SIZE] * image_number
        'min_spot_per_cluster' : [MIN_SPOT_PER_CLUSTER] * image_number
        'color_id'
        'visual_name'
        'threshold
        'image_path'
        'image_key'
        *location -> quantification.py

## Spots
        "spot_id" : ids
        ,"detection_id"
        ,"cluster_id" : cluster_index
        ,"z" : z
        ,"y" : y
        ,"x" : x
        ,"population"
        ,"intensity"
        ,*quantification.py* -> "drifted_z"
        ,*quantification.py* -> "drifted_y"
        ,*quantification.py* -> "drifted_x"

## Clusters
        "cluster_id" : cluster_index
        ,"detection_id"
        ,"z" : z
        ,"y" : y
        ,"x" : x
        ,"spot_number" : spots_number

## Drift
        'acquisition_id',
        'drift_type'
        'drift_z',
        'drift_y',
        'drift_x',
        'ref_bead_threshold',
        'drift_bead_threshold',
        'ref_bead_number',
        'drift_bead_number',
        'found_symetric_drift',
        'voxel_size',
        'bead_size'
        *location -> quantification.py

## Cell
        BigFish features
        acquisition_id
        detection_id
        cluster_number
        rna_number
        cell_coordinates
        label
        bbox
        nucleus_area_px
        nucleus_area_nm
        ["nucleus_mip_mean_signal","nucleus_mip_max_signal","nucleus_mip_min_signal","nucleus_mip_median_signal",
                           "nucleus_mean_mean_signal","nucleus_mean_max_signal","nucleus_mean_min_signal","nucleus_mean_median_signal"]



## Colocalisation
        location
        label
        spot1_number
        spot2_number
        object1 #renamed to detection_id1
        object2 #renamed to detection_id2
        spot1_total_number
        spot2_total_number
        population1
        population2
        count
        fraction
        sub_fraction
        cell_id

