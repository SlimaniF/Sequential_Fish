frameon = True
FILTER_RNA = ['POLR2B_20']

RENAME_RULE = {
    "Washout_18_1" : "Washout_19_1",
    "POLR2A_16" : "POLR2A_end",
    "POLR2A_0" : "POLR2A",
}

# Distributions
distribution_measures = [
    'rna_number', 
    'cluster_number', 
    'proportion_rna_in_foci', 
    'nb_rna_in_nuc', 
    'index_mean_distance_nuc', 
    'index_mean_distance_cell'
    ]

# Density analysis
min_diversity = 3
min_spots_number = 3
cluster_radius = 400 #nm

#Co-localization analysis
coloc_distance = 400 #nm
coloc_significance=1e-4