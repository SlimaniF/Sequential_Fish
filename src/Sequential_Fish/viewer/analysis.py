"""
Widgets for analyzing data in napari viewer
"""
from typing import cast

import pandas as pd

from napari.types import LayerDataTuple
from magicgui import magicgui

from ..customtypes import NapariWidget
from ..customtypes import table_dict_type


from ._density import multichannel_clustering, spot_count_map

##  Analysis tab
_ANALYSIS_WIDGETS = []
def register_analysis_widget(cls) :
    _ANALYSIS_WIDGETS.append(cls)
    return cls

def initiate_analysis_widgets(
        voxel_size : tuple,
        table_dict : table_dict_type,
        color_table : dict,

) :
    widget_list = []
    for cls in _ANALYSIS_WIDGETS :
        instance = cls(voxel_size=voxel_size, table_dict=table_dict, color_table=color_table)
        if hasattr(instance,"enabled") :
            if instance.enabled :
                widget_list.extend(instance.get_widgets())
        else :
            widget_list.extend(instance.get_widgets())
    return widget_list


## Analysis widgets
@register_analysis_widget
class MultichannelCluster(NapariWidget) :
    def __init__(
            self, 
            table_dict, 
            voxel_size,
            **_
            ):

        self.enabled=True
        self.Spots = table_dict.get('Spots')
        self.Detection = table_dict.get('Detection')

        if self.Spots is None or self.Detection is None :
            self.enabled=False
            print("Disabling MultichannelCluster(Detection results not found)")

            return None

        self.ref_Acquisition = table_dict['Acquisition']
        self.Detection = table_dict['Detection']
        self.Spots = table_dict['Spots']
        self.Gene_map = table_dict['Gene_map']
        self.voxel_size = voxel_size
        self.update(list(table_dict['Acquisition']['location'].unique()))
        super().__init__()

    def update(self, locations) :
        self.Acquisition = self.ref_Acquisition.loc[self.ref_Acquisition['location'].isin(locations)]

    def _create_widget(self) :
        @magicgui(
                cluster_radius = {
                    "widget_type" : "SpinBox",
                    "value" : max(self.voxel_size),
                    "min" : 0,
                    "max" : 100 * max(self.voxel_size),
                    "label" : "cluster radius (nm) :",
                },
                min_spot_number = {
                    "widget_type" : "SpinBox",
                    "min" : 0,
                    "max" : 100,
                    "value" : 4,
                    "label" : "min spots number :",
                },
                min_channel_number = {
                    "widget_type" : "SpinBox",
                    "min" : 1,
                    "max" : len(self.Gene_map['target'].unique()),
                    "value" : 1,
                    "label" : "min rna number :",
                },
                call_button= "multichannel DBSCAN"
        )
        def multichannel_DBSCAN(cluster_radius, min_spot_number, min_channel_number) -> LayerDataTuple :
            multichannel_clusters = multichannel_clustering(
                self.Acquisition,
                self.Detection,
                self.Spots,
                self.Gene_map,
                voxel_size= self.voxel_size,
                cluster_radius=cluster_radius,
                nb_min_spots= min_spot_number,
                no_filtering=False,
            )

            multichannel_clusters: pd.DataFrame = multichannel_clusters.loc[multichannel_clusters['unique_target_number'] >= min_channel_number]

            c = []
            for location_index, location in enumerate(multichannel_clusters.index.get_level_values(0).unique()) :
                c.extend([location_index]* len(multichannel_clusters.index.get_level_values(0)[multichannel_clusters.index.get_level_values(0) == location]))

            z = list(multichannel_clusters['z'].astype(int))
            y = list(multichannel_clusters['y'].astype(int))
            x = list(multichannel_clusters['x'].astype(int))
            single_number = multichannel_clusters['single_molecule_number']
            target_names = multichannel_clusters['target_names']

            clusters = list(zip(c,z,y,x))
            clusters = pd.array(clusters, dtype=int)



            name = "multichannel_clusters_r{0}_n{1}".format(cluster_radius, min_spot_number)

            layer_data = cast(LayerDataTuple,
            (clusters, 
                         {"scale" : self.voxel_size, 
                          "name" : name, 
                          'ndim' : 4, 
                          'face_color' : 'white',
                          'symbol' : 'square',
                          'features' : {'single_number' : single_number, 'target_names' : target_names}, 
                          'size' : 0.1, 
                          'blending' : 'additive'}
                          , 'Points')
            )
            return layer_data
            
        
        return multichannel_DBSCAN

@register_analysis_widget
class SpotCountMapper(NapariWidget) :
    def __init__(
            self, 
            table_dict, 
            voxel_size,
            **_,
            ):
        
        self.enabled=True
        self.Spots = table_dict.get('Spots')
        self.Detection = table_dict.get('Detection')

        if self.Spots is None or self.Detection is None :
            self.enabled=False
            print("Disabling spots heatmap (Detection results not found)")

            return None

        self.ref_Acquisition = table_dict['Acquisition']
        self.Detection = table_dict['Detection']
        self.Spots = table_dict['Spots']
        self.Gene_map = table_dict['Gene_map']
        self.voxel_size = voxel_size
        self.update(list(table_dict['Acquisition']['location'].unique()))
        super().__init__()

    def update(self, locations) :
        self.Acquisition = self.ref_Acquisition.loc[self.ref_Acquisition['location'].isin(locations)]

    def _create_widget(self) :
        @magicgui(
                targets ={
                    "widget_type" : "Select",
                    "choices" : list(self.Gene_map['target'].sort_values().unique()),
                    "value" : list(self.Gene_map['target'].unique()),
                    "label" : " ",
                },
                call_button= "spot heatmap"
        )
        def generate_spot_count_map(targets) -> LayerDataTuple :
            Gene_map = self.Gene_map[self.Gene_map['target'].isin(targets)]
            spot_count_array = spot_count_map(
                Acquisition=self.Acquisition,
                Detection=self.Detection,
                Spots=self.Spots,
                Gene_map=Gene_map,
                no_filtering=False
            )

            layer_data = cast(LayerDataTuple,
            (
                spot_count_array,
                {"name" : "spot_count_map",
                 "scale" : self.voxel_size,
                 "blending" : "additive",
                 "colormap" : "inferno"
                 },
                 'Image'
            )
            )
            return layer_data
        return generate_spot_count_map