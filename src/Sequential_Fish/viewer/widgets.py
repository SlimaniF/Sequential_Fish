"""
Submodule containing custom class for napari widgets
"""
import os
from typing import cast

import numpy as np
import pandas as pd
import napari

from napari.types import LayerDataTuple
from magicgui import magicgui

from ..tools.utils import open_all_locations_one_cycle, safe_merge_no_duplicates
from .utils import open_segmentation
from ..customtypes import NapariWidget
from ..customtypes import table_dict_type


from smfishtools.preprocessing.alignement import shift_array
from ._density import multichannel_clustering, spot_count_map



##  Load tab
_LOAD_WIDGETS : 'list[NapariWidget]' = []
def register_load_widget(cls) :
    _LOAD_WIDGETS.append(cls)
    return cls

def initiate_load_widgets(
        table_dict : table_dict_type, 
        voxel_size : tuple,
        color_table : dict,
        run_path : str
) -> 'list[NapariWidget]':
    widget_list = []
    for cls in _LOAD_WIDGETS :
        instance = cls(table_dict= table_dict, voxel_size=voxel_size, color_table=color_table, run_path=run_path)
        if hasattr(instance,"enabled") :
            if instance.enabled :
                widget_list.extend(instance.get_widgets())
        else :
            widget_list.extend(instance.get_widgets())

    return widget_list

##  Location tab
_LOCATION_WIDGETS = []
def register_location_widget(cls) :
    _LOCATION_WIDGETS.append(cls)
    return cls

def initiate_location_widgets(
        *,
        tables_dict,
        Viewer,
        linked_widgets,
) :
    widget_list = []
    for cls in _LOCATION_WIDGETS :
        widget_list.extend(cls(table_dict=tables_dict, Viewer=Viewer, linked_widgets=linked_widgets).get_widgets())
    return widget_list


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

#######
# INDIVIDUAL WIDGETS
#######

#Load data widgets
@register_load_widget
class SpotsLoader(NapariWidget) :
    """
    Allow user to load detected spots as layer points
    """
    def __init__(    
        self, 
        table_dict : table_dict_type,
        voxel_size :tuple, 
        color_table,
        **_
    ) :
        self.enabled=True
        self.Spots = table_dict.get('Spots')
        self.Detection = table_dict.get('Detection')

        if self.Spots is None or self.Detection is None :
            print("Disabling spots loader (Detection results not found)")
            self.enabled=False

            return None

        self.Acquisition = table_dict['Acquisition'].loc[:,['acquisition_id','location','cycle']]
        self.Detection = safe_merge_no_duplicates(
            self.Detection,
            self.Acquisition,
            on= "acquisition_id",
            keys= ["location", "cycle"]
        )
        self.Gene_map = table_dict['Gene_map'].loc[:,['cycle','color_id','target']]
        
        self.update(list(self.Acquisition['location'].unique()))

        self.voxel_size = voxel_size
        self.color_table = color_table
        super().__init__()

    def update(self, locations) :


        data = safe_merge_no_duplicates(
            self.Spots,
            self.Detection,
            on= 'detection_id',
            keys= ["location", "cycle", "color_id","acquisition_id"]
        )

        data = pd.merge(
            data,
            self.Acquisition.loc[self.Acquisition['location'].isin(locations), ['acquisition_id']],
            on= 'acquisition_id',
            validate='m:1',
        )

        data = pd.merge(
            data,
            self.Gene_map.loc[:,['cycle','color_id','target']],
            on= ['cycle','color_id'],
            how='left'
        )

        assert not any(data['target'].isna()), "Missing values for `target` in Spots. Merge is incomplete."
        self.data = data
        self.populations = ['all'] + list(data['population'].unique()) 
        self.target = list(data['target'].unique())

    def _create_widget(self) :
        @magicgui(
            target={"choices":self.target},
            population={"choices" : self.populations},
            drift_correction={
                    "widget_type" : "CheckBox",
                    "text" : "drift correction",
                    "value" : True,
                    },
            call_button= 'Load spots',
            auto_call=False
        )
        def load(target, population, drift_correction) -> LayerDataTuple :
            
            if drift_correction : 
                name = "{1}_{0}_spots_corrected".format(target, population)
                symbol = 'disc'
                z_indexer = 'z'
                y_indexer = 'y'
                x_indexer = 'x'
            else :
                name = "{1}_{0}_spots_drifted".format(target, population)
                symbol = 'x'
                z_indexer = 'drifted_z'
                y_indexer = 'drifted_y'
                x_indexer = 'drifted_x'

            if population == 'all' :
                sub_Detec = self.data.loc[self.data['target'] == target]
            else :
                sub_Detec = self.data.loc[(self.data['target'] == target) & (self.data['population'] == population)]

            #Fetch color
            color = self.color_table[self.color_table['target'] == target]['color']
            assert len(color) == 1, "Gene_map has non unique targets."
            color = color.iat[0]

            #Fetch spots
            spots_array = np.empty(shape=(0,4),dtype=int)
            for location_index, location  in enumerate(sub_Detec['location'].unique()) :

                spot_data = sub_Detec.loc[sub_Detec['location'] == location]
                C = [location_index] * len(spot_data)
                Z = spot_data[z_indexer]
                Y = spot_data[y_indexer]
                X = spot_data[x_indexer]

                spots = np.array(
                    list(zip(C,Z,Y,X)),
                    dtype=int,
                )

                spots_array = np.concatenate([spots_array, spots])

            layerdata = cast(LayerDataTuple, (spots_array, 
                         {
                             "scale" : self.voxel_size,
                             "size" : 10, 
                             "name" : name, 
                             'ndim' : 4, 
                             'face_color' : '#0000' ,
                             'border_color' : color, 
                             'blending' : 'additive',
                             'symbol' : symbol
                             },
                        'Points'))
            return layerdata
        return load

@register_load_widget    
class ClustersLoader(NapariWidget) :
    """
    Allow user to load detected cluster as Points layer.
    """

    def __init__(
            self, 
            table_dict : table_dict_type,
            voxel_size :tuple, 
            color_table,
            **_
            ):
        
        self.enabled=True
        self.Spots = table_dict.get('Spots')
        self.Detection = table_dict.get('Detection')
        self.Clusters = table_dict.get('Clusters')

        if self.Spots is None or self.Detection is None or self.Clusters is None :
            print("Disabling cluster loader : (Detection results not found)")
            self.enabled=False

            return None

        self.Detection = table_dict['Detection']
        self.Acquisition = table_dict['Acquisition'].loc[:,['acquisition_id','cycle','location']]
        self.Gene_map = table_dict['Gene_map'].loc[:,['cycle','color_id','target']]
        self.Detection = safe_merge_no_duplicates(
            self.Detection,
            self.Acquisition,
            on= "acquisition_id",
            keys= ["location", "cycle"]
        )
        self.update(list(self.Acquisition['location'].unique()))

        self.voxel_size = voxel_size
        self.color_table = color_table
        super().__init__()

    def update(self, locations) :

        data = safe_merge_no_duplicates(
            self.Clusters,
            self.Detection,
            on= 'detection_id',
            keys= ["location", "cycle", "color_id","acquisition_id"]
        )

        data = pd.merge(
            data,
            self.Acquisition.loc[self.Acquisition['location'].isin(locations),['acquisition_id']],
            on= 'acquisition_id',
            validate='m:1',
        )

        data = pd.merge(
            data,
            self.Gene_map.loc[:,['cycle','color_id','target']],
            on= ['cycle','color_id'],
            how='left'
        )

        assert not any(data['target'].isna()), "Missing values for `target` in Spots. Merge is incomplete."

        self.data = data
        self.target = list(self.data['target'].unique())


    def _create_widget(self) :
        @magicgui(
            target={"choices":self.target},
            drift_correction={
                    "widget_type" : "CheckBox",
                    "text" : "drift correction",
                    "value" : True,
                    },
            call_button= 'Load clusters',
            auto_call=False
        )
        def load(target, drift_correction) -> LayerDataTuple :
            
            if drift_correction : 
                name = "{0}_clusters_corrected".format(target)
                symbol = "diamond"
                z_indexer = 'z'
                y_indexer = 'y'
                x_indexer = 'x'
            else :
                name = "{0}_clusters_drifted".format(target)
                symbol = "clobber"
                z_indexer = 'drifted_z'
                y_indexer = 'drifted_y'
                x_indexer = 'drifted_x'

            sub_data = self.data.loc[self.data['target'] == target]
            
            #Fetch color
            color = self.color_table[self.color_table['target'] == target]['color']
            assert len(color) == 1, "Gene_map has non unique targets."
            color = color.iat[0]


            #Fetch cluster centers
            spots_array = np.empty(shape=(0,4),dtype=int)
            for location_index, location  in enumerate(sub_data['location'].unique()) :
                spots_data = sub_data.loc[sub_data['location'] == location]
                C = [location_index] * len(spots_data)
                Z = spots_data[z_indexer]
                Y = spots_data[y_indexer]
                X = spots_data[x_indexer]

                spots = np.array(
                    list(zip(C,Z,Y,X)),
                    dtype=int,
                )

                spots_array = np.concatenate([spots_array, spots])
            layerdata = cast(
                LayerDataTuple,
                (spots_array, 
                         {"scale" : self.voxel_size,
                          "size" : 12, 
                          "name" : name, 
                          'ndim' : 4, 
                          'face_color' : color,
                          'symbol' : symbol, 
                          'blending' : 'additive'}
                          , 'Points')
            )
            return layerdata
        return load

@register_load_widget
class SignalLoader(NapariWidget) :
    def __init__(
            self, 
            table_dict : table_dict_type,
            voxel_size :tuple, 
            color_table : pd.DataFrame,
            **_
            ):
        
        #Table
        self.Gene_map = table_dict['Gene_map']

        Drift = table_dict['Drift']
        self.Drift = Drift.loc[:,['acquisition_id', 'drift_z', 'drift_y', 'drift_x']]
        
        self.Acquisition = table_dict['Acquisition'].set_index(["location","cycle"], verify_integrity=True)
        
        self.update(list(self.Acquisition.index.get_level_values(0).unique()))

        self.color_table = color_table.set_index("target",verify_integrity=True)
        self.voxel_size = voxel_size
        self.has_beads = not cast(bool,self.Acquisition['bead_channel'].isna().all())
        super().__init__()

    def update(self, locations) :

        self.data = pd.merge(
            self.Acquisition[self.Acquisition.index.get_level_values(0).isin(locations)].reset_index(drop=False),
            self.Drift,
            on='acquisition_id'
        )
        self.data = self.data.set_index(['location',"cycle"]).sort_index()
        self.target = list(self.Gene_map['target'].unique())

    def _create_widget(self) :
        @magicgui(
                target = {'choices' : self.target},
                drift_correction={
                    "widget_type" : "CheckBox",
                    "text" : "drift correction",
                    "value" : True,
                    },
                call_button="Load signal",
                signal_type = {
                    "widget_type" : "RadioButtons",
                    "choices" : ['fish','dapi','beads'] if self.has_beads else ['fish', 'dapi'],
                    "value" : "fish"
                },
                auto_call=False
                )
        def load(target, drift_correction, signal_type) -> LayerDataTuple:
            data = self.Gene_map.loc[self.Gene_map['target'] == target].iloc[0]
            color = self.color_table.at[target, "colormaps"]
            cycle, color_id = data['cycle'], int(data['color_id'])

            if signal_type == "fish" :
                channel_index = color_id
                name = f"{target}_fish"
            elif signal_type == "dapi" :
                channel_index = self.Acquisition['dapi_channel'].iat[0]
                name = f"dapi_signal_cycle{cycle}"
            elif signal_type == "beads" :
                channel_index = self.Acquisition['bead_channel'].iat[0]
                name = f"beads_signal_cycle{cycle}"
            else :
                raise NotImplementedError("Unimplemented choice")
            

            if drift_correction :
                name += "_corrected"
            else :
                name += "_signal_drifted"


            array = open_all_locations_one_cycle(
                self.Acquisition.reset_index(drop=False),
                cycle=cycle,
            )
            array = array[..., channel_index]
                
            if drift_correction :
                location_list = self.Acquisition.index.get_level_values(0).unique().to_list()
                location_list.sort()
                assert len(location_list) == len(array)

                location_index=0
                for location, stack in zip(location_list, array) :
                    drift = self.data.loc[(location,cycle),["drift_z","drift_y","drift_x"]].astype(int)
                    array[location_index] = shift_array(stack, *drift)

            layerdata = cast(LayerDataTuple, (
                array,
                {
                    "scale" : self.voxel_size, 
                    "name" : name, 
                    'blending' : 'additive', 
                    'colormap' : color if signal_type != "dapi" else "blue"},
                'Image'
            )
            )

            return layerdata

        return load

@register_load_widget    
class SegmentationLoader(NapariWidget) :
    
    def __init__(
            self,
            run_path : str,
            voxel_size :tuple,
            table_dict : table_dict_type, 
            segmentation_folder_name:str = "/segmentation/",
            **_
            ):
        
        self.enabled=True


        self.Acquisition = table_dict['Acquisition']

        self.update(list(self.Acquisition['location'].unique()))

        self.segmentation_fullpath = run_path + segmentation_folder_name

        if not os.path.isdir(self.segmentation_fullpath) or len(os.listdir(self.segmentation_fullpath)) == 0 :
            print("Disabling segmentation masks (Segmentation results not found)")
            self.enabled=False

        self.voxel_size = voxel_size
        super().__init__()

    def update(self,locations) :
        self.data = self.Acquisition.loc[self.Acquisition['location'].isin(locations)]

    def _create_widget(self) :

        @magicgui(
                call_button= "Load segmentation",
                target={
                    "widget_type" : "RadioButtons",
                    "choices" : ["nucleus","cytoplasm"],
                    "orientation" : "horizontal",
                    "value" : "nucleus",
                    "label" : " ",
                },
                auto_call=False
        )
        def load_segmentation(target) -> LayerDataTuple:

            shape = np.array(list(self.Acquisition['fish_shape']),dtype=int)
            shape = np.max(shape, axis=0)
            z_size = shape[0]
            name = "{0}_mask".format(target)
            locations = list(self.data.sort_values('location')['location'].unique())
            masks = open_segmentation(
                self.segmentation_fullpath, 
                locations, 
                target=target, 
                z_repeat= z_size
                ) #masks list sorted on Acquisition['location']

            layerdata = LayerDataTuple((
                masks,
                {"scale" : self.voxel_size, "name" : name, "blending" : "additive"},
                'Labels'
            ))


            return layerdata
        return load_segmentation

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

#Location widget
@register_location_widget
class LocationSelector(NapariWidget) :
    def __init__(
            self, 
            table_dict: table_dict_type, 
            Viewer : napari.Viewer, 
            linked_widgets : list,
            **_
            ):
        self.Full_Acquisiton = table_dict['Acquisition'].copy()
        self.location_choices = list(self.Full_Acquisiton['location'].unique())
        self.selection = self.location_choices.copy()
        self.Viewer = Viewer
        self.linked_widgets = linked_widgets
        super().__init__()
        
    def update_location(self) :
        for layer in self.Viewer.layers.copy() :
            self.Viewer.layers.remove(layer)
        self.Viewer.reset_view()

        for widget in self.linked_widgets : 
            widget.update(self.selection)
            widget.widget.update()

    def _create_widget(self) :
        @magicgui(
            selected_location={
                "widget_type" : "Select",
                "choices" : self.location_choices,
                "value" : self.location_choices,
                "label" : ' ',
            },
            call_button= "select locations"
        )
        def select_location(selected_location) :
            print("Selected locations : ", selected_location)
            self.selection = selected_location
            self.update_location()
        
        
        return select_location
