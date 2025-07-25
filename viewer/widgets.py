"""
Submodule containing custom class for napari widgets
"""

import numpy as np
import pandas as pd
import napari
from typing import TypedDict
from tqdm import tqdm

from napari.types import LayerDataTuple
from magicgui import magicgui
from magicgui.widgets import Widget

from ..customtypes import NapariWidget
from ..tools.utils import open_image
from .utils import open_segmentation
from .utils import pad_to_shape
from .utils import reorder_image_stack
from .utils import reshape_stack
from .utils import correct_map


from pbwrap.preprocessing.alignement import shift_array
from ._density import multichannel_clustering, spot_count_map


class table_dict_type(TypedDict) :
    Acquisition : pd.DataFrame
    Detection : pd.DataFrame
    Spots : pd.DataFrame
    Clusters : pd.DataFrame
    Drift : pd.DataFrame
    Cell : pd.DataFrame
    Colocalisation : pd.DataFrame
    Gene_map : pd.DataFrame



##  Load tab
_LOAD_WIDGETS : 'list[NapariWidget]' = []
def register_load_widget(cls) :
    _LOAD_WIDGETS.append(cls)
    return cls

def initiate_load_widgets(
        table_dict : dict, 
        voxel_size : tuple,
        color_table : dict,
        run_path : str
) -> 'list[NapariWidget]':
    widget_list = []
    for cls in _LOAD_WIDGETS :
        widget_list.extend(cls(table_dict= table_dict, voxel_size=voxel_size, color_table=color_table, run_path=run_path).get_widgets())
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
        table_dict : dict,
        color_table : dict,

) :
    widget_list = []
    for cls in _ANALYSIS_WIDGETS :
        widget_list.extend(cls(voxel_size=voxel_size, table_dict=table_dict, color_table=color_table, ).get_widgets())
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
        **kwargs
    ) :

        self.Spots = table_dict['Spots']
        self.Detection = table_dict['Detection'].loc[:,['detection_id']]
        self.Acquisition = table_dict['Acquisition'].loc[:,['acquisition_id','location']]
        self.Gene_map = table_dict['Gene_map'].loc[:,['cycle','color_id','target']]
        
        self.update(list(self.Acquisition['location'].unique()))

        self.voxel_size = voxel_size
        self.color_table = color_table
        super().__init__()

    def update(self, locations) :


        data = pd.merge(
            self.Spots,
            self.Detection,
            on= 'detection_id',
            validate='m:1',
        )

        data = pd.merge(
            data,
            self.Acquisition[self.Acquisition['location'].isin(locations)],
            on= 'acquisition_id',
            validate='m:1',
        )

        data = pd.merge(
            data,
            self.Gene_map,
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

            layerdata = (spots_array, 
                         {
                             "scale" : self.voxel_size,
                             "size" : 0.1, 
                             "name" : name, 
                             'ndim' : 4, 
                             'face_color' : '#0000' ,
                             'edge_color' : color, 
                             'blending' : 'additive',
                             'symbol' : symbol
                             },
                        'Points')
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
            **kwargs
            ):
        
        self.Clusters = table_dict['Clusters']
        self.Detection = table_dict['Detection'].loc[:,['detection_id','color_id']]
        self.Acquisition = table_dict['Acquisition'].loc[:,['acquisition_id','cycle','location']]
        self.Gene_map = table_dict['Gene_map'].loc[:,['cycle','color_id','target']]

        self.update(list(self.Acquisition['location'].unique()))

        self.voxel_size = voxel_size
        self.color_table = color_table
        super().__init__()

    def update(self, locations) :

        data = pd.merge(
            self.Clusters,
            self.Detection,
            on= 'detection_id',
            validate='m:1',
        )

        data = pd.merge(
            data,
            self.Acquisition[self.Acquisition['location'].isin(locations)],
            on= 'acquisition_id',
            validate='m:1',
        )

        data = pd.merge(
            data,
            self.Gene_map,
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
            layerdata = (spots_array, 
                         {"scale" : self.voxel_size,
                          "size" : 0.2, 
                          "name" : name, 
                          'ndim' : 4, 
                          'face_color' : color,
                          'symbol' : symbol, 
                          'blending' : 'additive'}
                          , 'Points')
            return layerdata
        return load

@register_load_widget
class FishSignalLoader(NapariWidget) :
    def __init__(
            self, 
            table_dict : table_dict_type,
            voxel_size :tuple, 
            color_table : pd.DataFrame,
            **kwargs
            ):
        
        #Table
        self.Gene_map = table_dict['Gene_map']

        Drift = table_dict['Drift']
        self.Drift = Drift.loc[:,['acquisition_id', 'drift_z', 'drift_y', 'drift_x']]
        
        self.Acquisition = table_dict['Acquisition']
        
        self.update(list(self.Acquisition['location'].unique()))

        self.color_table = color_table
        self.voxel_size = voxel_size
        super().__init__()

    def update(self, locations) :
        self.data = pd.merge(
            self.Acquisition[self.Acquisition['location'].isin(locations)],
            self.Drift,
            on='acquisition_id'
        )
        self.data = self.data.sort_values(['location', 'full_path'])
        self.target = sorted(list(self.Gene_map['target'].unique()))

    def _create_widget(self) :
        @magicgui(
                target = {'choices' : self.target},
                drift_correction={
                    "widget_type" : "CheckBox",
                    "text" : "drift correction",
                    "value" : True,
                    },
                call_button="Load fish signal",
                auto_call=False
                )
        def load(target, drift_correction) -> LayerDataTuple:
            data = self.Gene_map.loc[self.Gene_map['target'] == target].iloc[0]
            color = self.color_table[self.color_table['target'] == target]['colormaps'].iat[0]
            cycle, color_id = data['cycle'], int(data['color_id'])
            
            if drift_correction :
                name = "{0}_fish_signal_corrected".format(target)
            else :
                name = "{0}_fish_signal_drifted".format(target)


            sub_Acqu = self.data.loc[self.data['cycle'] == cycle]

            if len(sub_Acqu) > 1 : 
                    shapes = np.array(list(sub_Acqu['fish_shape']))
                    max_shape = np.max(shapes, axis=0)
            else :
                max_shape = sub_Acqu['fish_shape'].iat[0]

            image_list = []
            for index in tqdm(sub_Acqu.index, desc="Opening fish signal ({0})".format(target)) :
                fullpath = sub_Acqu.at[index, "full_path"]
                shape = sub_Acqu.at[index, 'fish_shape']
                image_map = sub_Acqu.at[index, 'fish_map'].copy()
                image_map = correct_map(map=image_map)
                z = image_map['z']
                c = image_map['c']
                image_number = shape[z] * shape[c]
                
                image = open_image(fullpath, image_number= image_number)
                image = reshape_stack(image, image_map=image_map, im_shape=shape)
                image = image[...,color_id]
                
                if drift_correction :
                    drift = list(sub_Acqu.loc[index, ['drift_z','drift_y','drift_x']].astype(int))
                    image = shift_array(image, *drift)

                max_shape_no_channel = tuple(max_shape[:c]) + tuple(max_shape[(c+1):])
                if (image.shape != max_shape_no_channel) :
                    image = pad_to_shape(image, new_shape=max_shape_no_channel)

                image_list.append(image)
            array = np.stack(image_list)

            print(f"SCALE FOR FISH : {self.voxel_size}")

            layerdata = (
                array,
                {
                    "scale" : self.voxel_size, 
                    "name" : name, 
                    'blending' : 'additive', 
                    'colormap' : color},
                'Image'
            )

            return layerdata

        return load

@register_load_widget   
class DapiSignalLoader(NapariWidget) :
    def __init__(
            self, 
            table_dict : table_dict_type,
            voxel_size :tuple,
            **kwargs
            ):
    

        self.Gene_map = table_dict['Gene_map']

        Drift = table_dict['Drift']
        self.Drift = Drift.loc[:,['acquisition_id', 'drift_z', 'drift_y', 'drift_x']]
        
        self.Acquisition : pd.DataFrame = table_dict['Acquisition']
        self.update(list(self.Acquisition['location'].unique()))
        self.channel_indexer : int = self.Acquisition['dapi_channel'].iat[0]

        self.voxel_size = voxel_size
        super().__init__()

    def update(self, locations) :
        """
        Updates acquisition to loc on selected locations.
        """
        self.data = pd.merge(
            self.Acquisition[self.Acquisition['location'].isin(locations)],
            self.Drift,
            on='acquisition_id'
        )
        self.data = self.data.sort_values(['location', 'full_path'])
        self.target = sorted(list(self.Gene_map['target'].unique()))

    def _create_widget(self) :
        @magicgui(
                call_button="Load dapi signal",
                drift_correction={
                    "widget_type" : "CheckBox",
                    "text" : "drift correction",
                    "value" : True,
                    },
                auto_call=False
                )
        def load_dapi(radio_button, drift_correction) -> LayerDataTuple:

            name = "dapi_signal"

            if drift_correction :
                name += "_corrected"
            else :
                name += "_drifted"

            image_list = []
            if len(self.data) > 1 :
                max_shape = np.array(list(self.data["fish_shape"]), dtype=int).max(axis=0)
            else :
                max_shape = self.data['fish_shape'].iat[0]
            max_shape = (max_shape[0],max_shape[2], max_shape[3])

            for index in tqdm(self.data.index, desc= "opening {0}".format(radio_button)) :
                full_path = self.data.at[index, "full_path"]
                shape = self.data.at[index, 'fish_shape']
                image_map = self.data.at[index, 'fish_map']
                image_number = shape[0] * shape[1]
                
                image = open_image(full_path, image_number=image_number)
                image = image.reshape(*shape)
                image = reorder_image_stack(image, image_map)
                image = image[:,:,:,self.channel_indexer]

                if drift_correction :
                    drift = list(self.data.loc[index, ['drift_z','drift_y','drift_x']].astype(int))
                    image = shift_array(image, *drift)

                if image.shape != max_shape :
                    image = pad_to_shape(image, new_shape=max_shape)

                image_list.append(image)

            array = np.stack(image_list)


            layerdata = (
                array,
                {"scale" : self.voxel_size, "name" : name, 'blending' : 'additive', 'colormap' : 'blue'},
                'Image'
            )

            return layerdata

        return load_dapi

@register_load_widget   
class BeadsSignalLoader(NapariWidget) :
    def __init__(
            self, 
            table_dict : table_dict_type,
            voxel_size :tuple,
            **kwargs,
            ):

        self.Gene_map = table_dict['Gene_map']
        self.Acquisition = table_dict['Acquisition']
        
        
        self.Drift = table_dict['Drift'].loc[table_dict['Drift']['drift_type'] == 'fish']
        self.Drift = self.Drift.loc[:,['acquisition_id', 'drift_z', 'drift_y', 'drift_x']]

        self.update(list(self.Acquisition['location'].unique()))
        self.channel_indexer = self.Acquisition['bead_channel'].iat[0]


        
        self.voxel_size = voxel_size
        super().__init__()

        if self.Acquisition['bead_channel'].isna().all() :
            self.disable_widget()

    def update(self, locations) :
        self.data = pd.merge(
            self.Acquisition[self.Acquisition['location'].isin(locations)],
            self.Drift,
            on='acquisition_id'
        )
        
        Gene_map = pd.merge(
            self.Gene_map,
            self.data.loc[:,['cycle']],
            on = 'cycle'
        )
        
        self.data = self.data.sort_values(['location', 'full_path'])
        self.Gene_map_filtered = Gene_map
        self.target = sorted(list(Gene_map['target'].unique()))

        joined_names = self.Gene_map.groupby('cycle')['target'].apply('-'.join)
        targets_names = ["{0}({1})".format(cycle_num,target) for target,cycle_num in zip(joined_names, joined_names.index)]
        self.target = list(joined_names.index)
        self.target_names = targets_names

    def _create_widget(self) :
        @magicgui(
                target = {
                    "widget_type" : "ComboBox",
                    'choices' : self.target_names,
                    "label" : 'cycle',
                    },
                drift_correction={
                    "widget_type" : "CheckBox",
                    "text" : "drift correction",
                    "value" : True,
                    },
                call_button="Load beads",
                auto_call=False
                )
        def load(target, drift_correction) -> LayerDataTuple:
            target_index = self.target_names.index(target)
            target = self.target[target_index]
            gene_data = self.Gene_map_filtered.loc[self.Gene_map_filtered['cycle'] == target].iloc[0]
            shape = self.data[self.data['cycle'] == target].iloc[0]['fish_shape']
            image_map = self.data[self.data['cycle'] == target].iloc[0]['fish_map'].copy()
            image_map = correct_map(image_map)
            cycle = gene_data['cycle']

            if target == 0 :
                color = 'green'
            elif drift_correction :
                color = 'red'
            else :
                color = 'blue'

            if drift_correction :
                name = "{0}_beads_signal_corrected".format(target)
            else :
                name = "{0}_beads_signal_drifted".format(target)

            sub_Acqu = self.data.loc[self.data['cycle'] == cycle]

            z = image_map['z']
            c = image_map['c']
            image_number = shape[z] * shape[c]
            image_list = []

            if len(sub_Acqu) > 1 :
                max_shape = np.array(list(sub_Acqu["fish_shape"]), dtype=int).max(axis=0)
            else :
                max_shape = sub_Acqu['fish_shape'].iat[0]

            max_shape = tuple(max_shape[:c]) + tuple(max_shape[(c+1):])

            for index in tqdm(sub_Acqu.index, desc="Opening beads ({0})".format(target)) :
                full_path = sub_Acqu.at[index, "full_path"]

                image = open_image(full_path, image_number=image_number)
                image = reshape_stack(image, image_map=image_map, im_shape= shape)
                image = image[...,self.channel_indexer]
                
                
                if drift_correction :
                    drift = list(sub_Acqu.loc[index, ['drift_z','drift_y','drift_x']].astype(int))
                    image = shift_array(image, *drift)
                
                if (image.shape != max_shape) :
                    image = pad_to_shape(image, new_shape=max_shape)


                image_list.append(image)
            array = np.stack(image_list)

            layerdata = (
                array,
                {"scale" : self.voxel_size, "name" : name, 'blending' : 'additive', 'colormap' : color},
                'Image'
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
            **kwargs
            ):
        
        Drift = table_dict['Drift']
        self.Drift = Drift.loc[Drift['drift_type'] == 'dapi']
        self.Acquisition = table_dict['Acquisition']

        self.update(list(self.Acquisition['location'].unique()))

        self.segmentation_fullpath = run_path + segmentation_folder_name
        self.voxel_size = voxel_size
        super().__init__()

    def update(self,locations) :
        self.data = pd.merge(
            self.Acquisition[self.Acquisition['location'].isin(locations)],
            self.Drift,
            on= 'acquisition_id',
            how= 'left', #with Drift loc on dapi, only Acquisition on cycle 0 will have a drift, other will be Na
            suffixes= ('','_drift')
        )

    def _create_widget(self) :

        @magicgui(
                call_button= "Load segmentation",
                object={
                    "widget_type" : "RadioButtons",
                    "choices" : ["nucleus","cytoplasm"],
                    "orientation" : "horizontal",
                    "value" : "nucleus",
                    "label" : " ",
                },
                drift_correction={
                    "widget_type" : "CheckBox",
                    "text" : "drift correction (nucleus)",
                    "value" : True,
                    },
                auto_call=False
        )
        def load_segmentation(object, drift_correction) -> LayerDataTuple:

            shape_fish = np.array(list(self.Acquisition['fish_shape']),dtype=int)
            shape_fish = np.max(shape_fish, axis=0)
            shape_dapi = np.array(list(self.Acquisition['dapi_shape']),dtype=int)
            shape_dapi = np.max(shape_dapi, axis=0)
            shape = np.max([shape_fish,shape_dapi],axis=0)
            z_size = shape[0]
            name = "{0}_mask".format(object)
            locations = list(self.data.sort_values('location')['location'].unique())
            masks = open_segmentation(self.segmentation_fullpath, locations , object=object, z_repeat= z_size) #masks list sorted on Acquisition['location']

            if drift_correction and object == "nucleus" :

                name += '_corrected'

                data = self.data.sort_values("location")
                sub_acq = data[~data['drift_id'].isna()]
                drift_z = list(sub_acq['drift_z'].astype(int))
                drift_y = list(sub_acq['drift_y'].astype(int))
                drift_x = list(sub_acq['drift_x'].astype(int))
                drift_list = list(zip(drift_z,drift_y,drift_x))
                assert len(drift_list) == len(masks), "Didn't find drift correction for all masks."

                for location_index, drift in tqdm(enumerate(drift_list), desc= "correcting drift", total= len(drift_list)) :
                    masks[location_index] = shift_array(masks[location_index], *drift)

            elif object == "nucleus" :
                name += '_drifted'


            layerdata = (
                masks,
                {"scale" : self.voxel_size, "name" : name, "blending" : "additive"},
                'Labels'
            )


            return layerdata
        return load_segmentation

## Analysis widgets
@register_analysis_widget
class MultichannelCluster(NapariWidget) :
    def __init__(
            self, 
            table_dict, 
            voxel_size,
            **kwargs
            ):
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

            layer_data = (clusters, 
                         {"scale" : self.voxel_size, 
                          "name" : name, 
                          'ndim' : 4, 
                          'face_color' : 'white',
                          'symbol' : 'square',
                          'features' : {'single_number' : single_number, 'target_names' : target_names}, 
                          'size' : 0.1, 
                          'blending' : 'additive'}
                          , 'Points')
            return layer_data
            
        
        return multichannel_DBSCAN

@register_analysis_widget
class SpotCountMapper(NapariWidget) :
    def __init__(
            self, 
            table_dict, 
            voxel_size,
            **kwargs,
            ):
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

            layer_data = (
                spot_count_array,
                {"name" : "spot_count_map",
                 "scale" : self.voxel_size,
                 "blending" : "additive",
                 "colormap" : "inferno"
                 },
                 'Image'
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
            linked_widgets : 'list[Widget]',
            **kwargs
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