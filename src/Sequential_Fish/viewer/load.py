"""
Napari widgets used to load information in napari ()
"""

import os
from typing import cast

from napari import Viewer
import numpy as np
import pandas as pd

from napari.types import LayerDataTuple
from napari.qt.threading import thread_worker, create_worker
from magicgui import magicgui

from ..tools.utils import open_all_locations_one_cycle, safe_merge_no_duplicates
from .utils import open_segmentation, update_layer_from_LayerDataTuple
from .types import NapariWidget
from ..customtypes import table_dict_type


from smfishtools.preprocessing.alignement import shift_array


class LoadWidget(NapariWidget) :
    """
    Subclass of NapariWidget aimed for widgets returning a LayerDataTuple, it wraps the loading process into a worker thread to not pause the main GUI thread.
    """

    def __init__(self, *, viewer : Viewer):
        super().__init__()
        self.viewer = viewer
        self._wrap_with_threading()

    def _wrap_with_threading(self):
        """
        Wraps the widget's function with napari's threading.
        """
        if hasattr(self.widget, 'native'):
            # Get the underlying function from the magicgui widget
            original_func = self.widget._function
            def threaded_func(*args, **kwargs):
                worker = create_worker(original_func, *args, **kwargs)
                worker.returned.connect(lambda result: update_layer_from_LayerDataTuple(self.viewer, result))
                worker.start()
                return worker
            # Replace the widget's function with the threaded version
            self.widget._function = threaded_func

_LOAD_WIDGETS : 'list[NapariWidget]' = []
def register_load_widget(cls) :
    _LOAD_WIDGETS.append(cls)
    return cls

def initiate_load_widgets(
        table_dict : table_dict_type, 
        voxel_size : tuple,
        color_table : dict,
        run_path : str,
        viewer : Viewer,
) -> tuple[list[NapariWidget], list[NapariWidget]]:
    widget_list = []
    linked_widgets = []
    for cls in _LOAD_WIDGETS :
        instance = cls(table_dict= table_dict, voxel_size=voxel_size, color_table=color_table, run_path=run_path, viewer=viewer)
        if hasattr(instance,"enabled") :
            if instance.enabled :
                widget_list.extend(instance.get_widgets())
                linked_widgets.append(instance)
        else :
            widget_list.extend(instance.get_widgets())
            linked_widgets.append(instance)

    return widget_list, linked_widgets


@register_load_widget
class SpotsLoader(LoadWidget) :
    """
    Allow user to load detected spots as layer points
    """
    def __init__(    
        self, 
        table_dict : table_dict_type,
        voxel_size :tuple, 
        color_table,
        viewer : Viewer,
        **_
    ) :
        self.viewer = viewer
        self.enabled=True
        self.Spots = table_dict.get('Spots')
        self.Detection = table_dict.get('Detection')

        if self.Spots is None or self.Detection is None :
            print("Disabling spots loader (Detection results not found)")
            self.enabled=False

            return None

        self.Acquisition = table_dict['Acquisition'].loc[:,['acquisition_id','location','cycle']]
        self.location_list = self.Acquisition['location'].unique().tolist()
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
        super().__init__(viewer=viewer)

    def update(self, locations) :

        data = safe_merge_no_duplicates(
            self.Spots,
            self.Detection,
            on= 'detection_id',
            keys= ["location", "cycle", "color_id","acquisition_id", "voxel_size", "spot_size"]
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
        self.location_list = locations

        assert not any(data['target'].isna()), "Missing values for `target` in Spots. Merge is incomplete."
        self.data = data
        data.groupby(["target","location"])["spot_id"].count().to_excel("/media/SSD_floricslimani/Fish_seq/test.xlsx")
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
        def load(target, population, drift_correction) :
            
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
            voxel_sizes = []
            spots_sizes = []
            for location_index, location  in enumerate(self.location_list) :

                spot_data = sub_Detec.loc[sub_Detec['location'] == location]
                C = [location_index] * len(spot_data)
                Z = spot_data[z_indexer]
                Y = spot_data[y_indexer]
                X = spot_data[x_indexer]
                voxel_sizes.extend(spot_data["voxel_size"].to_list())
                spots_sizes.extend(spot_data["spot_size"].to_list())

                if len(spot_data) == 0 :
                    spots = np.empty(shape=(0,4),dtype=int)
                else :
                    spots = np.array(
                    list(zip(C,Z,Y,X)),
                    dtype=int,
                    )

                spots_array = np.concatenate([spots_array, spots])
            spots_sizes = pd.Series(spots_sizes).apply(tuple)
            voxel_sizes = pd.Series(voxel_sizes).apply(tuple)
            layerdata = cast(LayerDataTuple, (spots_array, 
                         {
                             "scale" : self.voxel_size,
                             "size" : 10, 
                             "name" : name, 
                             'ndim' : 4, 
                             'face_color' : '#0000' ,
                             'border_color' : color, 
                             'blending' : 'additive',
                             'symbol' : symbol,
                             'features' : {"voxel size" : voxel_sizes, "spot size" : spots_sizes}
                             },
                        'Points'))
            return layerdata
        return load

@register_load_widget    
class ClustersLoader(LoadWidget) :
    """
    Allow user to load detected cluster as Points layer.
    """

    def __init__(
            self, 
            table_dict : table_dict_type,
            voxel_size :tuple, 
            color_table,
            viewer : Viewer,
            **_
            ):
        
        self.viewer = viewer
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
        self.location_list = self.Acquisition['location'].unique().tolist()
        self.update(list(self.Acquisition['location'].unique()))

        self.voxel_size = voxel_size
        self.color_table = color_table
        super().__init__(viewer=viewer)

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
        self.location_list = locations
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
        def load(target, drift_correction) :
            
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
            for location_index, location  in enumerate(self.location_list) :
                spots_data = sub_data.loc[sub_data['location'] == location]
                C = [location_index] * len(spots_data)
                Z = spots_data[z_indexer]
                Y = spots_data[y_indexer]
                X = spots_data[x_indexer]

                if len(spots_data) == 0 :
                    spots = np.empty(shape=(0,4),dtype=int)
                else :
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
class SignalLoader(LoadWidget) :
    def __init__(
            self, 
            table_dict : table_dict_type,
            voxel_size :tuple, 
            color_table : pd.DataFrame,
            viewer : Viewer,
            **_
            ):
        
        self.viewer = viewer

        #Table
        self.Gene_map = table_dict['Gene_map']

        Drift = table_dict.get('Drift')
        if Drift is None :
            self.Drift = None
        else :
            self.Drift = Drift.loc[:,['acquisition_id', 'drift_z', 'drift_y', 'drift_x']]
        
        self.Acquisition = table_dict['Acquisition'].set_index(["location","cycle"], verify_integrity=True)
        
        self.update(list(self.Acquisition.index.get_level_values(0).unique()))

        self.color_table = color_table.set_index("target",verify_integrity=True)
        self.voxel_size = voxel_size
        self.has_beads = not cast(bool,self.Acquisition['bead_channel'].isna().all())
        super().__init__(viewer=viewer)

    def update(self, locations) :

        if not self.Drift is None :
            self.data = pd.merge(
                self.Acquisition[self.Acquisition.index.get_level_values(0).isin(locations)].reset_index(drop=False),
                self.Drift,
                on='acquisition_id'
            )
        else :
            self.data = self.Acquisition[self.Acquisition.index.get_level_values(0).isin(locations)].reset_index(drop=False)
        self.data = self.data.set_index(['location',"cycle"]).sort_index()
        self.target = list(self.Gene_map['target'].unique())

    def _create_widget(self) :
        @magicgui(
                target = {'choices' : self.target},
                drift_correction={
                    "widget_type" : "CheckBox",
                    "text" : "drift correction",
                    "value" : not self.Drift is None,
                    "enabled" : not self.Drift is None
                    },
                call_button="Load signal",
                signal_type = {
                    "widget_type" : "RadioButtons",
                    "choices" : ['fish','dapi','beads'] if self.has_beads else ['fish', 'dapi'],
                    "value" : "fish"
                },
                auto_call=False
                )
        def load(target, drift_correction, signal_type):
            data = self.Gene_map.loc[self.Gene_map['target'] == target].iloc[0]
            color = self.color_table.at[target, "colormaps"]
            cycle, color_id = data['cycle'], int(data['color_id'])

            if signal_type == "fish" :
                channel_index = color_id
                name = f"{target}_fish"
            elif signal_type == "dapi" :
                channel_index = self.data['dapi_channel'].iat[0]
                name = f"dapi_signal_cycle{cycle}"
            elif signal_type == "beads" :
                channel_index = self.data['bead_channel'].iat[0]
                name = f"beads_signal_cycle{cycle}"
            else :
                raise NotImplementedError("Unimplemented choice")


            if drift_correction :
                name += "_corrected"
            else :
                name += "_signal_drifted"


            array = open_all_locations_one_cycle(
                self.data.reset_index(drop=False),
                cycle=cycle,
            )
            array = array[..., channel_index]

            if drift_correction :
                location_list = self.data.index.get_level_values(0).unique().to_list()
                location_list.sort()
                assert len(location_list) == len(array)

                location_index=0
                for location, stack in zip(location_list, array) :
                    drift = self.data.loc[(location,cycle),["drift_z","drift_y","drift_x"]].astype(int)
                    array[location_index] = shift_array(stack, *drift)
                    location_index +=1

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
class SegmentationLoader(LoadWidget) :
    
    def __init__(
            self,
            run_path : str,
            voxel_size :tuple,
            table_dict : table_dict_type, 
            viewer : Viewer,
            segmentation_folder_name:str = "/segmentation/",
            **_
            ):
        
        self.enabled=True
        self.viewer = viewer

        self.Acquisition = table_dict['Acquisition']

        self.update(list(self.Acquisition['location'].unique()))

        self.segmentation_fullpath = run_path + segmentation_folder_name

        if not os.path.isdir(self.segmentation_fullpath) or len(os.listdir(self.segmentation_fullpath)) == 0 :
            print("Disabling segmentation masks (Segmentation results not found)")
            self.enabled=False

        self.voxel_size = voxel_size
        super().__init__(viewer=viewer)

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
        def load_segmentation(target):

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