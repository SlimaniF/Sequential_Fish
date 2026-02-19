
import napari
import pandas as pd
import os

from numpy.random import shuffle
from .widgets import fish_container
from .widgets import detection_container
from .widgets import dapi_container
from .widgets import segmentation_container
from .widgets import locations_container
from .widgets import analysis_container
from ..types import table_dict_type
from magicgui import widgets as wi
from typing import cast

from ..run_saves.gui import select_path

from smfishtools.plot.utils import get_colors_list, _get_blue_colors, _get_green_colors, _get_orange_colors, _get_red_colors, _get_yellow_colors, _get_pink_colors, _get_purple_colors


VOXEL_SIZE = (200,97,97)


def main() :
    
    run_path = select_path()
    if run_path is None : quit()
    
    TABLES = ['Acquisition', 'Detection', 'Spots', 'Clusters', 'Drift', 'Cell', 'Gene_map']
    tables_dict = cast(
        table_dict_type,
        {
        table : pd.read_feather(run_path + '/result_tables/' + table + '.feather')  for table in TABLES
    })

    #Init viewer
    Viewer = napari.Viewer(title=os.path.basename(run_path))
    color_table = create_color_table(tables_dict)

    #Loading panel
    fish_buttons, fish_widgets = fish_container(VOXEL_SIZE, tables_dict, color_table)
    dapi_buttons, dapi_widgets = dapi_container(VOXEL_SIZE, tables_dict)
    segmentation_buttons, segmentation_widget = segmentation_container(run_path, tables_dict, VOXEL_SIZE)
    detection_buttons, detection_widgets = detection_container(VOXEL_SIZE, tables_dict, color_table=color_table)
    
    load_data_tab = wi.Container(widgets=[fish_buttons, dapi_buttons, segmentation_buttons, detection_buttons], labels=False, layout='vertical', name= 'Load data')

    Viewer.window.add_dock_widget(
        load_data_tab, 
        name='Data', 
        area='right', 
        add_vertical_stretch=True, 
        tabify=True
        )

    #Analysis panel
    multichannel_DBSCAN_container, multichannel_DBSCAN_instance = analysis_container(tables_dict, VOXEL_SIZE)
    Viewer.window.add_dock_widget(multichannel_DBSCAN_container, name='Analysis', area='right', add_vertical_stretch=True, tabify=True)

    #Location panel
    location_table = locations_container(tables_dict, Viewer, *detection_widgets, *fish_widgets, *dapi_widgets, *segmentation_widget, *multichannel_DBSCAN_instance)
    Viewer.window.add_dock_widget(location_table, name='Locations', area='left', add_vertical_stretch=True,)


    #Scale bar
    Viewer.scale_bar.visible = True
    Viewer.scale_bar.unit = 'nm'

    napari.run()


def create_color_table(tables_dict) :
    color_table = tables_dict['Gene_map'].loc[:,['map_id','target']]
    target_number = len(color_table)
    colors = get_colors_list(target_number, remove_black=True, remove_grey=True, remove_brown=True)
    shuffle(colors)
    color_table['color'] = colors


    colormaps_dict = {
        "blue" : _get_blue_colors (),
        "red" : _get_red_colors(),
        "bop orange" : _get_orange_colors(),
        "magenta" : _get_pink_colors(),
        "green" : _get_green_colors(),
        "yellow" : _get_yellow_colors(),
        "bop purple" : _get_purple_colors() 
    }

    colormaps = []
    for color in color_table['color'] :
        for colormap, color_list in colormaps_dict.items() :
            if color in color_list : 
                colormaps.append(colormap)
                break
        
    color_table['colormaps'] = colormaps

    return color_table

if __name__ == "__main__" :
    main()