
import napari
import pandas as pd
import os

from numpy.random import shuffle
<<<<<<< HEAD
from .widgets import fish_container
from .widgets import detection_container
from .widgets import dapi_container
from .widgets import segmentation_container
from .widgets import locations_container
from .widgets import analysis_container
from ..types import table_dict_type
from magicgui import widgets as wi
from typing import cast
=======
from .widgets import initiate_analysis_widgets
from .widgets import initiate_load_widgets
from .widgets import initiate_location_widgets
from magicgui.widgets import Container
>>>>>>> Modification_drift_correction

from ..status.gui import select_path_for_analysis

from smfishtools.plot.utils import get_colors_list, _get_blue_colors, _get_green_colors, _get_orange_colors, _get_red_colors, _get_yellow_colors, _get_pink_colors, _get_purple_colors


def main() :
    
    run_path = select_path_for_analysis()
    if run_path is None : quit()
    
    TABLES = ['Acquisition', 'Detection', 'Spots', 'Clusters', 'Drift', 'Cell', 'Gene_map']
    tables_dict = cast(
        table_dict_type,
        {
        table : pd.read_feather(run_path + '/result_tables/' + table + '.feather')  for table in TABLES
<<<<<<< HEAD
    })
=======
    }
    voxel_size = tables_dict['Detection'].at[0,'voxel_size']
>>>>>>> Modification_drift_correction

    #Init viewer
    Viewer = napari.Viewer(title=os.path.basename(run_path))
    color_table = create_color_table(tables_dict)

    #Loading tab
    load_data_widgets = initiate_load_widgets(
        voxel_size=voxel_size,
        table_dict=tables_dict,
        color_table=color_table,
        run_path=run_path,
    )
    load_data_container = Container(
        widgets=load_data_widgets,
        labels=False
    )

    #Analysis tab
    analysis_widgets = initiate_analysis_widgets(
        table_dict = tables_dict, 
        voxel_size = voxel_size,
        color_table=color_table,
        )
    analysis_container = Container(
        widgets= analysis_widgets,
        labels=False
    )
    
    # Location tab
    linked_widgets = []
    linked_widgets.extend(load_data_widgets)
    linked_widgets.extend(analysis_widgets)
    location_widgets = initiate_location_widgets(
        tables_dict=tables_dict, 
        Viewer=Viewer,
        linked_widgets=linked_widgets 
        )
    location_container = Container(
        widgets= location_widgets,
        labels=False
    )


    #Docking tabs
    Viewer.window.add_dock_widget(
        load_data_container, 
        name='Data', 
        area='right', 
        add_vertical_stretch=True, 
        tabify=True
        )

    Viewer.window.add_dock_widget(
        analysis_container, 
        name='Analysis', 
        area='right', 
        add_vertical_stretch=True, 
        tabify=True
        )

    Viewer.window.add_dock_widget(
        location_container, 
        name='Locations', 
        area='left', 
        add_vertical_stretch=True,
        )

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