import os, platform
from typing import cast, Sequence
import napari
import pandas as pd

from numpy.random import shuffle

from .types import NapariWidget
from ..settings import PipelineParameters
from ..customtypes import table_dict_type
from .analysis import initiate_analysis_widgets
from .load import initiate_load_widgets
from .locations import initiate_location_widgets
from .thresholds import initiate_thresholds_widgets, ThreholdsFileEditor
from .segmentation import initiate_segmentation_widgets
from .organoids import initiate_organoid_wizards, autodetect_organoids
from .chromatic_abberations import initiate_chromatic_widgets
from magicgui.widgets import Container, Widget

from .utils import get_colors_list, _get_blue_colors, _get_green_colors, _get_orange_colors, _get_red_colors, _get_yellow_colors, _get_pink_colors, _get_purple_colors
from ..settings import get_settings


def main(run_path) :
    
    if run_path is None : quit()
    system_type = platform.system()
    if system_type == "Linux" :
        try :
            os.environ["QT_QPA_PLATFORM"] = "xcb"
        except Exception :
            pass
    
    TABLES = ['Acquisition', 'Detection', 'Spots', 'Clusters', 'Drift', 'Cell', 'Gene_map']
    for table in TABLES.copy() : 
        if not os.path.isfile(run_path + '/result_tables/' + table + '.feather') :
            if table == "Acquisition" or table == "Gene_map" :
                raise FileNotFoundError("Acquisiton or Gene_map was not found in result tables. Make sure at least pipeline input was run successfuly.")
            else :
                TABLES.remove(table)

    tables_dict = cast(
        table_dict_type,
        {
        table : pd.read_feather(run_path + '/result_tables/' + table + '.feather')  for table in TABLES
        })

    settings = cast(PipelineParameters, get_settings(run_path, settings_name="pipeline"))
    voxel_size = settings.VOXEL_SIZE

    #Init viewer
    Viewer = napari.Viewer(title=os.path.basename(run_path))
    linked_widgets = []
    color_table = create_color_table(tables_dict)

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

    #Segmentation tab
    segmentation_container, new_linked_widgets = create_segmentation_tab(viewer=Viewer, voxel_size=voxel_size)
    linked_widgets.extend(new_linked_widgets)
    Viewer.window.add_dock_widget(
        segmentation_container,
        name="Segmentation",
        area="right",
        add_vertical_stretch=True,
        tabify=True
    )

    #Threholds tab
    full_threshold_container = create_thresholds_tab(
        viewer=Viewer,
        voxel_size= settings.VOXEL_SIZE,
        spot_size= settings.SPOT_SIZE,
        run_path= run_path,
        map_filename=settings.MAP_FILENAME,
        gene_number=len(settings.GENES_NAMES_KEY)
    )
    Viewer.window.add_dock_widget(
        full_threshold_container,
        name="Thresholds",
        area="right",
        add_vertical_stretch=True,
        tabify=True
    )

    #Chromatic tab
    if not settings.WAVELENGTH_LIST is None :
        chromatic_container, new_linked_widgets = create_chromatic_tab(
            viewer=Viewer,
            wavelength_list=settings.WAVELENGTH_LIST,
            voxel_size=voxel_size
        )
        linked_widgets.extend(new_linked_widgets)
        Viewer.window.add_dock_widget(
            chromatic_container,
            name= "Chromatic aberrations",
            area="right",
            add_vertical_stretch=True,
            tabify=True
        )

    #Organoids
    organoids_locations = autodetect_organoids(run_path, Acquisition=tables_dict['Acquisition'])
    if organoids_locations is None :
        print("No organoids location found.")
    else :
        print("Organoids locations found")
        organoids_wizards, organoids_location_linked_widgets = initiate_organoid_wizards(viewer= Viewer, organoids_locations=organoids_locations)
        linked_widgets.extend(cast(list[NapariWidget], organoids_location_linked_widgets))

    # Location tab
    location_widgets = initiate_location_widgets(
        tables_dict=tables_dict, 
        Viewer=Viewer,
        linked_widgets=linked_widgets 
        )
    location_container = Container(
        widgets= location_widgets,
        labels=False
    )


    #TODO Update
    # Viewer.window.add_dock_widget(
    #     analysis_container, 
    #     name='Analysis', 
    #     area='right', 
    #     add_vertical_stretch=True, 
    #     tabify=True
    #     )

    #Load tab
    load_data_container, new_linked_widgets = create_load_tab(
        tables_dict=tables_dict,
        voxel_size=voxel_size,
        color_table=color_table,
        run_path=run_path,
        viewer=Viewer
    )
    linked_widgets.extend(new_linked_widgets)
    Viewer.window.add_dock_widget(
        load_data_container, 
        name='Data', 
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
    Viewer.canvas.overlays.scale_bar.visible = True

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


def create_load_tab(
    tables_dict : table_dict_type,
    voxel_size : tuple,
    color_table : dict,
    run_path : str,
    viewer : napari.Viewer
) :
    #Loading tab
    load_data_widgets, linked_widgets = initiate_load_widgets(
        voxel_size=voxel_size,
        table_dict=tables_dict,
        color_table=color_table,
        run_path=run_path,
        viewer = viewer,
    )
    load_data_container = Container(
        widgets= cast(Sequence[Widget], load_data_widgets),
        labels=False
    )

    return load_data_container, linked_widgets

def create_segmentation_tab(
    viewer : napari.Viewer,
    voxel_size : tuple,
    ) :

    widget_list, linked_widgets =  initiate_segmentation_widgets(viewer, voxel_size=voxel_size)
    segmentation_container = Container(
        widgets= widget_list,
        labels=False
    )

    return segmentation_container, linked_widgets


def create_thresholds_tab(
    viewer : napari.Viewer,
    voxel_size : tuple,
    spot_size : tuple,
    run_path : str,
    map_filename : str,
    gene_number : int,
    ) :

    thresholds_widgets = initiate_thresholds_widgets(
        viewer= viewer,
        default_threshold = 0,
        default_spot_size = spot_size,
        voxel_size = voxel_size,
        )

    cyclefile_path= os.path.join(run_path, map_filename)
    threshold_editor = ThreholdsFileEditor(cyclefile_path, color_number = gene_number)

    thresholds_container = Container(
        widgets= thresholds_widgets,
        labels=True
    )
    threshold_editor_container = Container(
        widgets= threshold_editor.widget,
        labels=False
    )

    full_threshold_container = Container(
        widgets = [thresholds_container, threshold_editor_container],
        labels=False
    )

    return full_threshold_container

def create_chromatic_tab(
    viewer : napari.Viewer,
    wavelength_list : list[int],
    voxel_size : tuple,
    ) :

    chromatic_widgets, linked_widgets = initiate_chromatic_widgets(viewer=viewer, wavelength_list=wavelength_list, voxel_size=voxel_size)
    chromatic_container = Container(
        widgets= chromatic_widgets,
        labels=False
    )

    return chromatic_container, linked_widgets