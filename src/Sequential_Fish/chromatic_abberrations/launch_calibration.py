"""
This script launch napari to allow user to calibrate chromatic abberration correction using fluorescent beads.
"""

from magicgui import widgets
from .widgets import initiate_all_calibration_widgets
from ..viewer.chromatic_abberations import SignalCorrector
from ..viewer.load import SignalLoader, SpotsLoader
from ..settings import get_settings, PipelineParameters
from ..customtypes import table_dict_type

import napari,os, platform
import pandas as pd
from typing import cast

def main(run_path) :
    system_type = platform.system()
    if system_type == "Linux" :
        try :
            os.environ["QT_QPA_PLATFORM"] = "xcb"
        except Exception :
            pass
    Viewer = napari.Viewer(title= "SequentialFish - Chromatic abberration calibration")

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

    signal_loader = SignalLoader(
        table_dict=tables_dict,
        voxel_size=voxel_size,
        viewer=Viewer
    )
    spots_loader = SpotsLoader(
        table_dict=tables_dict,
        voxel_size=voxel_size,
        viewer=Viewer
    )

    location_list = tables_dict['Acquisition']["location"].unique().tolist()
    signal_loader.update(location_list[:1])
    spots_loader.update(location_list[:1])

    signal_corrector  = SignalCorrector(viewer=Viewer, voxel_size=(200,97,97), wavelength_list=[640,555])

    calibration_widgets = [signal_loader.widget, spots_loader.widget]
    calibration_widgets.extend(initiate_all_calibration_widgets())
    calibration_widgets.append(signal_corrector.widget)

    right_container = widgets.Container(widgets=calibration_widgets, labels=False)
    Viewer.window.add_dock_widget(right_container, name='Calibration tools', area='right')
    napari.run()


if __name__ == "__main__" :
    print("This script can be launched using command 'python -m Sequential_Fish calibration {run_path}'")