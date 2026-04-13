"""
This script launch napari to allow user to calibrate chromatic abberration correction using fluorescent beads.
"""

import napari,os, platform
from magicgui import widgets
from .widgets import initiate_all_calibration_widgets
from ..viewer.chromatic_abberations import SignalCorrector

def main() :
    system_type = platform.system()
    if system_type == "Linux" :
        try :
            os.environ["QT_QPA_PLATFORM"] = "xcb"
        except Exception :
            pass
    Viewer = napari.Viewer(title= "SequentialFish - Chromatic abberration calibration")

    signal_corrector  = SignalCorrector(viewer=Viewer, voxel_size=(200,97,97), wavelength_list=[640,555])

    calibration_widgets = initiate_all_calibration_widgets()
    calibration_widgets.append(signal_corrector.widget)

    right_container = widgets.Container(widgets=calibration_widgets, labels=False)
    Viewer.window.add_dock_widget(right_container, name='Calibration tools', area='right')
    napari.run()


if __name__ == "__main__" :
    os.environ["QT_QPA_PLATFORM"] = "xcb"
    main()
    print("calibration closed")