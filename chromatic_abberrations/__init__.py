"""
Module to handle chromatic abberrations corrections and its calibration.
"""
from pathlib import Path

CALIBRATION_FOLDER = str(Path(__file__).parent.resolve()) + "/saved_calibrations/"

def run(*args) :
    from .launch_calibration import main
    main()