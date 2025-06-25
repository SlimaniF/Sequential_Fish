"""
Module to handle chromatic abberrations corrections and calibration.
"""
from pathlib import Path

CALIBRATION_FOLDER = str(Path(__file__).parent.resolve()) + "/chromatic_abberrations/saved_calibrations/"