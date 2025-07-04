"""
Module to handle chromatic abberrations corrections and its calibration.
"""

from .calibration import load_calibration

from .correction import apply_polynomial_transform_3d_spots
from .correction import apply_polynomial_transform_3d_to_signal
from .correction import get_polynomial_features

from .constant import CALIBRATION_FOLDER

def run(*args) :
    from .launch_calibration import main
    main()