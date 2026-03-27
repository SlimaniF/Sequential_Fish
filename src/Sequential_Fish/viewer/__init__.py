"""
Entry point to launch Sequential Fish viewer
"""
from .main import main

from .load import initiate_load_widgets
from . analysis import initiate_analysis_widgets
from . locations import initiate_location_widgets
from .organoids import initiate_organoid_wizards
from . thresholds import initiate_thresholds_widgets
from .thresholds import ThreholdsFileEditor


def run(run_path : str) :
    main(run_path)