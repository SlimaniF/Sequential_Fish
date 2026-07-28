"""
Entry point to launch Sequential Fish viewer
"""
from .main import main

import logging

from .load import initiate_load_widgets
from . analysis import initiate_analysis_widgets
from . locations import initiate_location_widgets
from .organoids import initiate_organoid_wizards
from . thresholds import initiate_thresholds_widgets
from .thresholds import ThreholdsFileEditor
from .utils import update_layer_from_LayerDataTuple


def run(run_path : str) :
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler()]
    )
    logging.info("Starting viewer module")
    main(run_path)