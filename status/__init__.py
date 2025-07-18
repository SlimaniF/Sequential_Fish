"""
Module handling cache. Enabling selection of runs in analysis.
"""

from .gui import select_path_for_pipeline
from .gui import select_path_for_analysis

from .parameters import (
    load_pipeline_parameters,
    write_pipeline_parameters,
    get_raw_pipeline_parameters,
    exists_pipeline_parameters,
    load_analysis_parameters,
    write_analysis_parameters,
    get_raw_analysis_parameters,
    exists_analysis_parameters,
    )

from .cache import read_cache

from .parameters import load_pipeline_parameters

def run(*args) :
    from .open_status import main
    main()