"""
Module handling cache. Enabling selection of runs in analysis.
"""

from .settings import SETTINGS_NAMES
from .settings import ALLOWED_SETTINGS
from .settings import get_settings
from .settings import write_settings
from .settings import modify_settings as main

def run(run_path, settings_name : SETTINGS_NAMES) :
    main(run_path, settings_name=settings_name)