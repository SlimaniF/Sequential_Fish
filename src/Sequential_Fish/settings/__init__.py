"""
Module handling cache. Enabling selection of runs in analysis.
"""

from .settings import get_settings
from .settings import write_settings

def run(*_) :
    from .open_status import main
    main()