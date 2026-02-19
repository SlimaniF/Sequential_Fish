import os

__version__ = 0.1
__script_dir__ = os.path.dirname(os.path.abspath(__file__))
__run_cache_path__ = os.path.join(__script_dir__, "run_saves", "Run_cache.feather")