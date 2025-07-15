import os, json
from typing import Dict

def read_cache() -> Dict[str,str]:
    script_path = __file__
    if os.path.isfile(script_path + "/run_cache.json") :
        return json.load(script_path + "/run_cache.json")
    else :
        return dict()