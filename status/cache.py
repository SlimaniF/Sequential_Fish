import os, json
from typing import Dict

def read_cache() -> Dict[str,str]:
    script_path = os.path.dirname(__file__)
    if os.path.isfile(script_path + "/run_cache.json") :
        with open(script_path + "/run_cache.json","r") as run_cache:
            return json.load(run_cache)
    else :
        return dict()