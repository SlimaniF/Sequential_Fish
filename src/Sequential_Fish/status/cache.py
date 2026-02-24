import os, json
from typing import Dict

def read_cache() -> Dict[str,str]:
    script_path = os.path.dirname(__file__)
    if os.path.isfile(script_path + "/run_cache.json") :
        with open(script_path + "/run_cache.json","r") as run_cache:
            try :
                return json.load(run_cache)
            except Exception as e :
                print("Exception raised when loading cache. Ignoring cache")
                print(e)
                return dict()
    else :
        return dict()

def write_cache(path_map : dict) :
    script_path = os.path.dirname(__file__)
    with open(script_path + "/run_cache.json", 'w') as cache_file:
        json.dump(path_map, cache_file, indent=2)