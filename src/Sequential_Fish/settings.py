import os
import json

from pydantic import ValidationError
from .types import parameters_dict
from .pipeline_parameters import get_default_settings


def get_settings(run_path : str) -> parameters_dict :

    if os.path.isfile(run_path + "/pipeline_settings") :
        return _load_settings(run_path + "/pipeline_settings")
    else :
        settings = get_default_settings()
        write_settings(settings, run_path)
        return settings

def _load_settings(run_path : str) :
    with open(run_path, "r") as f:
        settings = json.load(f)
    
    try : settings = parameters_dict(**settings)

    except ValidationError as e :
        print(f"Incorrect settings, using default settings \n{e}")
        settings = get_default_settings()
    
    return settings

def get_settings_path() :
    return os.path.join(os.path.dirname(__file__) , "settings.json")

def write_settings(settings : parameters_dict, run_path : str) :
    if not isinstance(settings, parameters_dict) :
        raise TypeError("Expected SettingsDict type, got {}".format(type(settings)))
    else :
        with open(run_path + "/pipeline_settings.json", mode="w") as f:
            json.dump(settings.model_dump(), f, indent=4)