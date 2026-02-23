import os
import json

from pydantic import ValidationError
from .types import pipeline_parameters


def get_settings(run_path : str) -> pipeline_parameters :

    if os.path.isfile(run_path + "/pipeline_settings") :
        return _load_settings(run_path + "/pipeline_settings")
    else :
        settings = pipeline_parameters.from_default_parameters()
        write_settings(settings, run_path)
        return settings

def _load_settings(run_path : str) :
    with open(run_path, "r") as f:
        settings = json.load(f)
    
    try : settings = pipeline_parameters(**settings)

    except ValidationError as e :
        print(f"Incorrect settings, using default settings \n{e}")
        settings = pipeline_parameters.from_default_parameters()
    
    return settings

def get_settings_path() :
    return os.path.join(os.path.dirname(__file__) , "settings.json")

def write_settings(settings : pipeline_parameters, run_path : str) :
    if not isinstance(settings, pipeline_parameters) :
        raise TypeError("Expected SettingsDict type, got {}".format(type(settings)))
    else :
        with open(run_path + "/pipeline_settings.json", mode="w") as f:
            json.dump(settings.model_dump(), f, indent=4)