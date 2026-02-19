import os
import json

from pydantic import BaseModel, ValidationError
from .types import parameters_dict
from .pipeline_parameters import get_default_settings


def get_settings() -> parameters_dict :

    setting_path = get_settings_path()

    if os.path.isfile(setting_path) :
        return _load_settings()
    else :
        settings = get_default_settings()
        write_settings(settings)
        return settings

def _load_settings() :
    settings_path = get_settings_path()
    with open(settings_path, "r") as f:
        settings = json.load(f)
    
    try : settings = parameters_dict(**settings)

    except ValidationError as e :
        print(f"Incorrect settings, using default settings \n{e}")
        settings = get_default_settings()
    
    return settings

def get_settings_path() :
    return os.path.join(os.path.dirname(__file__) , "settings.json")

def write_settings(settings : parameters_dict) :
    if not isinstance(settings, parameters_dict) :
        raise TypeError("Expected SettingsDict type, got {}".format(type(settings)))
    else :
        settings_path = get_settings_path()
        with open(settings_path, mode="w") as f:
             json.dump(settings.dict(), f, indent=4)