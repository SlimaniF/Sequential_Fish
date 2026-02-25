import os
import json
from typing import Literal

from pydantic import ValidationError
from ..customtypes import PipelineParameters, AnalysisParameters



def get_settings(
    run_path : str,
    settings_name : Literal['pipeline','analysis'] = 'pipeline'
    ) -> PipelineParameters | AnalysisParameters :

    match settings_name :
        case 'pipeline' :
            model = PipelineParameters
        case 'analysis' :
            model = AnalysisParameters
        case _ :
            raise ValueError(f"Incorrect value for settings :{settings_name}; expected 'pipeline' or 'analysis'.")

    if os.path.isfile(run_path + f"/{settings_name}_settings.json") :
        print("open settings")
        with open(run_path + f"/{settings_name}_settings.json", mode="r") as setting_file:
            try :
                settings = model.model_validate_json(setting_file.read())
            except ValidationError as validation_error:
                for error in validation_error.errors() :
                    print(f"{error["loc"]} :\n{error["msg"]}")
                print(f"Uncorrect settings found in {settings_name}_settings.json")
                raise ValidationError() from validation_error
        
        return settings

    else :
        print("default settings")
        settings = model.from_default_parameters()
        write_settings(settings, run_path)
        return settings

def get_settings_path() :
    return os.path.join(os.path.dirname(__file__) , "settings.json")

def write_settings(
    settings : PipelineParameters | AnalysisParameters, 
    run_path : str
    ) :
    
    if not isinstance(settings, PipelineParameters) :
        raise TypeError("Expected SettingsDict type, got {}".format(type(settings)))
    filename = settings.get_filename()
    with open(run_path + f"/{filename}", mode="w") as f:
        json.dump(settings.model_dump(), f, indent=4)