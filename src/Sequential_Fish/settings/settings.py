import os
import json
import warnings
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
        with open(run_path + f"/{settings_name}_settings.json", mode="r") as setting_file:
            try :
                saved_settings = setting_file.read()
                settings = model.model_validate_json(saved_settings)
                saved_settings = json.loads(saved_settings)

            except ValidationError as validation_error:
                for error in validation_error.errors() :
                    print(f"{error["loc"]} :\n{error["msg"]}")
                raise ValidationError(f"Uncorrect settings found in {settings_name}_settings.json") from validation_error

            filled_settings = set(model.model_fields.keys()) - set(saved_settings.keys()) #True if any missing value were added by default
            if filled_settings :
                filled_settings_display = {setting_name : getattr(settings,setting_name) for setting_name in filled_settings}
                warnings.warn(f"Saved parameters were found but some values are missing. Autofilling and saving with default values :\n {filled_settings_display}")
                write_settings(settings, run_path)

        return settings

    else :
        settings = model.from_default_parameters()
        write_settings(settings, run_path)
        return settings


def write_settings(
    settings : PipelineParameters | AnalysisParameters, 
    run_path : str
    ) :
    
    if not isinstance(settings, PipelineParameters) :
        raise TypeError("Expected SettingsDict type, got {}".format(type(settings)))
    filename = settings.get_filename()
    with open(run_path + f"/{filename}", mode="w") as f:
        json.dump(settings.model_dump(), f, indent=4)