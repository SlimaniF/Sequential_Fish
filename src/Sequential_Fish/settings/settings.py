import os
import sys
import json
from pathlib import Path
import warnings
from typing import Literal, cast
from PyQt5.QtWidgets import QApplication
from pydantic import ValidationError

from ..customtypes import PipelineParameters, AnalysisParameters
from . gui import ParametersModifier

ALLOWED_SETTINGS = ['pipeline', 'analysis']
SETTINGS_NAMES = Literal['pipeline', 'analysis']

def get_settings(
    run_path : str | Path,
    settings_name : Literal['pipeline','analysis'] = 'pipeline'
    ) -> PipelineParameters | AnalysisParameters :

    match settings_name :
        case 'pipeline' :
            model = PipelineParameters
        case 'analysis' :
            model = AnalysisParameters
        case _ :
            raise ValueError(f"Incorrect value for settings :{settings_name}; expected 'pipeline' or 'analysis'.")

    settings_path = os.path.join(run_path , f"{settings_name}_settings.json")
    if os.path.isfile(settings_path) :
        with open(settings_path, mode="r") as setting_file:
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
    run_path : str | Path
    ) :
    
    if not isinstance(settings, PipelineParameters) :
        raise TypeError("Expected SettingsDict type, got {}".format(type(settings)))
    filename = settings.get_filename()
    with open(os.path.join(run_path , f"{filename}"), mode="w") as f:
        json.dump(settings.model_dump(), f, indent=4)


def modify_settings(
    run_path : str | Path,
    settings_name  : Literal["pipeline","analysis"]
    ) :

    match settings_name :
        case "pipeline" :
            data_model = PipelineParameters
        case "analysis" :
            data_model = AnalysisParameters

    settings = get_settings(run_path, settings_name)

    app = QApplication(sys.argv[:1])
    gui_modifier = ParametersModifier(data_model, **dict(settings))

    if gui_modifier.exec() :
        new_parameters_dict = gui_modifier.get_parameters()
        try : 
            settings = data_model.model_validate(new_parameters_dict)
        except ValidationError as e :
            print("Uncorrect parameter set.")
            print(e)

            wrong_attributes = [error_dict["loc"][0] for error_dict in e.errors()]
            wrong_attributes = cast(list[str], wrong_attributes)
            print(f"Parameters {wrong_attributes} were wrongly set, reverting to previous values.")
            for att in wrong_attributes :
                model_attribute = getattr(data_model.model_fields, att)
                if hasattr(model_attribute, "default_factory") :
                    new_parameters_dict[att] = model_attribute.default_factory()
                elif hasattr(model_attribute, "default") :
                    new_parameters_dict[att] = model_attribute.default
                else :
                    raise ValidationError(f"{att} has no default value, and was set wrongly ({new_parameters_dict[att]}) could not revert to default.") from e                
                
            settings = AnalysisParameters(**new_parameters_dict)
            write_settings(
                settings= settings,
                run_path=run_path
            )

        else :
            write_settings(
                settings= settings,
                run_path=run_path
            )