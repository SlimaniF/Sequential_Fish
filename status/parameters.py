import os
from ..customtypes import PipelineParameters, AnalysisParameters
import Sequential_Fish.default_pipeline_parameters as pipeline_parameters
import Sequential_Fish.default_analysis_parameters as analysis_parameters


_PIPELINE_PARAMETERS_FILENAME = "/pipeline_parameters.json"
_ANALYSIS_PARAMETERS_FILENAME = "/analysis_parameters.json"

def get_raw_pipeline_parameters() :
    run_parameters = _get_defined_variable(pipeline_parameters)
    return PipelineParameters(**run_parameters)


def write_pipeline_parameters(run_path, parameters : PipelineParameters) :
    if not isinstance(parameters, PipelineParameters) : raise TypeError("Invalid type for parameters argument : {}; custom PipelineParameters type is mandatory".format(type(parameters)))
    if not os.path.isdir(run_path) : raise ValueError("Passed run_path is not a folder")

    with open(run_path + _PIPELINE_PARAMETERS_FILENAME, "w") as f:
        f.write(parameters.model_dump_json(indent=2))

def load_pipeline_parameters(run_path : str) : 
    with open(run_path + _PIPELINE_PARAMETERS_FILENAME,"r") as parameters_file :
        parameters = PipelineParameters.model_validate_json(parameters_file.read())

    return parameters

def exists_pipeline_parameters(run_path) :
    return os.path.isfile(run_path + _PIPELINE_PARAMETERS_FILENAME)

def get_raw_analysis_parameters() :
    run_parameters = _get_defined_variable(analysis_parameters)
    return AnalysisParameters(**run_parameters)

def write_analysis_parameters(run_path, parameters : AnalysisParameters) :
    if not isinstance(parameters, AnalysisParameters) : raise TypeError("Invalid type for parameters argument : {}; custom AnalysisParameters type is mandatory".format(type(parameters)))
    if not os.path.isdir(run_path) : raise ValueError("Passed run_path is not a folder")
    
    with open(run_path + _ANALYSIS_PARAMETERS_FILENAME, "w") as f:
            f.write(parameters.model_dump_json(indent=2))

def load_analysis_parameters(run_path) :
    with open(run_path + _ANALYSIS_PARAMETERS_FILENAME,"r") as parameters_file :
        parameters = AnalysisParameters.model_validate_json(parameters_file.read())

    return parameters

def exists_analysis_parameters(run_path) :
    return os.path.isfile(run_path + _ANALYSIS_PARAMETERS_FILENAME)

def _get_defined_variable(script) :
    variables = {var: val for var, val in vars(script).items() if not var.startswith("__")}
    
    return variables