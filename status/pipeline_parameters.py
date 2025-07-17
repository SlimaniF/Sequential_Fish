from ..customtypes import PipelineParameters
import Sequential_Fish.default_pipeline_parameters as parameters

def get_raw_pipeline_parameters() :
    run_parameters = _get_defined_variable(parameters)
    return PipelineParameters(**run_parameters)


def write_pipeline_parameters(run_path, parameters : PipelineParameters) :
    if not isinstance(parameters, PipelineParameters) : raise TypeError("Invalid type for parameters argument : {}; custom PipelineParameters type is mandatory".format(type(parameters)))

    with open(run_path + "/pipeline_parameters.json", "w") as f:
        f.write(parameters.model_dump_json(indent=2))

def load_pipeline_parameters(parameters_path) : 
    parameters = PipelineParameters.model_validate_json(parameters_path)

    return parameters

def _get_defined_variable(script) :
    variables = {var: val for var, val in vars(script).items() if not var.startswith("__")}
    
    return variables