import pandas as pd
import numpy as np
import warnings
import os
import default_pipeline_parameters as parameters
from Sequential_Fish.tools import get_datetime
from Sequential_Fish import __run_cache_path__
from Sequential_Fish import __version__

COLUMNS = [
    "run_id",
    "creation_date",
    "last_modification_date",
    "pipeline_version"
]

from .._pipeline_scripts import PIPELINE_SCRIPTS

def _get_defined_variable(script) :
    variables = {var: val for var, val in vars(script).items() if not var.startswith("__")}
    for key, value in variables.copy().items() :
        if isinstance(value,dict) :
            variables.update(value)
            del variables[key]
    
    return variables

def create_run_dataframe() :
    """
    Create empty dataframe with correct columns and dtypes.
    """
    run_dataframe = pd.DataFrame(columns=COLUMNS, dtype='string')
    run_dataframe['run_id'] = run_dataframe['run_id'].astype(int)
    run_parameters = _get_defined_variable(parameters)

    for para in run_parameters.keys() :
        para_type = type(run_parameters[para]).__name__
        
        if para_type in ['tuple', 'list'] : 
            run_dataframe[para.upper()] = pd.Series([],dtype=object)
        elif para_type == 'dict' :
            raise AssertionError("Shouldn't use dict of dict in RUN PARAMETERS.")
        elif para_type == 'float' :
            run_dataframe[para.upper()] = pd.Series([],dtype= float)
        elif para_type == 'int' :
            run_dataframe[para.upper()] = pd.Series([],dtype=int)
        elif para_type == 'bool' :
            run_dataframe[para.upper()] = pd.Series([],dtype=bool)
        elif para_type == 'str' :
            run_dataframe[para.upper()] = pd.Series([],dtype='string')
        else :
            run_dataframe[para.upper()] = pd.Series([],dtype=object)

    for script in PIPELINE_SCRIPTS :
        run_dataframe[script] = pd.Series([],dtype=bool)

    return run_dataframe

def check_run_dataframe(run_dataframe: pd.DataFrame) :
    """
    Check for new parameters to be saved in Run_dataframe such as new scripts, or new user parameters.
    Also drop duplicates and na on RUN_PATH.
    """
    run_dataframe = run_dataframe.drop_duplicates(subset='RUN_PATH')
    run_dataframe = run_dataframe.dropna(subset='RUN_PATH')
    
        #Check run_id
    if 'index' in run_dataframe.columns : 
        run_dataframe = run_dataframe = run_dataframe.drop('index',axis=1)
    if len(run_dataframe) != len(run_dataframe['run_id'].unique()) or run_dataframe['run_id'].isna().any(): 
        warnings.warn("run_id was not unique in cache or contained na, cleaning run_ids...")
        run_dataframe = run_dataframe.drop('run_id', axis=1)
        run_dataframe = run_dataframe.reset_index(drop=True).reset_index(drop=False, names='run_id')
    
    model = create_run_dataframe()
    for col in model.columns :
        if not col in run_dataframe.columns :
            if model[col].dtype == 'bool' :
                run_dataframe[col] = False
            elif model[col].dtype == 'object' :
                run_dataframe[col] = pd.Series([tuple()] * len(run_dataframe))
            else :
                run_dataframe[col] = np.NaN
    
    return run_dataframe

def create_new_run(index=None) -> pd.DataFrame:
    """
    Create a new run.  
    Initialize a file into cache using parameters from **Sequential_Fish.pipeline_parameters.py**, then write parameters into **parameters.csv** in run_path directory.  
    **parameters.csv** is then used for pipeline parameters.
    """
    
    version = __version__
    date = get_datetime()
    run_parameters = _get_defined_variable(parameters)
    run_dataframe = get_run_cache()

    RUN_PATH = run_parameters['RUN_PATH']
    if RUN_PATH.endswith('/') : RUN_PATH = RUN_PATH[:-1]
    run_name = os.path.basename(RUN_PATH)

    if index is None :
        if len(run_dataframe) == 0 :
            index = pd.Index([0])
        else :
            index = pd.Index([run_dataframe.index.max() + 1])
    else :
        assert len(index) == 1, f"index of len != 1 given to create new run : {index}"

    if run_dataframe['run_id'].empty :
        new_run_id = 0
    elif index[0] in run_dataframe.index :
        new_run_id = run_dataframe.at[index[0], 'run_id']
    else :
        new_run_id = run_dataframe['run_id'].max() + 1
    
    new_run = pd.DataFrame({
        'run_id' : [new_run_id],
        'run_name' : [run_name],
        'creation_date' : [date],
        'last_modification_date' : [date],
        'pipeline_version' : [version]
    }, dtype='string')

    new_run['run_id'] = new_run['run_id'].astype(int)

        
    for para in run_parameters.keys() :
        para_type = type(run_parameters[para]).__name__
        
        if para_type in ['tuple', 'list'] : 
            new_run[para.upper()] = pd.Series([run_parameters[para]],dtype=object)
        elif para_type == 'dict' :
            raise AssertionError("Shouldn't use dict of dict in RUN PARAMETERS.")
        elif para_type == 'float' :
            new_run[para.upper()] = pd.Series([run_parameters[para]],dtype= float)
        elif para_type == 'int' :
            new_run[para.upper()] = pd.Series([run_parameters[para]],dtype=int)
        elif para_type == 'bool' :
            new_run[para.upper()] = pd.Series([run_parameters[para]],dtype=bool)
        elif para_type == 'str' :
            new_run[para.upper()] = pd.Series([run_parameters[para]],dtype='string')
        else :
            new_run[para.upper()] = pd.Series([run_parameters[para]],dtype=object)

    for script in PIPELINE_SCRIPTS :
        new_run[script] = pd.Series([False]).astype(bool)
    
    new_run.index = index
    
    return new_run

def add_new_run(run_dataframe : pd.DataFrame) :
    """
    Create a new line for Run dataframe. Check before in 'RUN_PATH' is already in df because it will erase it.
    """
    new_run = create_new_run()

    run_dataframe = pd.concat([
        run_dataframe,
        new_run
    ], axis= 0)
    
    save_run_cache(run_dataframe)

def get_run_cache() :
    
    if not os.path.isfile(__run_cache_path__) :
        print(f"Creating run_cache at {__run_cache_path__}")
        run_dataframe = create_run_dataframe()
        save_run_cache(run_dataframe)
    else :
        try :
            run_dataframe = pd.read_feather(__run_cache_path__)
        except Exception as e :
            warnings.warn("Could not open cache, reseting cache...")
            run_dataframe = create_run_dataframe()
    
    return run_dataframe.sort_values(['run_id'])

def save_run_cache(run_dataframe : pd.DataFrame) :
    """
    Safely save cache by converting objects to tuples.
    """
    
    obj_col = list(run_dataframe.select_dtypes(object,exclude="string").columns)
    for col in obj_col :
        try :
            run_dataframe[col] = run_dataframe[col].apply(tuple)
        except Exception :
            raise Exception(f"Could not convert col {col} to tuple.\n{run_dataframe[col]}")
    
    run_dataframe.reset_index(drop=True).to_feather(__run_cache_path__)
    return True

def _get_run_path_index(run_dataframe : pd.DataFrame, run_path : str) :
    
    assert len(run_dataframe['RUN_PATH'].unique()) == len(run_dataframe), "RUN_PATH is not unique in run_dataframe"
    
    index = run_dataframe.loc[run_dataframe['RUN_PATH'] == run_path].index
    
    if len(index) == 1 :
        pass
    else :
        index = None
    
    return index
    
def get_parameter(
    run_dataframe : pd.DataFrame,
    run_path : str,
    parameter : str
) :
    index = _get_run_path_index(run_dataframe, run_path)
    if len(index) == 1 : index = index[0]
    else : raise AssertionError()
    
    parameter = run_dataframe.at[index, parameter.upper()]
    
    if isinstance(parameter, (list, np.ndarray)) :
        parameter = list(parameter)
    
    return parameter

def get_parameter_dict(
    run_path : str,
    parameters : 'list[str]'
) :
    run_dataframe = get_run_cache()
    
    return {parameter : get_parameter(run_dataframe, run_path, parameter) for parameter in parameters}