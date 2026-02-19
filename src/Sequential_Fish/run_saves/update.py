import os, warnings
import pandas as pd
import numpy as np
from Sequential_Fish import __run_cache_path__
from ..tools import get_datetime
from .dataframe import add_new_run, _get_run_path_index, create_new_run, get_run_cache, save_run_cache, create_run_dataframe
from .._pipeline_scripts import PIPELINE_SCRIPTS


def validate_script(RUN_PATH, script: str) :
    """
    Opens and write in Run_cache that passed script was successfully runned and update datetime.
    """
    
    if script.endswith('.py') : script = script[:-3]
    script = os.path.basename(script)
        
    date = get_datetime()
    run_dataframe = get_run_cache()

    assert script in run_dataframe.columns, f"Script {script} was not found in run_dataframe table. Update run_dataframe or check script name matches _pipeline_scripts.py"
    
    target_index = run_dataframe[run_dataframe['RUN_PATH'] == RUN_PATH].index
    if len(target_index) == 1 :
        pass
    else :
        warnings.warn(f'Coul not validate script {script}, run path was not unique in cache for {RUN_PATH}\n {target_index}')
        return None
    
    run_dataframe.loc[target_index, [script]] = pd.Series([True], index=target_index, dtype=bool)
    run_dataframe.loc[target_index, ['last_modification_date']] = pd.Series([date], index=target_index, dtype='string')
    save_run_cache(run_dataframe)

def fail_script(RUN_PATH, script) :
    """
    Opens and write in Run_cache that passed script was successfully runned.
    
    """
    if script.endswith('.py') : script = script[:-3]
    script = os.path.basename(script)
    run_dataframe = get_run_cache()
    assert script in run_dataframe.columns, f"Script {script} was not found in run_dataframe table. Update run_dataframe or check script name matches _pipeline_scripts.py"
    
    target_index = run_dataframe[run_dataframe['RUN_PATH'] == RUN_PATH].index
    if len(target_index) == 1 :
        pass
    else :
        warnings.warn(f'Coul not validate script {script}, run path was not unique in cache for {RUN_PATH}')
        return None
    
    run_dataframe.loc[target_index, [script]] = pd.Series([False], index=target_index, dtype=bool)
    save_run_cache(run_dataframe)

def update_run(run_path) :
    """
    Update in cache parameters from pipeline_parameters.py for run
    """
    run_dataframe = get_run_cache()
    run_index =_get_run_path_index(run_dataframe, run_path)
    
    if run_index is None :
        run_index = pd.Index([0])
    else :
        run_dataframe = run_dataframe.drop(run_index,axis=0)
    new_run = create_new_run(index=run_index)
    print(new_run)
    
    run_dataframe = pd.concat([
        run_dataframe,
        new_run,
    ],axis = 0)
    
    save_run_cache(run_dataframe)
    
def add_path_to_cache(run_path) :
    """
    add run path to cache without specifying any parameters (has to go through update during pipeline)
    """
    
    run_dataframe = get_run_cache()
    
    if run_path in run_dataframe['RUN_PATH'] :
        return False
    
    if len(run_dataframe) == 0 :
        new_run_id = 0
        new_index = pd.Index([0])
    else :
        new_run_id = run_dataframe['run_id'].max() + 1
        new_index = pd.Index([run_dataframe['run_id'].index.max() + 1])
    
    model = create_run_dataframe()
    new_run = create_new_run(new_index)
    for col in model.columns : # Not using concat to avoid messing with dtypes
        if model[col].dtype == 'bool' :
            new_run[col] = pd.Series([False], dtype=bool, index=new_index)
        elif model[col].dtype == 'string' :
            new_run[col] = pd.Series([''],dtype='string', index=new_index)
        elif model[col].dtype == 'object' :
            new_run[col] = pd.Series([list()] * len(new_run), index=new_index)
        else :
            new_run[col] = pd.Series([np.nan], dtype=float, index=new_index)
    new_run.index = (new_index)
    
    new_run['run_id'] = pd.Series([new_run_id], dtype=int, index=new_index)
    new_run['RUN_PATH'] = pd.Series([run_path], dtype='string', index=new_index)
    new_run['run_name'] = pd.Series([os.path.basename(run_path)], dtype='string', index=new_index)
    
    
    run_dataframe = pd.concat([
        run_dataframe,
        new_run,
    ], axis=0)
    
    save_run_cache(run_dataframe)
    
    return True

def check_run(run_path) :
    """
    Check if folder is already in cached runs. If not initiate a new run with current parameters and scripts to False.
    """
    run_dataframe = get_run_cache()
    
    if run_path in list(run_dataframe['RUN_PATH']) : 
        print(f"Updating RUN : {run_path}")
        update_run(run_path)
    else : #
        print(f"Initiating new RUN : {run_path}")
        add_new_run(run_dataframe)
        
def run_status() :
    run_dataframe = get_run_cache()
    run_dataframe['folder'] = run_dataframe['RUN_PATH'].apply(os.path.basename)
    print(run_dataframe.loc[:,['run_id', 'folder','RUN_PATH',] + PIPELINE_SCRIPTS])
    