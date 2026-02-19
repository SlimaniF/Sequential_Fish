import sys
import os
import pandas as pd
import warnings

from Sequential_Fish import __run_cache_path__ as run_cache_path
from Sequential_Fish import viewer, pipeline, analysis
from Sequential_Fish.run_saves import create_run_dataframe, check_run_dataframe, check_run, run_status, get_run_cache
from Sequential_Fish._pipeline_scripts import PIPELINE_SCRIPTS
from Sequential_Fish.pipeline_parameters import RUN_PATH



def main():

    MODULES = ['viewer', 'pipeline', 'analysis', 'status']

    run_dataframe = get_run_cache()
    run_dataframe = check_run_dataframe(run_dataframe)
    run_dataframe.reset_index(drop=True).to_feather(run_cache_path)
    
    #CALL TO MODULES
    if len(sys.argv) < 2:
        print("Usage: python -m my_package <module> [args...]")
        print("Available modules: {0}".format(MODULES))
        sys.exit(1)

    module = sys.argv[1]
    submodules = sys.argv[2:]

    if module == "viewer":

        if len(submodules) > 0 :
            print(f"No argument for viewer. Ignoring passed arguments : {submodules}")

        viewer.run()
        
    elif module == "pipeline":
        
        check_run(RUN_PATH)
        
        if len(submodules) == 0 :
            pipeline.run() # This loads RUN_PATH from pipeline parameters and fix it for all scripts
        else :
            from Sequential_Fish.pipeline.runner import launch_script, script_folder # This loads RUN_PATH from pipeline parameters and fix it for all scripts
            if not all([script in PIPELINE_SCRIPTS for script in submodules]) :
                print(f"Unknown pipeline scripts. \nChoose from : {PIPELINE_SCRIPTS}")
            else :
                error_count = 0
                for script in submodules : 
                    if not script.endswith('.py') : script += ".py"
                    sucess = launch_script(script_folder + '/' + script, run_path=RUN_PATH)

                    if not sucess : error_count +=1
                
                print(f"all scripts ended with {error_count} errors. If any, error can be checked in run_log file.")
    
    elif module == "analysis" :
        
        if len(submodules) == 0 : 
            submodules = ['all']
            print("Starting all analysis modules")
        else :
            print("Starting selected analysis modules")
            
        analysis.run(*submodules)
        
        print("Done.")
    
    elif module == "status" :
        run_status()
    
    else:
        print(f"Unknown module: {module}")
        print("Available modules: {0}".format(MODULES))
        sys.exit(1)

if __name__ == "__main__":
    main()
