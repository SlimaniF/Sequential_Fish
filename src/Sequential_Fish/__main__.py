import sys
<<<<<<< HEAD
import os
import logging

from . import viewer, pipeline, analysis
from ._pipeline_scripts import PIPELINE_SCRIPTS
from .settings import write_settings
from .types import pipeline_parameters
from .pipeline.runner import launch_script 

=======

from Sequential_Fish import viewer, pipeline, analysis, chromatic_abberrations, status
from Sequential_Fish.status import select_path_for_pipeline, select_path_for_analysis

from .pipeline import PIPELINE_SCRIPTS
>>>>>>> Modification_drift_correction


def main():

<<<<<<< HEAD
    MODULES = ['viewer', 'pipeline', 'analysis']

    #CALL TO MODULES
    if len(sys.argv) < 3:
        print("Usage: python -m my_package <module> [args...] <run path>")
=======
    MODULES = ['viewer', 'pipeline', 'analysis', 'status', 'calibration']

    
    #CALL TO MODULES
    if len(sys.argv) < 2:
        print("Usage: python -m my_package <module> [args...]")
>>>>>>> Modification_drift_correction
        print("Available modules: {0}".format(MODULES))
        sys.exit(1)

    module = sys.argv[1]
<<<<<<< HEAD
    print("module : ", module)
    submodules = sys.argv[2:-1]
    print("submodules : ", submodules)
    RUN_PATH = sys.argv[-1]
    print("RUN_PATH : ", RUN_PATH)
    if not os.path.isdir(os.path.join(os.getcwd(), RUN_PATH)) :
        raise FileNotFoundError(f"Couldn't find a directory at {RUN_PATH}")
    else : RUN_PATH = str(os.path.join(os.getcwd(), RUN_PATH))
     
    if module == "viewer":
        if len(submodules) > 0 :
            print(f"No argument for viewer. Ignoring passed arguments : {submodules}")
        viewer.run()
        
    elif module == "pipeline":

        log_file = RUN_PATH + "/run_log.log"
        logging.basicConfig(
            filename=log_file,
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
        )

        if not os.path.isfile(RUN_PATH + "/pipeline_settings.json") :
            logging.info("No settings found, initializing settings.json")
            settings = pipeline_parameters.from_default_parameters()
            settings.RUN_PATH = RUN_PATH
            settings.SAVE_PATH = RUN_PATH + "/visuals/"
            write_settings(settings, RUN_PATH)
        else :
            logging.info("Loading parameters")
        
        if len(submodules) == 0 :
            pipeline.run(RUN_PATH)
        else :
            if not all([script in PIPELINE_SCRIPTS for script in submodules]) :
                print(f"Unknown pipeline scripts. \nChoose from : {PIPELINE_SCRIPTS}")
            else :
                error_count = 0
                for script in submodules : 
                    sucess = launch_script(script, run_path=RUN_PATH)
=======
    submodules = sys.argv[2:]

    if module == "viewer":

        if len(submodules) > 0 :
            print(f"No argument for viewer. Ignoring passed arguments : {submodules}")

        viewer.run()
        
    elif module == "pipeline":
        
        if '-p' in submodules : # for run in command line only : run path is passed after -p flag
            run_path_index = submodules.index('-p')
            run_path = submodules[run_path_index+1]
            submodules.pop(run_path_index)
            submodules.pop(run_path_index)
        else :
            run_path = select_path_for_pipeline() 
            if run_path is None : quit()

        if len(submodules) == 0 :
            pipeline.run(run_path) # This loads RUN_PATH from pipeline parameters and fix it for all scripts
        else :
            from Sequential_Fish.pipeline.runner import launch_script, script_folder # This loads RUN_PATH from pipeline parameters and fix it for all scripts
            if not all([script in PIPELINE_SCRIPTS.keys() for script in submodules]) :
                print(f"Unknown pipeline scripts. \nChoose from : {PIPELINE_SCRIPTS.keys()}")
            else :
                error_count = 0
                for script in submodules : 
                    sucess = launch_script(script_name=script, script=PIPELINE_SCRIPTS[script], run_path=run_path)
>>>>>>> Modification_drift_correction

                    if not sucess : error_count +=1
                
                print(f"all scripts ended with {error_count} errors. If any, error can be checked in run_log file.")
    
    elif module == "analysis" :
<<<<<<< HEAD
=======

        if '-p' in submodules : # for run in command line only : run path is passed after -p flag
            run_path_index = submodules.index('-p')
            run_path = submodules[run_path_index+1]
            submodules.pop(run_path_index)
            submodules.pop(run_path_index)
        else :
            run_path = select_path_for_analysis() 
>>>>>>> Modification_drift_correction
        
        if len(submodules) == 0 : 
            submodules = ['all']
            print("Starting all analysis modules")
        else :
            print("Starting selected analysis modules")
            
<<<<<<< HEAD
        analysis.run(*submodules)
        
        print("Done.")
    
    else:
        print(f"Unknown module: {module}")
        print("Available modules: {0}".format(MODULES))
=======
        analysis.run(run_path, *submodules)
        
        print("Done.")
    
    elif module == "status" :
        status.run()

    elif module == "calibration" :
        chromatic_abberrations.run()
    
    else:
        print(f"Unknown module: {module}")
        print("Available modules: {0}".format(MODULES))
        print("To skip path selection, use '-p' flag followed with fullpath to your directory.")
>>>>>>> Modification_drift_correction
        sys.exit(1)

if __name__ == "__main__":
    main()
