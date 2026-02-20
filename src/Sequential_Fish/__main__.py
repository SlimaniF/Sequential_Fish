import sys
import os
import logging

from Sequential_Fish import viewer, pipeline, analysis
from ._pipeline_scripts import PIPELINE_SCRIPTS
from .settings import write_settings
from .pipeline_parameters import get_default_settings



def main():

    MODULES = ['viewer', 'pipeline', 'analysis']

    #CALL TO MODULES
    if len(sys.argv) < 3:
        print("Usage: python -m my_package <module> [args...] <run path>")
        print("Available modules: {0}".format(MODULES))
        sys.exit(1)

    module = sys.argv[1]
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

        if not os.path.isfile(RUN_PATH + "/settings.json") :
            logging.info("No settings found, initializing settings.json")
            settings = get_default_settings()
            settings.RUN_PATH = RUN_PATH
            settings.SAVE_PATH = RUN_PATH + "/visuals/"
            write_settings(settings, RUN_PATH)
        else :
            logging.info("Loading parameters")
        
        if len(submodules) == 0 :
            pipeline.run(RUN_PATH) # This loads RUN_PATH from pipeline parameters and fix it for all scripts
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
    
    else:
        print(f"Unknown module: {module}")
        print("Available modules: {0}".format(MODULES))
        sys.exit(1)

if __name__ == "__main__":
    main()
