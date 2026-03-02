import sys
import os
import logging

from . import viewer, pipeline, analysis, chromatic_abberrations
from .customtypes.parameters import PipelineParameters
from .settings import write_settings
from .pipeline import launch_script

def main():

    MODULES = ['viewer', 'pipeline', 'analysis', 'settings', 'calibration']

    #CALL TO MODULES
    if len(sys.argv) < 2:
        print("Usage: python -m my_package <module> [args...]")
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
        viewer.run(RUN_PATH)
        
    elif module == "pipeline":

        log_file = RUN_PATH + "/run_log.log"
        logging.basicConfig(
            filename=log_file,
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
        )

        logging.info("\n\nSequential Fish is starting\n")


        if not os.path.isfile(RUN_PATH + "/pipeline_settings.json") :
            logging.info("No settings found, initializing settings.json")
            settings = PipelineParameters.from_default_parameters()
            settings.RUN_PATH = RUN_PATH
            settings.SAVE_PATH = RUN_PATH + "/visuals/"
            write_settings(settings, RUN_PATH)
        else :
            logging.info("Loading parameters")
        
        if len(submodules) == 0 :
            pipeline.run(RUN_PATH)
        else :
            error_count = 0
            for script in submodules : 
                sucess = launch_script(script, run_path=RUN_PATH)

                if not sucess : 
                    print(f"Error caught during computing. Aborting script {script}")
                    error_count +=1
            
            print(f"all scripts ended with {error_count} errors. If any, error can be checked in run_log file.")
    
    elif module == "analysis" :
        
        if len(submodules) == 0 : 
            submodules = ['all']
            print("Starting all analysis modules")
        else :
            print("Starting selected analysis modules")
            
        analysis.run(*submodules)
        
        print("Done.")

    elif module == "chromatic_abberations" :
        chromatic_abberrations.run()
    
    elif module == "settings" :
        raise NotImplementedError()
    
    else:
        print(f"Unknown module: {module}")
        print("Available modules: {0}".format(MODULES))
        sys.exit(1)

if __name__ == "__main__":
    main()
