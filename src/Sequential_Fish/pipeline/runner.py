import warnings
import logging
import os, traceback
from datetime import datetime

from ..status import get_raw_pipeline_parameters, write_pipeline_parameters

from .input import main as input
from .detection import main as detection
from .segmentation import main as segmentation
from .drift import main as drift
from .alignement import main as alignement
from .washout import main as washout
from .quantification import main as quantification

"""
Main script to launch pipeline

Script are listed in `scripts_rounds` in corresponding order.  
All scripts of a round are launched even if an error is returned within the round but next round will be launched only if all previous scripts were sucessful.

"""
script_folder = os.path.abspath(__file__)
script_folder = os.path.dirname(script_folder)


scripts_rounds = [
    {'input' : input},
    {'detection' : detection, 'drift' : drift},
    {'segmentation' : segmentation},
    {'alignement' : alignement},
    {'washout' : washout},
    {'quantification' : quantification}
    ]

def log_was_initiated() :
    return logging.getLogger().hasHandlers()

def initiate_log_config(run_path : str) :
    log_file = run_path + "/run_log.log"
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        )


### script launcher function
def launch_script(script, script_name, run_path):
    """Launch script from pipeline."""
    
    print("launching : ", script_name)

    if not log_was_initiated() : 
        print("run log initiated at : {}".format(run_path))
        initiate_log_config(run_path)

    if not os.path.isfile(run_path + "/pipeline_parameters.json") :
        print("No existing pipeline parameters configuration, creating one from default_pipeline_parameters.py.")
        pipeline_parameters = get_raw_pipeline_parameters()
        write_pipeline_parameters(run_path, pipeline_parameters)
    
    try:
        logging.info(f"Exécution du script : {script_name}")
        script_start = datetime.now()
        script(run_path)
        script_end = datetime.now()
        run_duration = (script_end - script_start).total_seconds()
        result = True
        

    except Exception as e:
        script_end = datetime.now()
        run_duration = (script_end - script_start).total_seconds()

        logging.info("script duration : {0}".format(run_duration))
        logging.error(f"script failed {script_name}:\n{traceback.format_exc()}")

        return False
    
    else :
        
        logging.info("script duration : {0}".format(run_duration))
        logging.info(f"script succeed {script_name}:\n")

        return True

def main(run_path : str):
    start_time = datetime.now()

    if not log_was_initiated() : 
        print("run log initiated at : {}".format(run_path))
        initiate_log_config(run_path)


    logging.info("NEW RUN")
    
    for round in scripts_rounds :
        sucess = []
        
        for script_name, script in round.items() :
            result = launch_script(script, script_name, run_path=run_path)
            sucess.append(result)

        if all(sucess) :
            logging.info("Step {0} succeed.".format(round.keys()))
        else :
            logging.error("Step {0} failed, ending run.".format(round.keys()))
            quit()

    end_time = datetime.now()
    duration = (end_time - start_time).total_seconds()
    logging.info(f"Run ends. Total duration : {duration:.2f} seconds")

if __name__ == "__main__":
    main()