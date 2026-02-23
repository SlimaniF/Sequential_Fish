"""
Main script to launch pipeline

Script are listed in `scripts_rounds` in corresponding order.  
All scripts of a round are launched even if an error is returned within the round but next round will be launched only if all previous scripts were sucessful.

"""
import logging
import os, traceback
from datetime import datetime

#pipeline scripts
from .input import main as run_input
from .alignement import main as run_alignement
from .detection import main as run_detection
from .drift import main as run_drift
from .quantification import main as run_quantification
from .segmentation import main as run_segmentation
from .washout import main as run_washout

script_folder = os.path.abspath(__file__)
script_folder = os.path.dirname(script_folder)

scripts_rounds = {
    'input' : ['input.py'],
    'analysis' : ['detection.py', 'segmentation.py', 'drift.py'],
    'alignement' : ['alignement.py'],
    'filtering' : ['washout.py'],
    'quantification' : ['quantification.py']
}

### script launcher function
def launch_script(script_name, run_path):
    """Launch script from pipeline."""
    logging.info(f"Exécution du script : {script_name}")
    script_start = datetime.now()
    try:
        match script_name :
            case "input" : run_input(run_path)
            case "detection" : run_detection(run_path)
            case "segmentation" : run_segmentation(run_path)
            case "drift" : run_drift(run_path)
            case "alignement" : run_alignement(run_path)
            case "washout" : run_washout(run_path)
            case "quantification" : run_quantification(run_path)
            case _ : raise ValueError(f"Unsupported script name : {script_name}")

        script_end = datetime.now()
        run_duration = (script_end - script_start).total_seconds()
        

    except Exception as e:
        script_end = datetime.now()
        run_duration = (script_end - script_start).total_seconds()

        logging.info(f"script duration : {run_duration}")
        logging.error(f"script failed {script_name}:\n{traceback.format_exc()}")
        logging.error(f"CalledProcessError {script_name}:\n{str(e)}")

        return False
    
    else :
        
        logging.info(f"script duration : {run_duration}")
        logging.info(f"script succeed {script_name}:\n")

        return True

def main(
    run_path : str
    ):

    start_time = datetime.now()
    logging.info("NEW RUN")
    
    for round_key in scripts_rounds.keys() :
        scripts = scripts_rounds[round_key]
        sucess = []
        for script in scripts:
            if os.path.exists(script_folder + '/' + script):
                result = launch_script(script_folder + '/' + script, run_path=run_path)
                sucess.append(result)
            else:
                logging.warning(f"File {script} not found.")
                sucess.append(False)

        if all(sucess) :
            logging.info(f"Step {round_key} succeed.")
        else :
            logging.error(f"Step {round_key} failed, ending run.")
            quit()

    end_time = datetime.now()
    duration = (end_time - start_time).total_seconds()
    logging.info(f"Run ends. Total duration : {duration:.2f} seconds")