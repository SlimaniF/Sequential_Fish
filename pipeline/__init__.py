"""
Entry point to launch Sequential Fish analysis pipeline.
"""

from .input import main as input
from .detection import main as detection
from .segmentation import main as segmentation
from .drift import main as drift
from .alignement import main as alignement
from .washout import main as washout
from .quantification import main as quantification

PIPELINE_SCRIPTS = {
    "input" : input,
    "detection" : detection,
    "segmentation" : segmentation,
    "drift" : drift,
    "alignement" : alignement,
    "washout" : washout,
    "quantification" : quantification,
}

def run(run_path : str, *args) :
    from .runner import main
    main(run_path)