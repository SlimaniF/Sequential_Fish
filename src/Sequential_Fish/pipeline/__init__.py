"""
Entry point to launch Sequential Fish analysis pipeline.
"""

__version__ = 0.1

def run(RUN_PATH, *_) :
    from .runner import main
    main(RUN_PATH)