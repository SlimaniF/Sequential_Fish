"""
Entry point to launch Sequential Fish analysis pipeline.
"""

def run(run_path : str, *args) :
    from .runner import main
    main(run_path)