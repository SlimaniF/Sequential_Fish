"""
Entry point to launch Sequential Fish analysis pipeline.
"""

__version__ = 0.1

def run(*args) :
    from .runner import main
    main()