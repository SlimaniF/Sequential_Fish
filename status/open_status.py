import os
from .cache import read_cache
from PyQt5.QtWidgets import QApplication
from .gui import CacheUpdater

def main():
    """
    Open graphical interface to select a run. Will suggest runs that are saved in cache.
    """
    
    run_cached = read_cache()
    
    path_list = list(run_cached['RUN_PATH'].unique())
    
    selection_path = {}
    for path in path_list :
        if path.endswith('/') : path = path[:-1]
        selection_path[os.path.basename(path)] = path
    
    app = QApplication([])
    dialog = CacheUpdater(selection_path)
    if dialog.exec():  # Show dialog and check if OK was pressed
        print(dialog.get_selected_paths())
        return dialog.get_selected_paths()
    else : 
        print("Canceled")
    return None