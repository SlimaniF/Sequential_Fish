from .cache import read_cache
from PyQt5.QtWidgets import QApplication
from .gui import PipelineCacheDialog

def main():
    """
    Open graphical interface to select a run. Will suggest runs that are saved in cache.
    """
    
    selection_path = read_cache()
    
    app = QApplication([])
    dialog = PipelineCacheDialog(selection_path)
    if dialog.exec():  # Show dialog and check if OK was pressed
        print(dialog.get_selected_path())
        return dialog.get_selected_path()
    else : 
        print("Canceled")
    return None