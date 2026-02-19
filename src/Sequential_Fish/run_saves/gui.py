"""
Submodule handling data opening and merging for viewer
"""
import os
from PyQt5.QtWidgets import QApplication, QDialog, QListWidget, QVBoxLayout, QPushButton, QHBoxLayout, QFileDialog
from .update import get_run_cache, add_path_to_cache

TABLES = ['Acquisition', 'Detection', 'Spots', 'Clusters', 'Drift', 'Cell', 'Gene_map']

class PathSelector(QDialog):
    """
    Window for path selection
    """
    def __init__(self, path_map):
        super().__init__()
        self.setWindowTitle("Select Paths")
        self.setGeometry(100, 100, 400, 300)

        self.path_map = path_map  # Dictionary {custom_name: real_path}

        # List Widget
        self.list_widget = QListWidget(self)
        self.list_widget.addItems(path_map.keys())  # Display custom names

        # OK Button
        self.ok_button = QPushButton("OK", self)
        self.ok_button.clicked.connect(self.accept)
        # Load button
        self.load_button = QPushButton("Load", self)
        self.load_button.clicked.connect(self.load_folder)
        
        #Buttons horizontal layout
        button_layout = QHBoxLayout()
        button_layout.addWidget(self.load_button)
        button_layout.addWidget(self.ok_button)
        
        # Layout
        layout = QVBoxLayout()
        layout.addWidget(self.list_widget)
        layout.addLayout(button_layout)
        self.setLayout(layout)

    def load_folder(self) :
        folder = QFileDialog.getExistingDirectory(self, "Select Folder")
        if folder:
            if self.check_folder(folder) :
                name = os.path.basename(folder)
                self.path_map[name] = folder
                self.list_widget.addItem(name)
            else :
                print("Could not load folder.")

    def check_folder(self,folder) :
        dirlist = os.listdir(folder)
        
        if 'result_tables' in dirlist :
            data_tables_dirlist = os.listdir(folder + '/result_tables/')
            TABLES = ['Acquisition', 'Spots', 'Clusters', 'Drift', 'Detection', 'Gene_map']
            tables_test = all([(table + '.feather') in data_tables_dirlist] for table in TABLES)
            if not tables_test :
                print(f"Couln'd find all required dataframes in .feather.\nRequired list : {TABLES}")
        else :
            print("Couldn't find result_tables folder in selected folder. Please select main folder of a Sequential_Fish acquisition.")
            tables_test = False
            
        if 'segmentation' in dirlist :
            segmentation_test = True
        else :
            print("Couldn't find segmentation folder in selected folder. Please select main folder of a Sequential_Fish acquisition.")
            segmentation_test = False
        
        res = tables_test and segmentation_test
        if res : add_path_to_cache(folder)
        
        return res

    def get_selected_paths(self):
        selected_names = [item.text() for item in self.list_widget.selectedItems()]
        
        res = [self.path_map[name] for name in selected_names]
        if len(res) == 1 :
            return res[0]
        elif len(res) > 1 :
            raise AssertionError("size should be 1 or 0 not ", len(res), res)
        else : 
            return None

def select_path():
    """
    Open graphical interface to select a run. Will suggest runs that are saved in cache.
    """
    
    run_cached = get_run_cache()
    
    path_list = list(run_cached['RUN_PATH'].unique())
    
    selection_path = {}
    for path in path_list :
        if path.endswith('/') : path = path[:-1]
        selection_path[os.path.basename(path)] = path
    
    app = QApplication([])
    dialog = PathSelector(selection_path)
    if dialog.exec():  # Show dialog and check if OK was pressed
        return dialog.get_selected_paths()
    return None