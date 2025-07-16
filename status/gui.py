"""
Submodule handling data opening and merging for viewer
"""
import os
import pandas as pd

from PyQt5.QtWidgets import QApplication, QDialog, QListWidget, QVBoxLayout, QPushButton, QHBoxLayout, QFileDialog, QListWidgetItem, QLabel
from PyQt5.QtGui import QColor
from PyQt5.QtCore import Qt

from .pipeline_parameters import get_raw_pipeline_parameters, write_pipeline_parameters
from .cache import read_cache

TABLES = ['Acquisition', 'Detection', 'Spots', 'Clusters', 'Drift', 'Cell', 'Gene_map']

TOOLTIPS = {

    "check"     :   
            """Check that runs in cache are found on hardrive.
            red : Folder not compatible with SequentialFish software.
            orange : Parameters and images found. Can be used in pipeline
            green : Segmentation and detection results found. Can be used in analysis.""",

    "rm"        :   """Delete run from cache (not from hardrive).""",
    "add"       :   """Load run from hardrive into cached runs.""",
}

class CacheDialog(QDialog) :
    """
    PyQt5 Dialog basis for cache management
    """
    def __init__(self, path_map : dict):
        super().__init__()
        self.setWindowTitle("Select Paths")
        self.setGeometry(100, 100, 400, 300)

        self.path_map = path_map  # Dictionary {custom_name: real_path}

        # List Widget
        self.list_widget = QListWidget(self)
        self.list_widget.addItems(path_map.keys())  # Display custom names
        
        # Layout
        layout = QVBoxLayout()
        layout.addWidget(self.list_widget)
        self.setLayout(layout)

    def get_selected_path(self):
        selected_names = [item.text() for item in self.list_widget.selectedItems()]
        
        res = [self.path_map[name] for name in selected_names]
        if len(res) == 1 :
            return res[0]
        elif len(res) > 1 :
            raise AssertionError("size should be 1 or 0 not ", len(res), res)
        else : 
            return None

class PathSelector(CacheDialog):
    """
    Window for path selection
    """
    def __init__(self, path_map : dict):
        super().__init__(path_map)

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
        self.layout().addLayout(button_layout)

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
        # if res : add_path_to_cache(folder)
        
        return res

class CacheUpdater(CacheDialog):
    """
    Class for update run in cache (adding, removing or checking that on drive and in Acquisition match.
    """

    def __init__(self, path_map : dict):
        super().__init__(path_map)
        self.check_cache()

        header = QLabel("Sequential Fish runs", self)
        header.setAlignment(Qt.AlignCenter)
        header.setStyleSheet("""
            font-size: 16pt;
            font-weight: bold;
            margin-bottom: 10px;
        """)
        self.layout().insertWidget(0,header)

        # Buttons init
        self.check_button = QPushButton("Modify parameters", self)

        self.add_button = QPushButton("Add",self)
        self.add_button.setToolTip(TOOLTIPS['add'])
        self.add_button.clicked.connect(self.add_to_cache)
        
        self.rm_button = QPushButton("Remove",self)
        self.rm_button.setToolTip(TOOLTIPS['rm'])
        self.rm_button.clicked.connect(self.rm_from_cache)

        buttons = [self.check_button, self.add_button, self.rm_button]

            ## Add to layout
        button_Hbox = QHBoxLayout()
        for button in buttons :
            button_Hbox.addWidget(button)
        self.layout().addLayout(button_Hbox)

        #Ok / Cancel
        self.ok_button = QPushButton("Ok", self)
        self.ok_button.clicked.connect(self.accept)
        self.cancel_button = QPushButton("Cancel",self)
        self.cancel_button.clicked.connect(self.reject)

            ## Add to layout
        OkCancel_Hbox = QHBoxLayout()
        OkCancel_Hbox.addWidget(self.ok_button)
        OkCancel_Hbox.addWidget(self.cancel_button)
        self.layout().addLayout(OkCancel_Hbox)

    
    def check_cache(self) :
        run_number = len(self.list_widget)
        for run in range(run_number) :
            
            item = self.list_widget.item(run)
            run_path = self.path_map[item.text()]

            if not os.path.isdir(run_path) :
                print("Couldn't find folder at {}".format(run_path))
                item.setForeground(QColor("red"))
            
            elif not check_Acquisition_run_path(run_path) :
                item.setForeground(QColor("red"))

            else :
                item.setForeground(QColor("green"))


    def add_to_cache(self) :
        folder = QFileDialog.getExistingDirectory(self, "Select main folder of run to add")
        if folder :
            color = "red"
            if check_Acquisition_run_path(folder) :
                color = "orange"
            if check_analysis_completed(folder) :
                color = "green"
            name = os.path.basename(folder)
            if name in self.path_map.keys() :
                print("Already a run with folder name {}, delete or rename folder first.".format(name))
                return None
            else : 
                self.path_map[name] = folder

                if "pipeline_parameters.json" in os.listdir(folder + "/") :
                    pass
                else :
                    pipeline_parameters = get_raw_pipeline_parameters()
                    write_pipeline_parameters(pipeline_parameters)

                new_run_item = QListWidgetItem(name)
                new_run_item.setForeground(QColor(color))
                self.list_widget.addItem(new_run_item)

    def rm_from_cache(self) :
        for item in self.list_widget.selectedItems() :
            self.list_widget.takeItem(self.list_widget.row(item))
            del self.path_map[item.text()]



def check_Acquisition_run_path(folder) :
    """
    Checks that Acquisition.feather is found in result_tables and that path written inside is found on hardrive
    """
    files = os.listdir(folder)
    if not "result_tables" in files :
        print("Couldn't find result tables folder at {}".format(folder))
        return False
    
    result_tables = os.listdir(folder + "/result_tables/")
    if "Acquisition.feather" not in result_tables :
        print("Couldn't find Acquisition table at {}".format(folder))
        print(result_tables)
        return False
    
    Acquisition = pd.read_feather("{0}/result_tables/Acquisition.feather".format(folder), columns=["acquisition_id", "full_path"])
    for path_to_check in Acquisition['full_path'].unique().tolist() :
        if path_to_check is None : continue
        if not os.path.isfile(path_to_check) :
            print("Path written in cached Acquisition not found on hardrive : {}".format(path_to_check))
    return True



def check_analysis_completed(folder) :
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
        # if res : add_path_to_cache(folder)
        
        return res

def select_path_for_pipeline():
    """
    Open graphical interface and cached runs. User can add a new folder and modify pipeline parameters.  
    Adding a new folder will create a default parameters configuration if not modified by user.
    """
    selection_path = read_cache()

    app = QApplication([])
    dialog = CacheUpdater(selection_path)
    if dialog.exec():  # Show dialog and check if OK was pressed
        return dialog.get_selected_path()
    return None


def select_path_for_analysis():
    """
    Open graphical interface to select a run. Will suggest runs that are saved in cache.
    """
 
    selection_path = read_cache()
    
    app = QApplication([])
    dialog = PathSelector(selection_path)
    if dialog.exec():  # Show dialog and check if OK was pressed
        return dialog.get_selected_path()
    return None