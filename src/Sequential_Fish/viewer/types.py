"""
Super class for napari usage.
"""
from .utils import update_layer_from_LayerDataTuple

from napari import Viewer
from napari.qt.threading import create_worker
from napari.components import ViewerModel
from magicgui.widgets import FunctionGui
from abc import ABC, abstractmethod
from typing import Any

class NapariWidget(ABC) :
    """
    Common super class for custom widgets added to napari interface during run
    Each sub class as a specific function, but the widget can be acess with attribute .widget
    """
    def __init__(self):
        self.widget = self._create_widget()
        self.widgets = []
        self.register_widget(self.widget)

    @abstractmethod
    def _create_widget(self) -> FunctionGui:
        """
        This should return a widget you can add to the napari (QWidget)
        """

    def __call__(self, *args: Any, **kwds: Any) -> Any:
        return None

    def register_widget(self, attr) :
        self.widgets.append(attr)

    def get_widgets(self) -> 'list[FunctionGui]' :
        return self.widgets
    
    def disable_widget(self) :
        for w in self.widgets :
            w.enabled = False

class NapariWizard(ABC) :
    """
    Commong super class for wizards to launch only when viewer is open in the context of an organoid run.
    """
    def __init__(self, napari_window : ViewerModel) -> None:
        pass

    @abstractmethod
    def start_listening(
        self, 
        ) :
        """
        Principal action of widget is launched here.
        """

class ThreadedWidget(NapariWidget) :
    """
    NapariWidgets inheriting from this class will launch the widget function in a concurent thread.
    Use only if widget returnsa LayerDataTuple
    """
    
    def __init__(self, *, viewer : Viewer):
        super().__init__()
        self.viewer = viewer
        self._wrap_with_threading()

    def _wrap_with_threading(self):
        """
        Wraps the widget's function with napari's threading.
        """
        if hasattr(self.widget, 'native'):
            # Get the underlying function from the magicgui widget
            original_func = self.widget._function
            def threaded_func(*args, **kwargs):
                worker = create_worker(original_func, *args, **kwargs)
                worker.returned.connect(lambda result: update_layer_from_LayerDataTuple(self.viewer, result))
                worker.start()
                return worker
            # Replace the widget's function with the threaded version
            self.widget._function = threaded_func


class UserInputError(Exception) :
    """
    Exception raised when a widget cannot function due to an error in its usage by the user.
    """