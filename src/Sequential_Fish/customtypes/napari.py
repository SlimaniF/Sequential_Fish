"""
Super class for napari usage.
"""

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

class OrganoidWizard(ABC) :
    """
    Commong super class for wizards to launch only when viewer is open in the context of an organoid run.
    """

    @abstractmethod
    def start_listening(self) :
        """
        Principal action of widget is launched here.
        """