"""
Super class for napari usage.
"""

from magicgui.widgets import FunctionGui
from abc import ABC, abstractmethod

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
    def _create_widget(self) :
        """
        This should return a widget you can add to the napari (QWidget)
        """
        pass

    def register_widget(self, attr) :
        self.widgets.append(attr)

    def get_widgets(self) -> 'list[FunctionGui]' :
        return self.widgets