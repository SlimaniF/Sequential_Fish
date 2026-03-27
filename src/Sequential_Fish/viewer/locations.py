"""
Widgets handlig location selection in Napari viewer
"""

from typing import cast

import pandas as pd
import napari

from magicgui import magicgui
from tqdm import tqdm

from ..customtypes import NapariWidget
from ..customtypes import table_dict_type



##  Location tab
_LOCATION_WIDGETS = []
def register_location_widget(cls) :
    _LOCATION_WIDGETS.append(cls)
    return cls

def initiate_location_widgets(
        *,
        tables_dict,
        Viewer,
        linked_widgets,
) :
    widget_list = []
    for cls in _LOCATION_WIDGETS :
        widget_list.extend(cls(table_dict=tables_dict, Viewer=Viewer, linked_widgets=linked_widgets).get_widgets())
    return widget_list

#Location widget
@register_location_widget
class LocationSelector(NapariWidget) :
    def __init__(
            self, 
            table_dict: table_dict_type, 
            Viewer : napari.Viewer, 
            linked_widgets : list,
            **_
            ):
        self.Full_Acquisiton = table_dict['Acquisition'].copy()
        location_choices = list(self.Full_Acquisiton['location'].unique())
        self.location_choices = pd.Series(location_choices)
        self.selection = pd.Series(self.location_choices)
        self.linked_widgets = linked_widgets
        self.Viewer = Viewer
        super().__init__()
        
    def update_location(self) :
        for layer in self.Viewer.layers.copy() :
            self.Viewer.layers.remove(layer)
        self.Viewer.reset_view()

        for widget in tqdm(self.linked_widgets, desc= "Updating locations") : 
            widget.update(self.selection)
            if hasattr(widget,"widget") :
                widget.widget.update()

    def _create_widget(self) :
        @magicgui(
            selected_location={
                "widget_type" : "Select",
                "choices" : list(self.location_choices),
                "value" : list(self.location_choices),
                "label" : ' ',
            },
            call_button= "select locations"
        )
        def select_location(selected_location) :
            print("Selected locations : ", selected_location)
            self.selection = cast(pd.Series, self.location_choices[cast(slice, self.location_choices.isin(selected_location))])
            self.selection.sort_values()
            self.update_location()
        
        
        return select_location