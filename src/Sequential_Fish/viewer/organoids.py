import os
import warnings
from typing import cast

from tqdm import tqdm
import pandas as pd
from napari.components import ViewerModel
from napari.layers import Image, Points, Labels
from qtpy.QtCore import QTimer


from ..customtypes import NapariWizard
from ..customtypes.organoids import OrganoidLocations, load_organoid_locations, get_locations_export, LocationsDataStructureError


def autodetect_organoids(
    run_path: str,
    Acquisition : pd.DataFrame
    ) -> OrganoidLocations | None :

    location_list = list(Acquisition["location"].unique())

    organoids_locations = None
    for filename in os.listdir(run_path):
        if not filename.endswith(".json") :
            continue
        
        for location_export in get_locations_export().keys() :
            try :
                organoids_locations = load_organoid_locations(
                    json_data= os.path.join(run_path,filename),
                    location_export=location_export
                )
                if not organoids_locations.validate(location_list) :
                    warnings.warn("Organoid locations metadata was found but could not find position for all locations.")
                    return None

            except LocationsDataStructureError :
                continue
            else :
                break
        
        if not organoids_locations is None : break

    return organoids_locations

class OrganoidWizard(NapariWizard) :
    pass


_ORGANOID_WIZARDS : list = []
def register_organoid_wizard(cls) :
    _ORGANOID_WIZARDS.append(cls)
    return cls

def initiate_organoid_wizards(
    organoids_locations : OrganoidLocations
) -> list:
    wizard_list = []
    for cls in _ORGANOID_WIZARDS :
        instance = cls(organoids_locations)
        if hasattr(instance,"enabled") :
            if instance.enabled :
                wizard_list.extend(instance.get_wizards())
        else :
            wizard_list.extend(instance.get_wizards())

    return wizard_list


@register_organoid_wizard
class OrganoidMapper(OrganoidWizard) :

    def __init__(self, viewer : ViewerModel) -> None:
        self.viewer = viewer
        self.listener = None
        self.start_listening()

    def start_listening(self) :
        if self.listener is None :
            print("Start listening")
            self.listener = self.viewer.layers.events.inserted.connect(self.place_layers_on_map)

    def stop_listening(self) :
        if not self.listener is None :
            print("Stop listening")
            self.viewer.layers.events.inserted.disconnect(self.listener)
            self.listener = None
    
    def place_layers_on_map(self, event) :

        try :
            self.stop_listening()
            layer = event.value
            print(layer.name)
            print(event, "from ", layer,"\n Starting scattering")
            print(layer.data.shape)

            location_index = 0
            for layer_data in tqdm(layer.data, desc="Mapping field of views") :
                if isinstance(layer, Image) :
                    new_fov = self.viewer.add_image(layer_data)
                    new_fov = cast(Image, new_fov)

                elif isinstance(layer, Points) :
                    new_fov = self.viewer.add_points(layer_data)
                    new_fov = cast(Points, new_fov)
                elif isinstance(layer, Labels) :
                    new_fov = self.viewer.add_labels(layer_data)
                    new_fov = cast(Labels, new_fov)
                else :
                    print("Unsported layer type for mapping, ignoring mapping.")
                    break
                
                
                new_image.translate = translate_coords[location_index]
                location_index += 1
            
            self.viewer.layers.selection.clear()
            QTimer.singleShot(0, lambda: self.viewer.layers.remove(layer)) #Delays removing after viewer internal state is updated.
            

        finally :
            self.start_listening()
