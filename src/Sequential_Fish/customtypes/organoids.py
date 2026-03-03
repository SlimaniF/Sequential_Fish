from pydantic import BaseModel
from abc import ABC, abstractmethod
from typing import Iterator, TypedDict

class Location(BaseModel) :
    """
    Expected data for each location. This class is used to enforce data structure passed to OrganoidLocations
    """
    name : str
    x : float
    y : float
    z_begin : float
    z_end : float

class Location_dict(TypedDict) :
    x : float
    y : float
    z_begin : float
    z_end : float

class OrganoidLocations :
    """
    Target class for organoids location data, validate that loaded data can be exploited by viewer module. 
    If you construct a class to load your own data structure it needs to output data as an instance of this class.
    """

    def __init__(
        self,
        locations : list[Location]
        ) -> None:
        self._locations = {}
        for loc in locations :
            self._locations[loc.name] = Location.model_validate(loc)

    def __getitem__(
        self,
        location : str
        ) -> Location_dict :
        loc = self._locations[location]
        return {"x" : loc.x, "y" : loc.y, "z_begin" : loc.z_begin, "z_end" : loc.z_end}

    def __iter__(self) -> Iterator[Location_dict] :
        return iter([self[name] for name in self._locations.keys()])
    
    def __contains__(self, name) -> bool :
        return name in self._locations

    def validate(self, locations : list[str]) :
        """
        Returns True if all locations are in OrganoidLocations.
        """
        return all([location in self for location in locations])

class LocationsExport(ABC) :

    @property
    @abstractmethod
    def name(self) -> str:
        """
        Every subclass must define a unique name.
        """

    @abstractmethod
    def get_locations(self, json_data) -> OrganoidLocations :
        """
        This method should return an OrganoidLocations instance enforcing your data can be used by viewer module.
        """
    @abstractmethod
    def set_name(self) :
        """
        set name attribute for LocationExport
        """


def load_organoid_locations(json_data, location_export : str) -> OrganoidLocations:
    if location_export in _LOCATIONSEXPORTS.keys() :
        return _LOCATIONSEXPORTS[location_export].get_organoid_locations(json_data)
    else :
        raise NotImplementedError(f"There is no implementation of json data matching name {location_export}. Available implementations are {_LOCATIONSEXPORTS}")


### LocationsExport
# 
# If you are using a differently shaped data for organoids locations on your microscope you will need to develop
# a class able to retrieve data for x,y and z locations. Then all you need to do is set a name and get_locations methods 
# yielding above defined class. Add the @register_locations_exports to make it usable in load_organoid_locations.

_LOCATIONSEXPORTS = {}
def register_locations_exports(cls) :
    if issubclass(cls, LocationsExport) :
        _LOCATIONSEXPORTS[cls.name] = cls

@register_locations_exports
class StageSelectionExport(LocationsExport) :
    """
    Model class for structured load of organoids locations file produced from Vutara SRX 7.6.07. Default configuration for this package.
    """