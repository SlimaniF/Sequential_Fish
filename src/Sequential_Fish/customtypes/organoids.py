import json
from pathlib import Path
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
    """
    If you are using a differently shaped data for organoids locations on your microscope you will need to develop
    a class able to retrieve data for x,y and z locations. Then all you need to do is set a name and get_locations methods 
    yielding above defined class. Add the @register_locations_exports to make it usable in load_organoid_locations.
    """
    name : str # All LocationsExport class must define a unique class attribute "name"

    @abstractmethod
    def get_locations(self) -> OrganoidLocations :
        """
        This method should return an OrganoidLocations instance enforcing your data can be used by viewer module.
        """


def load_organoid_locations(json_data, location_export : str) -> OrganoidLocations:
    if location_export in _LOCATIONSEXPORTS.keys() :
        return _LOCATIONSEXPORTS[location_export].get_organoid_locations(json_data)
    else :
        raise NotImplementedError(f"There is no implementation of json data matching name {location_export}. Available implementations are {_LOCATIONSEXPORTS}")


### LocationsExport
# 


_LOCATIONSEXPORTS = {}
def register_locations_exports(cls) :
    if issubclass(cls, LocationsExport) :
        if cls.name in _LOCATIONSEXPORTS.keys() :
            raise LocationDataNameError(f"Cannot register {cls.name} as location export as another class is already registered with same name")
        else :
            _LOCATIONSEXPORTS[cls.name] = cls

class LocationsDataStructureError(ValueError) :
    """
    Error raised when trying to open an organoid location file with a method that is not fit for the data structure found, indicating user that he needs to adjust his metadata format or 
    implement a new LocationsExport subclass matching his pattern.
    """

class LocationDataNameError(ValueError) :
    """
    Error raised when two subclass of LocationsExport share the same class attribute 'name'.
    """


@register_locations_exports
class VutaraLocationsExport(LocationsExport) :
    """
    Model class for structured load of organoids locations file produced from Vutara SRX 7.6.07. Default configuration for this package.
    """

    name = "Vutara_SRX_CaptureLocations"

    def __init__(self, location_fullpath : str | Path) -> None:
        super().__init__()

        with open(location_fullpath, "r") as location_rawdata :
            self.rawdata = json.load(location_rawdata)
        self.location_list : list[Location] = self.validate_structure()

    def validate_structure(self) :
        """
        Validates that the raw data has the expected structure.
        Raises LocationsDataStructureError if validation fails.
        """
        # Step 1: Check if rawdata is a dictionary
        if not isinstance(self.rawdata, dict):
            raise LocationsDataStructureError("Raw data must be a dictionary.")

        # Step 2: Check if "value" key exists
        if "value" not in self.rawdata:
            raise LocationsDataStructureError('Raw data must contain a "value" key.')

        # Step 3: Check if "CaptureLocations" key exists in "value"
        if "CaptureLocations" not in self.rawdata["value"]:
            raise LocationsDataStructureError('"value" must contain a "CaptureLocations" key.')

        # Step 4: Check if "value" key exists in "CaptureLocations"
        if "value" not in self.rawdata["value"]["CaptureLocations"]:
            raise LocationsDataStructureError('"CaptureLocations" must contain a "value" key.')

        # Step 5: Check if "value" is a list
        capture_locations = self.rawdata["value"]["CaptureLocations"]["value"]
        if not isinstance(capture_locations, list):
            raise LocationsDataStructureError('"CaptureLocations.value" must be a list.')

        # Step 6: Validate each location in the list
        try:
            Location_list = [Location.model_validate(loc) for loc in capture_locations]
        except Exception as e:
            raise LocationsDataStructureError(f"Could not validate data structure from Locations dicts found.\n{e}") from e
        else :
            return Location_list

    def get_locations(self) -> OrganoidLocations:
        return OrganoidLocations(self.location_list)