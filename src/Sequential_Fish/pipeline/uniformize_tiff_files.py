"""
Aim : Sanitize experiments where cycle number is uneven through locations. Such a case can sucess if acquisition computer crashes during experiement..;
"""

WARNING = """
This script aims at sanitizing experiments where cycle number is uneven through locations. Such a case can sucess if acquisition computer crashes during experiement...
"""

from .input import UnevenCycleNumber, InputIntegrityError
from .input import main as launch_input
from ..settings import get_settings, PipelineParameters

import os
import uuid
import re
import logging

import tifffile
import ome_types
from ome_types import OME
from typing import cast
from ome_types import from_xml
from tqdm import tqdm
import numpy as np
from pandas import read_excel, read_csv

def main(run_path) :
    logging.info("\n\nStarting uniformization script.\n")
    
    try :
        print("Launching input script...")
        logging.info("Run input script.")
        launch_input(run_path)
    
    except InputIntegrityError as error :

        print(str(error))

        if isinstance(error, UnevenCycleNumber) :
            expected_cycle_number = error.expected_cycle
            found_cycle_number = error.found_cycle
        else :
            expected_cycle_number = get_expected_cycle_number(run_path)
            found_cycle_number = "ambiguous"

        print(f"Script input returns {type(error)}.    \nexpected cycle : {expected_cycle_number}  \nfound : {found_cycle_number}")
        logging.info(f"Script input returns {type(error)}.    \nexpected cycle : {expected_cycle_number}  \nfound : {found_cycle_number}")
        
        while True :
            answ = input(f"Uniformize all locations to have {expected_cycle_number} cycles ?(y/n) This will modify images in run folder, make sure you keep a copy of your raw data.\n")

            if answ == "y" :
                logging.info(f"Uniformize all locations to have {expected_cycle_number} cycles.")

                pipeline_parameters = get_settings(run_path, settings_name="pipeline")
                pipeline_parameters = cast(PipelineParameters, pipeline_parameters)
                FISH_FOLDER = pipeline_parameters.FISH_FOLDER
                location_list = sorted(os.listdir(os.path.join(run_path,FISH_FOLDER)))
                fish_path = (run_path + "/{0}/{1}/".format(FISH_FOLDER, location)for location in location_list)
                tifffile_list = (os.path.join(folder, sorted(os.listdir(folder))[0]) for folder in fish_path)
                
                for masterfile_fullpath in tqdm(tifffile_list, desc="Uniformizing locations", total= len(location_list)) :

                    location = os.path.dirname(masterfile_fullpath)
                    logging.info(f"Location : {location}")
                    masterfilename = os.path.basename(masterfile_fullpath)
                    found_cycle_number = get_cycle_number_from_metadata(masterfile_fullpath)
                    
                    logging.info(f"Found {found_cycle_number} in metadata, expected {expected_cycle_number}.")
                    if found_cycle_number > expected_cycle_number :
                        logging.info(f"Removing cycle(s).")
                        
                        remove_page(
                            location_path= location,
                            masterfilename=masterfilename,
                            expected_cycle_number=expected_cycle_number
                        )
                    
                    elif found_cycle_number < expected_cycle_number :
                        logging.info(f"Adding cycle(s).")
                        
                        add_page(
                            location_path=location,
                            masterfilename=masterfilename,
                            expected_cycle_number=expected_cycle_number
                        )
                    
                    #Else continue

                break

            elif answ == "n" :
                logging.info(f"Abborted.")
                quit()
    else :
        print("No error when running input step, you can proceed to next steps.")


def get_cycle_number_from_metadata(
    masterfile_fullpath : str,
    ) -> int:
    if ".ome." in masterfile_fullpath :
        fish_im = None
        with tifffile.TiffFile(masterfile_fullpath) as main_tif :
            #Axes map
            description = getattr(main_tif.pages[0], "description")
            metadata = from_xml(description)

            if not metadata.binary_only is None :
                raise ValueError(f"{os.path.basename(masterfile_fullpath)} is not tiff serie masterfile.")

            return metadata.images[0].pixels.size_t


    return 0


def add_page(
    location_path : str,
    masterfilename : str,
    expected_cycle_number : int
    ) :

    with tifffile.TiffFile(os.path.join(location_path,masterfilename)) as tif :

        main_metadata = from_xml(
            getattr(tif.pages[0], "description")
        )
        dtype = tif.pages[0].dtype
        main_uuid = cast(str,main_metadata.uuid)
        photometric = getattr(tif.pages[0], "photometric")
        pages_number = len(tif.pages)
        compression = getattr(tif.pages[0], "compression")
        full_stack = tif.asarray()
    
    found_cycle_number, *_,Y,X = full_stack.shape
    empty_stack_to_add_number = expected_cycle_number - found_cycle_number
    assert empty_stack_to_add_number > 0

    new_tiff_data_block = []
    empty_stack = np.zeros(shape=(pages_number, Y,X), dtype= full_stack.dtype)
    for i in range(empty_stack_to_add_number) :
        
        #Cycle metadata and filename
        new_cycle_metadata = OME()
        new_cycle_uuid = f'urn:uuid:{uuid.uuid4()}'
        new_cycle_metadata.binary_only = OME.BinaryOnly(metadata_file=masterfilename, uuid=main_uuid)
        new_filename = re.sub(r"\d{3}", f"{found_cycle_number + i}".rjust(3,"0"), masterfilename,count=1)
        
        #TiffData building
        new_tiff_data_block.append(
            ome_types.model.TiffData(
                first_t = found_cycle_number + i,
                plane_count = pages_number,
                uuid=ome_types.model.TiffData.UUID(value=new_cycle_uuid, file_name=new_filename)
            )
        )

        #Updating mainfilemetadata
        main_metadata.images[0].pixels.size_t += 1

        #Writing new stacks
        with tifffile.TiffWriter(os.path.join(location_path, new_filename), bigtiff=True) as tif :
            tif.write(
                empty_stack,
                photometric=photometric,
                dtype=dtype,
                description=ome_types.to_xml(new_cycle_metadata),
                compression=compression
            )

    #Re-writing main stack to update main_metadata
    with tifffile.TiffWriter(os.path.join(location_path, masterfilename), bigtiff=True, shaped=False) as tif:

        #Updating main metadata to link new cycles
        main_metadata.images[0].pixels.tiff_data_blocks.extend(new_tiff_data_block)
        tif.write(
            full_stack[0],
            photometric=photometric,
            dtype = dtype,
            compression = compression,
            description = ome_types.to_xml(main_metadata)
        )


def remove_page(
    location_path : str,
    masterfilename : str,
    expected_cycle_number : int
    ) :

    with tifffile.TiffFile(os.path.join(location_path,masterfilename)) as tif :

        main_metadata = from_xml(
            getattr(tif.pages[0], "description")
        )
        dtype = tif.pages[0].dtype
        photometric = getattr(tif.pages[0], "photometric")
        compression = getattr(tif.pages[0], "compression")
        full_stack = tif.asarray()
    
    found_cycle_number, *_,Y,X = full_stack.shape
    stack_to_remove_number = found_cycle_number - expected_cycle_number
    assert stack_to_remove_number > 0

    for i in range(stack_to_remove_number) :

        #Updating mainfilemetadata
        main_metadata.images[0].pixels.tiff_data_blocks.pop()
        main_metadata.images[0].pixels.size_t -= 1

    #Re-writing main stack to update main_metadata
    with tifffile.TiffWriter(os.path.join(location_path, masterfilename), bigtiff=True, shaped=False) as tif:
        tif.write(
            full_stack[0],
            photometric=photometric,
            dtype = dtype,
            compression = compression,
            description = ome_types.to_xml(main_metadata)
        )

def get_expected_cycle_number(run_path : str) :

    pipeline_parameters = cast(PipelineParameters, get_settings(run_path=run_path, settings_name="pipeline"))
    map_filename = os.path.join(run_path, pipeline_parameters.MAP_FILENAME)
    if map_filename.endswith(".xlsx") :
        return len(read_excel(map_filename))
    elif map_filename.endswith(".csv") :
        return len(read_csv(map_filename))
    else :
        raise NotImplementedError("Expected xlsx or csv extension for cycle map.")