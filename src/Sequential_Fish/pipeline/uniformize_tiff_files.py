"""
Aim : Sanitize experiments where cycle number is uneven through locations. Such a case can sucess if acquisition computer crashes during experiement..;
"""

WARNING = """
This script aims at sanitizing experiments where cycle number is uneven through locations. Such a case can sucess if acquisition computer crashes during experiement...
"""

from .input import UnevenCycleNumber, AxisMap
from .input import map_from_metadata, auto_map_channels, get_ordered_shape_from_metadata, get_memory_shape_from_metadata
from .input import main as launch_input
from ..settings import get_settings, PipelineParameters
from ..tools._folder_integrity import assert_run_folder_integrity

import os
import tifffile
from ome_types import OME
from typing import cast
from ome_types import from_xml
from tqdm import tqdm
import numpy as np

def main(run_path) :
    
    try :
        print("Launching input script...")
        launch_input(run_path)
    
    except UnevenCycleNumber as error :
        expected_cycle_number = error.expected_cycle
        found_cycle_number = error.found_cycle

        print(f"Script input returns uneven cycle error.    \nexpected cycle : {expected_cycle_number}  \nfound : {found_cycle_number}")
        
        while True :
            answ = input(f"Uniformize all locations to have {expected_cycle_number} cycles ?(y/n) This will modify images in run folder, make sure you keep a copy of your raw data.\n")

            if answ == "y" :
                pipeline_parameters = get_settings(run_path, settings_name="pipeline")
                pipeline_parameters = cast(PipelineParameters, pipeline_parameters)
                FISH_FOLDER = pipeline_parameters.FISH_FOLDER
                color_number = len(pipeline_parameters.GENES_NAMES_KEY)
                has_beads = not pipeline_parameters.BEAD_CHANNEl is None
                location_list = get_location_list(FISH_FOLDER, run_path)
                fish_path = (run_path + "/{0}/{1}/".format(FISH_FOLDER, location)for location in location_list)
                tifffile_list = (os.path.join(folder, sorted(os.listdir(folder))[0]) for folder in fish_path)
                
                for file in tqdm(tifffile_list, desc="Uniformizing locations", total= len(location_list)) :
                    found_cycle_number = get_cycle_number_from_metada(file, color_number=color_number, expected_cycle_number=expected_cycle_number, has_bead_channel=has_beads)
                    
                    if found_cycle_number > expected_cycle_number :
                        os.makedirs(os.path.dirname(file.replace(FISH_FOLDER,"save_modified")),exist_ok=True)
                        print("removing. Found : ", found_cycle_number)
                        remove_page(
                            input_path=file,
                            output_path=file.replace(FISH_FOLDER,"save_modified"),
                            cycle_number=expected_cycle_number
                        )
                    elif found_cycle_number < expected_cycle_number :
                        os.makedirs(os.path.dirname(file.replace(FISH_FOLDER,"save_modified")),exist_ok=True)
                        print("Adding. Found : ", found_cycle_number)
                        add_page(
                            input_path=file,
                            output_path=file.replace(FISH_FOLDER,"save_modified"),
                            cycle_number=expected_cycle_number
                        )
                break

            elif answ == "n" :
                quit()
    else :
        print("No error when running input step, you can proceed to next steps.")


def get_location_list(fish_folder, run_path) -> list[str]:
    
    file_dict = assert_run_folder_integrity(
    run_path=run_path,
    fish_folder=fish_folder,
    nucleus_folder=fish_folder
    )
    location_list = list(file_dict.keys())
    location_list.sort()

    return location_list

def get_cycle_number_from_metada(
    file_fullpath : str,
    color_number : int,
    expected_cycle_number : int,
    has_bead_channel : bool,
    ) -> int:
    if ".ome." in file_fullpath :
        fish_im = None
        with tifffile.TiffFile(file_fullpath) as main_tif :
            #Axes map
            fish_map = map_from_metadata(main_tif.series[0].axes)
            if fish_map is None : #Could not infer from metadata
                print("Could not infer axes map from metadata opening image.")
                fish_im = main_tif.asarray()
                fish_map = cast(AxisMap, auto_map_channels(fish_im, color_number=color_number, cycle_number=expected_cycle_number, has_bead_channel=has_bead_channel))
            
            #Shape
            fish_shape = None
            fish_reodered_shape = None
            #Infer from metadata
            if not main_tif.ome_metadata is None :
                ome_metadata = from_xml(main_tif.ome_metadata)
                fish_shape = get_memory_shape_from_metadata(ome_metadata, fish_map)
                
                if not fish_shape is None :
                    found_cycle_number = fish_shape[fish_map['cycles']]
                    return found_cycle_number

    return 0

def remove_page(input_path :str, output_path :str, cycle_number : int):
    # Read the OME-TIFF file
    with tifffile.TiffFile(input_path) as tif:
        # Get the OME metadata
        ome_metadata = OME.from_xml(tif.ome_metadata)
        images = ome_metadata.images

        # Get the image data and metadata for the first (and only) image
        print("images : ",len(images))
        image = images[0]
        pixels = image.pixels

        # Read all pages (3D stacks)
        all_pages = tif.pages
        num_pages = len(all_pages)
        print("num_pages : ",num_pages)

        pages_to_keep = all_pages[:cycle_number] 

        # Read the data for the pages to keep
        data_to_keep = [page.asarray() for page in pages_to_keep]

        # Write the new OME-TIFF file
        # with tifffile.TiffWriter(output_path, bigtiff=True) as tif_out:
        #     for i, data in enumerate(data_to_keep):
        #         # Write each page as a new subfile
        #         tif_out.write(
        #             data,
        #             photometric='minisblack',
        #             description=ome_metadata.to_xml().encode('utf-8')
        #         )

def add_page(input_path :str, output_path :str, cycle_number : int) :
    
    with tifffile.TiffFile(input_path) as tif:
        # Get the OME metadata
        ome_metadata = OME.from_xml(tif.ome_metadata)
        # Get the first (and only) image's metadata
        image = ome_metadata.images[0]
        pixels = image.pixels

        # Read all pages
        all_pages = tif.pages
        # Get the shape and dtype of the first page
        first_page = all_pages[0]
        page_number = len(all_pages)
        shape = first_page.shape
        dtype = first_page.dtype

        # Create an empty page with the same shape and dtype
        empty_page = np.zeros(shape, dtype=dtype)

        # Write all original pages + the new empty page to a new file
        with tifffile.TiffWriter(output_path, bigtiff=True) as tif_out:
            for page in all_pages:
                print(page.name)
                tif_out.write(
                    page.asarray(),
                    photometric=page.photometric,
                    description=tif.ome_metadata.encode('utf-8')
                )
            # Add the empty page
            for i in range(cycle_number - page_number) :
                tif_out.write(
                    empty_page,
                    photometric=first_page.photometric,
                    description=tif.ome_metadata.encode('utf-8')
                )