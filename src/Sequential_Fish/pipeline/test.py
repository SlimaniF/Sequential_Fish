import os
import ome_types
from ome_types.model import TiffData
import tifffile
import numpy as np
from tqdm import tqdm

RUN_PATH = "/media/SSD_floricslimani/Fish_seq/POLR2/2026-04-09 - HeLa-Puro_POLR2_Run18/FISH_Z-stacks"
output_path = os.path.join(RUN_PATH, "Location-01")
filename = "img000_000_000000_0000000000.ome.tif"
os.makedirs(output_path, exist_ok=True)


print("Opening image...")
file_path_main = os.path.join(RUN_PATH, "Location-01", "img001_000_000000_0000000000.ome.tif")
file_path_other = os.path.join(RUN_PATH, "Location-01", "img000_000_000000_0000000000.ome.tif")

with tifffile.TiffFile(file_path_other) as tif :
    metadata = ome_types.from_xml(tif.pages[0].description)
    axes_order = tif.series[0].axes
    photometric = tif.pages[0].photometric
    compression = tif.pages[0].compression
    dtype = tif.pages[0].dtype

with tifffile.TiffFile(file_path_main) as tif :
    image_stack = tif.asarray()



print(f"\nFILENAME : {filename}")
print("images len : ", len(metadata.images))
if len(metadata.images) > 0 :
    print("Image[0] : ", dict(metadata.images[0]).keys())
    print("Image[0] pixels datablocks : ", metadata.images[0].pixels.tiff_data_blocks)
    print("Image[0] pixels keys : ", dict(metadata.images[0].pixels).keys())
    print("Image[0] size_t : ", metadata.images[0].pixels.size_t)
    print("Image[0] size_z : ", metadata.images[0].pixels.size_z)
    print("Image[0] size_c : ", metadata.images[0].pixels.size_c)
    print("Image[0] size_y : ", metadata.images[0].pixels.size_y)
    print("Image[0] size_x : ", metadata.images[0].pixels.size_x)
    print("Image[0] dimension_order : ", metadata.images[0].pixels.dimension_order)
    print("Image[0] physical size z : ", metadata.images[0].pixels.physical_size_z, metadata.images[0].pixels.physical_size_z_unit)
    print("Image[0] physical size y : ", metadata.images[0].pixels.physical_size_y, metadata.images[0].pixels.physical_size_y_unit)
    print("Image[0] physical size x : ", metadata.images[0].pixels.physical_size_x, metadata.images[0].pixels.physical_size_x_unit)




if cycle_number == 16 :
    print(f"found {cycle_number} cycles. Adding an empty cycle")
    new_filename = filename.replace("img000", f"img0{cycle_number}")
    with tifffile.TiffWriter(os.path.join(output_path, new_filename), bigtiff=True) as tif :
        
        metadata.images[0].pixels.tiff_data_blocks.append(
            TiffData(
                first_c = 0,
                first_z = 0,
                first_t = cycle_number +1,
                plane_count = len(empty_cycle),
                uuid=TiffData.UUID(value=tiff_uuid, file_name=filename)
            )
        )

        for plane in tqdm(empty_cycle, desc="Writing planes") :
            tif.write(
                plane,
                photometric=photometric,
                metadata={'axes' : axes_order},
                dtype=dtype,
                compression=compression
            )

else :
    print(f"found {cycle_number} cycles. Leaving")