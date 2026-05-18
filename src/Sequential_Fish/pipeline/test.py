import os
import ome_types
from ome_types.model import TiffData
import tifffile
import numpy as np
from tqdm import tqdm
import uuid

RUN_PATH = "/home/floric/Documents/SeqfFish/2026-03-24 - HeLa_POLR2_Run17_tiff/FISH_Z-stacks/"
output_path = os.path.join(RUN_PATH, "output")
filename = "img000_000_000000_0000000000.ome.tif"
os.makedirs(output_path, exist_ok=True)


print("Opening image...")
file_path_main = os.path.join(RUN_PATH, "Location-01", "img000_000_000000_0000000000.ome.tif")
file_path_other = os.path.join(RUN_PATH, "Location-01", "img001_000_000000_0000000000.ome.tif")

with tifffile.TiffFile(file_path_other) as tif :
    axes_order = tif.series[0].axes
    photometric = tif.pages[0].photometric
    compression = tif.pages[0].compression
    dtype = tif.pages[0].dtype

    print("\nOther")
    print("axes order : ",axes_order)
    print("photometric : ",photometric)
    print("compression : ",compression)
    print("dtype : ",dtype)

with tifffile.TiffFile(file_path_main) as tif :
    metadata = ome_types.from_xml(tif.pages[0].description)
    main_axes_order = tif.series[0].axes
    main_dtype = tif.pages[0].dtype
    main_photometric = tif.pages[0].photometric
    main_compression = tif.pages[0].compression
    image_number = len(tif.pages)

    print("\nMain")
    print("axes order : ",main_axes_order)
    print("photometric : ",main_photometric)
    print("compression : ",main_compression)
    print("dtype : ",main_dtype)
    print("image number : ", image_number)

    image_stack = tif.asarray()

cycle_number, *stack_shape = image_stack.shape
print("stack shape : ", stack_shape)

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
    empty_cycle = np.zeros(shape=(image_number, stack_shape[-2], stack_shape[-1]), dtype= image_stack.dtype)
    new_filename = filename.replace("img000", f"img0{cycle_number}")
    with tifffile.TiffWriter(os.path.join(output_path, new_filename), bigtiff=True) as tif :
        
        tiff_uuid = f'urn:uuid:{uuid.uuid4()}'
        metadata.images[0].pixels.tiff_data_blocks.append(
            TiffData(
                first_t = cycle_number +1,
                plane_count = len(empty_cycle),
                uuid=TiffData.UUID(value=tiff_uuid, file_name=filename)
            )
        )

        print(f"Writing new_cycle of shape : {empty_cycle.shape} with order : {axes_order}")
        tif.write(
                empty_cycle,
                photometric=photometric,
                metadata={'axes' : "ZYX"},
                dtype=dtype,
                compression=compression
            )

    with tifffile.TiffWriter(os.path.join(output_path, filename), bigtiff=True) :
        print(f"Writing main tiff of shape : {image_stack[0].shape} with order : {main_axes_order}")
        tif.write(
            image_stack[0],
            photometric=main_photometric,
            metadata={"axes" : main_axes_order},
            dtype = main_dtype,
            compression = main_compression,
            description = ome_types.to_xml(metadata)
        )


else :
    print(f"found {cycle_number} cycles. Leaving")