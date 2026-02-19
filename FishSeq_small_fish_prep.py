RUN_PATH = "/media/SSD_floricslimani/Fish_seq/Davide/2026-02-05 - M-FISH - HeLa - 13 genes_Run1/"


import os
import numpy as np
from smfishtools.utils import open_image
from tqdm import tqdm
import bigfish.stack as stack

im_path = RUN_PATH + "/Location-01/img000_000_000000_0000000000.ome.tif"
os.makedirs(RUN_PATH + "/small_fish_sample/",exist_ok=True)

image_stack = open_image(im_path)

cycle = 0
print(image_stack.shape)
for image in tqdm(image_stack) :
    image = np.moveaxis(image,[0,1,2,3], [1,0,2,3])
    print("reodered shape : ", image.shape)
    assert image.shape[0] <= 5
    for index, color in enumerate(image) :
        stack.save_image(image=color, path= RUN_PATH + "/small_fish_sample/cycle_{0}_color_{1}.tiff".format(cycle,index))
    cycle +=1
