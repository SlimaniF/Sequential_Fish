from torch.cuda import is_available
from ..customtypes import NapariWidget
from .utils import update_layer_from_LayerDataTuple

import torch

from typing import cast
import numpy as np

from napari import Viewer
from napari.layers import Image
from napari.types import LayerDataTuple
from napari.qt.threading import create_worker


from magicgui import magicgui
import cellpose.models as models
from cellpose.core import use_gpu

available_models = ["cpsam"] + models.get_user_models()
gpu_is_available = use_gpu()


def estimate_batch_size() :
    
    if torch.cuda.is_available() :
        device = torch.cuda.current_device()
        total_memory = torch.cuda.get_device_properties(device).total_memory
        print(f"  Total Memory: {total_memory / (1024**3):.2f} GB")

        if (total_memory / (1024**3)) > 8 :
            batch_size = 16
        else : 
            batch_size = 8

    else :
        batch_size = 2


    return batch_size


class SegmentationWidget(NapariWidget) :
    
    def __init__(self, *, viewer : Viewer):
        super().__init__()
        self.viewer = viewer
        self._wrap_with_threading()

    def _wrap_with_threading(self):
        """
        Wraps the widget's function with napari's threading.
        """
        if hasattr(self.widget, 'native'):
            # Get the underlying function from the magicgui widget
            original_func = self.widget._function
            def threaded_func(*args, **kwargs):
                worker = create_worker(original_func, *args, **kwargs)
                worker.returned.connect(lambda result: update_layer_from_LayerDataTuple(self.viewer, result))
                worker.start()
                return worker
            # Replace the widget's function with the threaded version
            self.widget._function = threaded_func



_SEGMETATION_WIDGETS : 'list[NapariWidget]' = []
def register_segmentation_widget(cls) :
    _SEGMETATION_WIDGETS.append(cls)
    return cls

def initiate_segmentation_widgets(
        viewer : Viewer,
        voxel_size : tuple,
) -> tuple[list[NapariWidget], list[NapariWidget]]:
    widget_list = []
    linked_widgets = []

    for cls in _SEGMETATION_WIDGETS :
        instance = cls(voxel_size=voxel_size, viewer=viewer)
        if hasattr(instance,"enabled") :
            if instance.enabled :
                widget_list.extend(instance.get_widgets())
        else :
            widget_list.extend(instance.get_widgets())

        if hasattr(instance, "update") and callable(getattr(instance, "update")):
            linked_widgets.append(instance)

    return widget_list, linked_widgets


@register_segmentation_widget
class SegmentationTester(SegmentationWidget) :
    def __init__(self, voxel_size : tuple[int,int,int], viewer : Viewer):
        self.voxel_size = voxel_size
        self.anisotropy : float = voxel_size[0] / voxel_size[1]
        self.model : models.CellposeModel | None = None
        self.model_name = "cpsam"
        self.batch_size = estimate_batch_size()

        super().__init__(viewer=viewer)

    def _create_widget(self):
        
        @magicgui(
            use_gpu = {"widget_type" : "CheckBox", "enabled" : gpu_is_available},
            model = {"choices" : available_models},
            anisotropy = {"enabled" : False}, # computed from voxel size, just to show user
            min_size = {"max" : 2**32}
        )
        def segment_image(
            image : Image,
            model : str = "cpsam",
            diameter : int = 30,
            anisotropy : float = self.anisotropy,
            cellprob_threshold : float = 0.,
            flow_threshold : float = 0.4,
            min_size : int = 15,
            do_3D : bool = False,
            use_gpu : bool = gpu_is_available,
        ) :

            # Init model
            if model != self.model_name or self.model is None:
                self.model = models.CellposeModel(
                    gpu=use_gpu,
                    pretrained_model=model
                )
                self.model_name = model

            image_data = cast(np.ndarray, image.data)
            #Getting number of z slice:
            if image_data.ndim == 4 :
                z_slice = image_data.shape[1]
            else :
                z_slice = image_data.shape[0]

            # if segmentation in 2D do mean projection
            if image_data.ndim == 4 and not do_3D: # case cycles-z-yx
                image_data = np.mean(image_data, axis=1)
            elif image_data.ndim == 3 and not do_3D : #case z-yx
                image_data = np.mean(image_data, axis=0)

            if image_data.ndim == 3 and not do_3D : # case of several cycles to segment
                image_data = [image for image in image_data]
            elif image_data.ndim == 4 and do_3D :
                image_data = [image for image in image_data]
            elif image_data.ndim == 3 and do_3D : 
                pass
            elif image_data.ndim == 2 and not do_3D : #only one cycle to segment
                pass
            else :
                raise AssertionError("Unforseen dimension number in segmentation")

            mask, *_ =self.model.eval(
                image_data,
                diameter= diameter,
                cellprob_threshold=cellprob_threshold,
                flow_threshold=flow_threshold,
                min_size=min_size,
                do_3D=do_3D,
                z_axis= 0 if do_3D else None,
                anisotropy=anisotropy if do_3D else None,
                batch_size= self.batch_size if use_gpu else 2 # Recommendation from cellpose doc can be increase to 16 for large memory GPU
            )
            
            if isinstance(mask,list) :
                mask = np.stack(mask)
            else : raise AssertionError(f"Expected list type got {type(mask)}")
            
            if mask.ndim == 4 and do_3D : #Segmented a list of 3D
                voxel_size = (1,) + self.voxel_size
            
            elif mask.ndim == 3 and not do_3D : #Segmented a list of 2D
                voxel_size = (1,) + self.voxel_size
                mask = mask.reshape(mask.shape[0],1,*mask.shape[1:])
                mask = np.repeat(mask,z_slice, axis=1)

            else :
                raise AssertionError("Unforseen dimension number in segmentation")
            

            layer_name = str(image.name) + "_segmentation"
            res = LayerDataTuple((
                mask,
                {"scale" : voxel_size, "name" : layer_name, "blending" : "additive"},
                'Labels'
            ))

            return res
            
        return segment_image