import numpy as np
from skimage import filters
from skimage import morphology
from skimage import segmentation


def detect_tissue(
    frame: np.ndarray,
    sigma: float = 20.0,
    min_size: int = 400,
    fill_holes: bool = False,
) -> np.ndarray:
    smooth = filters.gaussian(frame, sigma=sigma)

    thresh_value = filters.threshold_otsu(smooth)
    thresh = smooth > thresh_value

    im = morphology.remove_small_objects(thresh, min_size=min_size)

    if fill_holes:
        from scipy import ndimage as ndi

        fill = ndi.binary_fill_holes(im)
        clear = segmentation.clear_border(fill)
        erode = morphology.binary_erosion(clear)
        im = morphology.binary_dilation(erode)

    return im
