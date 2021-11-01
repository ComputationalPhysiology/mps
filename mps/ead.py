from typing import List
from typing import Optional
from typing import Tuple
from typing import Union

import numpy as np


try:
    import matplotlib.pyplot as plt

    has_mpl = True
except ImportError:
    has_mpl = False


def detect_ead(
    y: Union[List[float], np.ndarray],
    sigma: float = 3,
    prominence_level: float = 0.04,
) -> Tuple[bool, Optional[int]]:
    """Detect (Early afterdepolarizations) EADs
    based on peak prominence.

    Arguments
    ---------
    y : np.ndarray
        The signal that you want to detect EADs
    sigma : float
        Standard deviation in the gaussian smoothing kernal
        Default: 3.0
    prominence_level: float
        How prominent a peak should be in order to be
        characterized as an EAD. This value shold be
        between 0 and 1, with a greater value being
        more prominent. Defaulta: 0.04

    Notes
    -----
    Given a signal :math:`y` we want to determine wether we have
    an EAD present in the signal. `EADs <https://en.wikipedia.org/wiki/Afterdepolarization>`_
    are abnormal depolarizations happening after the upstroke in an action potential.

    We assume that an EAD occurs betweeen the maximum value of the signal
    (i.e the peak) and the next minimum value (i.e when the signal is at rest)

    To remove noisy patterns we first smooth the signal
    with a gaussian filter. Then we take out only the part
    of the signal that is between its maximum and the next
    minimum values. Then we find the peaks with a
    `Topographic Prominence <https://en.wikipedia.org/wiki/Topographic_prominence>`_
    greather than the given prominence level

    Returns
    -------
    bool:
        Flag indicating if an EAD is found or not
    int or None:
        Index where we found the EAD. If no EAD is found then
        this will be None. I more than one peaks are found then
        only the first will be returned.

    """
    from scipy.ndimage import gaussian_filter1d
    from scipy.signal import find_peaks

    y = np.array(y)
    idx_max = int(np.argmax(y))
    idx_min = idx_max + int(np.argmin(y[idx_max:]))

    y_tmp = y[idx_max:idx_min] - y[idx_min]
    if len(y_tmp) == 0:
        return False, None

    y_smooth = gaussian_filter1d(y_tmp / np.max(y_tmp), sigma)
    peaks, props = find_peaks(y_smooth, prominence=prominence_level)

    return len(peaks) > 0, None if len(peaks) == 0 else int(peaks[0] + idx_max)


def analyze_eads(
    chopped_data: List[List[float]],
    chopped_times: List[List[float]],
    sigma: float = 3,
    prominence_threshold: float = 0.04,
    plot=True,
    fname: str = "",
) -> int:
    """
    Loop over all beats and check for EADs.

    Arguments
    ---------
    chopped_data : list
        List of the amplitude of each beat
    chopped_data : list
        List of the time stamps of each beat
    sigma: float
        Standard deviation in the gaussian smoothing kernal used for
        EAD detection. Default: 3.0
    prominence_threshold: float
        How prominent a peak should be in order to be
        characterized as an EAD. This value shold be
        between 0 and 1, with a greater value being
        more prominent. Defaulta: 0.04
    plot : bool
        If True we plot the beats with the potential EAD marked
        as a red dot. Default: True
    fname : str
        Path to file that you want to save the plot.

    Returns
    -------
    int:
        The number of EADs
    """

    num_eads = 0
    peak_inds = {}
    for i, (c, t) in enumerate(zip(chopped_data, chopped_times)):
        has_ead, peak_index = detect_ead(
            c,
            sigma=sigma,
            prominence_level=prominence_threshold,
        )
        if has_ead and peak_index is not None:
            num_eads += 1
            peak_inds[i] = peak_index

    if plot and has_mpl:
        plot_ead(fname, chopped_data, chopped_times, peak_inds)

    return num_eads


def plot_ead(fname, chopped_data, chopped_times, peak_inds):

    fig, ax = plt.subplots()
    num_eads = 0
    for i, (c, t) in enumerate(zip(chopped_data, chopped_times)):
        ax.plot(t, c)
        if i in peak_inds:
            ax.plot([t[peak_inds[i]]], [c[peak_inds[i]]], "ro")
    ax.set_title(f"EADs are marked with red dots. Found EADs in {num_eads} beats")
    if fname != "":
        fig.savefig(fname)
    plt.close()
