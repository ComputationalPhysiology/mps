#!/usr/bin/env python3
__author__ = "Henrik Finsberg (henriknf@simula.no), 2017--2019"
__maintainer__ = "Henrik Finsberg"
__email__ = "henriknf@simula.no"
__license__ = """
c) 2001-2020 Simula Research Laboratory ALL RIGHTS RESERVED

END-USER LICENSE AGREEMENT
PLEASE READ THIS DOCUMENT CAREFULLY. By installing or using this
software you agree with the terms and conditions of this license
agreement. If you do not accept the terms of this license agreement
you may not install or use this software.

Permission to use, copy, modify and distribute any part of this
software for non-profit educational and research purposes, without
fee, and without a written agreement is hereby granted, provided
that the above copyright notice, and this license agreement in its
entirety appear in all copies. Those desiring to use this software
for commercial purposes should contact Simula Research Laboratory AS:
post@simula.no

IN NO EVENT SHALL SIMULA RESEARCH LABORATORY BE LIABLE TO ANY PARTY
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE
"MPS" EVEN IF SIMULA RESEARCH LABORATORY HAS BEEN ADVISED
OF THE POSSIBILITY OF SUCH DAMAGE. THE SOFTWARE PROVIDED HEREIN IS
ON AN "AS IS" BASIS, AND SIMULA RESEARCH LABORATORY HAS NO OBLIGATION
TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
SIMULA RESEARCH LABORATORY MAKES NO REPRESENTATIONS AND EXTENDS NO
WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESSED, INCLUDING, BUT
NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS
"""
import concurrent.futures
import itertools as it
import json
import logging
from collections import namedtuple
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.ndimage import gaussian_filter1d

from . import average, bin_utils, plotter, utils

logger = utils.get_logger(__name__)

chopped_data = namedtuple("chopped_data", "data, times, pacing, parameters")
chopping_parameters = namedtuple(
    "chopping_parameters", "use_pacing_info, extend_end, extend_front"
)
mps_prevalence = namedtuple(
    "mps_prevalence", "prevalence, tissue_covered_area, is_beating, is_tissue"
)
Upstroke = namedtuple(
    "Upstroke",
    ["index", "value", "upstroke", "dt", "x0", "sigmoid", "time_APD20_to_APD80"],
)
APDAnalysis = namedtuple(
    "APDAnalysis",
    [
        "apds",
        "capds",
        "apd_points",
        "is_significant",
        "triangulation",
        "slope_APD80",
        "slope_cAPD80",
        "const_APD80",
        "const_cAPD80",
        "apd_dt",
    ],
)


def snr(y):
    return np.mean(y) / np.std(y)


def active_pixel_mask(frames, cutoff_factor=0.5, *args, **kwargs):

    logger.debug("Get active pixel mask")

    avg, inds = average.get_temporal_average(
        frames, alpha=cutoff_factor, return_indices=True
    )

    mask = np.ones_like(frames.T[0], dtype=bool).reshape(-1)
    mask[inds] = False
    mask = mask.reshape(frames.T[0].T.shape)
    return mask


def remove_points(x, y, t_start, t_end, normalize=True):
    """
    Remove points in x and y between start and end.
    Also make sure that the new x starts a zero if
    normalize = True
    """

    print(max(x), min(x), t_start, t_end)
    start = next(i for i, t in enumerate(x) if t > t_start) - 1
    try:
        end = next(i for i, t in enumerate(x) if t > t_end)
    except StopIteration:
        end = len(x) - 1

    logger.debug(
        ("Remove points for t={} (index:{}) to t={} (index:{})" "").format(
            t_start, start, t_end, end
        )
    )
    x0 = x[:start]
    x1 = np.subtract(x[end:], x[end] - x[start])
    x_new = np.concatenate((x0, x1))

    if normalize:
        x_new -= x_new[0]
    y_new = np.concatenate((y[:start], y[end:]))

    return x_new, y_new


def remove_spikes(y, pacing, spike_duration=7):
    """
    Remove spikes from signal

    Arguments
    ---------
    y : array
        The signal where you want to remove spikes
    pacing : array
        The pacing amplitude of same length as y
    """

    msg = (
        "Pacing and signal must be of equal length, "
        "got len(y) = {}, len(pacing) = {}"
    ).format(len(y), len(pacing))
    assert len(y) == len(pacing), msg

    # Find time of pacing
    (inds,) = np.where(np.diff(np.array(pacing, dtype=float)) > 0)

    if len(inds) == 0:
        logger.warning("No pacing found. Spike removal not possible.")
        return y

    spike_points = np.concatenate(
        [np.arange(i, i + spike_duration) for i in inds]
    ).astype(int)

    return np.delete(y, spike_points)


def beating_frequency(times, unit="ms"):
    """
    Returns the approximate beating frequency
    in Hz
    """
    if len(times) == 0:
        return np.nan
    # Get chopped data
    # Find the average lenght of each beat in time
    t_mean = np.mean([ti[-1] - ti[0] for ti in times])
    # Convert to seconds
    if unit == "ms":
        t_mean /= 1000.0
    # Return the freqency
    return 1.0 / t_mean


def beating_frequency_modified(data, times, unit="ms"):
    """
    Returns the approximate beating frequency
    in Hz
    """

    t_maxs = [t[np.argmax(c)] for c, t in zip(data, times)]

    dt = np.diff(t_maxs)
    if unit == "ms":
        dt = np.divide(dt, 1000.0)
    return np.divide(1, dt)


def find_pacing_period(pacing_data, pacing_amp=5):
    """
    Find period of pacing in pacing data

    Arguments
    ---------
    pacing_data : array
        The pacing data
    pacing_amp : float
        The stimulus amplitude

    """

    periods = []
    pace = np.where(pacing_data == pacing_amp)[0]
    for i, j in zip(pace[:-1], pace[1:]):
        if not j == i + 1:
            periods.append(j - i)
    if len(periods) > 0:
        return int(np.median(periods))
    else:
        return None


def apd(V, percent, t=None, v_r=None, return_coords=False, rule=0, use_spline=True):
    r"""
    Return the action potential duration at the given
    percent repolarization, so that percent = 0
    would be zero, and percent = 1 given the time from triggering
    to potential is down to resting potential.

    Arguments
    ---------

    V : array
        The signal
    percent: float
        The percentage (value between zero and 1)
    t : array
        List with times
    v_r : float
        The resting potential(optional). If not provided the minimum value
        will be used
    return_coords:
        If True, return the coordinates of the start and stop
        of the APD, and not the duration itself
    rule : int
        If 0 APD_X is computed as the difference between
        the intersection points of the line at X percent from the max
        If 1 APD_X is computed from the maximum derivative to the
        intersection line. (Default: 0)
    use_spline : bool
        Use spline interpolation to get exact localization of the zeros
        (Default : True)

    Returns
    -------
    duration : float
        The percentage action potential duration

    .. Note::

        If the signal has more intersection than two with the
        APX_X line, then the first and last occurence will be used.


    Notes
    -----
    If the signal represent voltage, we would like to compute the action
    potential duration at a given percentage of the signals. Formally,
    Let :math:`p \in (0, 100)`, and define

    .. math::
        y_p = \tilde{y} - \left(1- \frac{p}{100} \right).

    Now let :math:`\mathcal{T}` denote the set of solutions of :math:`y_p = 0`,
    i.e

    .. math::
        y_p(t) = 0, \;\; \forall t \in \mathcal{T}

    Then there are three different scenarios; :math:`\mathcal{T}` contains one,
    two or more than two elements. Note that :math:`\mathcal{T}` cannot be
    empty (because of the intermediate value theorem).
    The only valid scenario is the case when :math:`\mathcal{T}` contains
    two (or more) elements. In this case we define

    .. math::
        \mathrm{APD} \; p = \max \mathcal{T} - \min \mathcal{T}

    for :math:`p < 0.5` and

    .. math::
        \mathrm{APD} \; p = \min \mathcal{T} / \min \mathcal{T}  - \min \mathcal{T}

    """
    if t is None:
        t = range(len(V))

    msg = "The signal and time are not of same lenght"
    assert len(t) == len(V), msg

    y = normalize_signal(V, v_r) - (1 - percent)

    if use_spline:

        k = 3 if rule == 0 else 5
        try:
            f = UnivariateSpline(t, y, s=0, k=k)
        except Exception as ex:
            logger.warning(
                (
                    "Unable to compute APD {}. " "Please change your settings, {}" ""
                ).format(percent * 100, ex)
            )
            return None

        inds = f.roots()
    else:
        inds = t[np.where(np.diff(np.sign(y)))[0]]

    if len(inds) == 0:
        logger.warning("Warning: no root was found for APD {}".format(percent))
        x1 = x2 = 0
    if len(inds) == 1:
        x1 = x2 = inds[0]
        logger.warning("Warning: only one root was found for APD {}" "".format(percent))
    else:
        start_index = np.argmax(np.diff(inds))
        x1 = sorted(inds)[start_index]
        x2 = sorted(inds)[start_index + 1]
        # if percent < 0.5:
        #     x2 = sorted(inds)[-1]
        # else:
        #     x2 = sorted(inds)[1]

    if rule == 1:

        if use_spline:
            h = f.derivative()
            x1 = np.interp(np.argmax(h(h.derivative().roots())), range(len(t)), t)
        else:
            x1 = np.argmax(np.diff(y))

    if return_coords:
        g = UnivariateSpline(t, V, s=0)
        val_x1 = g(x1)
        val_x2 = g(x2)

        val_th = np.min(V) + (1 - percent) * (np.max(V) - np.min(V))

        return np.array([x1, x2]), np.array([val_x1, val_x2, val_th])

    else:
        return x2 - x1


def tau(x, y, a=0.75):
    """
    Decay time. Time for the signal amplitude to go from maxium to
    (1 - a) * 100 % of maximum

    Arguments
    ---------
    x : array
        The time stamps
    y : array
        The signal
    a : The value for which you want to estimate the time decay
    """
    Y = UnivariateSpline(x, normalize_signal(y) - a, s=0, k=3)
    t_max = x[np.argmax(y)]
    r = Y.roots()
    if len(r) >= 2:
        t_a = r[1]
    elif len(r) == 1:
        logger.warning(
            (
                "Only one zero was found when computing tau{}. " "Result might be wrong"
            ).format(int(a * 100))
        )
        t_a = r[0]
    else:
        logger.warning(
            (
                "No zero found when computing tau{}. "
                "Return the value of time to peak"
            ).format(int(a * 100))
        )
        t_a = x[0]

    return t_a - t_max


def time_to_peak(x, y, pacing=None):
    """
    Computed the time to peak from pacing is
    triggered to maximum amplitude. Note, if pacing info
    is not provided it will compute the time from
    the beginning of the trace (which might not be consistent)
    to the peak.

    Arguments
    ---------
    x : array
        The time stamps
    y : array
        The signal
    pacing : array
        The pacing amplitude
    """

    if pacing is None:
        return x[np.argmax(y)]

    t_max = x[np.argmax(y)]
    if pacing is None:
        t_start = x[0]
    else:
        try:
            start_idx = (
                next(i for i, p in enumerate(np.diff(pacing.astype(float))) if p > 0)
                + 1
            )
        except StopIteration:
            start_idx = 0
    t_start = x[start_idx]

    return t_max - t_start


def upstroke(x, y, a=0.8):
    """
    Compute the time from (1-a)*100 % signal
    amplitude to peak. For example if if a = 0.8
    if will compute the time from the starting value
    of APD80 to the upstroke.
    """

    Y = UnivariateSpline(x, normalize_signal(y) - (1 - a), s=0, k=3)
    t_max = x[np.argmax(y)]
    r = Y.roots()
    if len(r) >= 1:
        if len(r) == 1:
            logger.warning(
                (
                    "Only one zero was found when computing upstroke{}. "
                    "Result might be wrong"
                ).format(int(a * 100))
            )
        t_a = r[0]
    else:
        logger.warning(
            (
                "No zero found when computing upstroke{}. "
                "Return the value of time to peak"
            ).format(int(a * 100))
        )
        t_a = x[0]

    return t_max - t_a


def interpolate(t: np.ndarray, trace: np.ndarray, dt: float = 1.0):

    f = UnivariateSpline(t, trace, s=0, k=1)
    t0 = np.arange(t[0], t[-1] + dt, dt)
    return t0, f(t0)


def find_upstroke_values(
    t: np.ndarray, y: np.ndarray, upstroke_duration: int = 50, normalize: bool = True
) -> np.ndarray:

    # Find intersection with APD50 line
    y_mid = (np.max(y) + np.min(y)) / 2
    f = UnivariateSpline(t, y - y_mid, s=0, k=3)
    zeros = f.roots()
    if len(zeros) == 0:
        return np.array([])
    idx_mid = next(i for i, ti in enumerate(t) if ti > zeros[0])

    # Upstroke should not be more than 50 ms
    N = upstroke_duration // 2
    upstroke = y[idx_mid - N : idx_mid + N + 1]

    if normalize and len(upstroke) > 0:
        upstroke_smooth = gaussian_filter1d(upstroke, 2.0)

        y_min = np.min(upstroke_smooth)
        y_max = np.max(upstroke_smooth)

        upstroke_normalized = (upstroke - y_min) / (y_max - y_min)

        return upstroke_normalized

    return upstroke


def time_between_APDs(t, y, from_APD, to_APD):
    """Find the duration between first intersection
    of two APD lines

    Arguments
    ---------
    t : np.ndarray
        Time values
    y : np.ndarray
        The trace
    from_APD: int
        First APD line
    to_APD: int
        Second APD line

    Returns
    -------
    float:
        The time between `from_APD` to `to_APD`

    Example
    -------

    .. code:: python

        # Compute the time between APD20 and APD80
        t2080 = time_between_APDs(t, y, 20, 80)

    """
    y_norm = normalize_signal(y)

    if not (0 < from_APD < 1):
        from_APD /= 100
    y_from = UnivariateSpline(t, y_norm - from_APD, s=0, k=3)
    t_from = y_from.roots()[0]

    if not (0 < to_APD < 1):
        to_APD /= 100
    y_to = UnivariateSpline(t, y_norm - to_APD, s=0, k=3)
    t_to = y_to.roots()[0]

    return t_to - t_from


def max_relative_upstroke_velocity(
    t: np.ndarray, y: np.ndarray, upstroke_duration: int = 50, sigmoid_fit: bool = True
) -> Upstroke:
    """Estimate maximum relative upstroke velocity

    Arguments
    ---------
    t : np.ndarray
        Time values
    y : np.ndarray
        The trace
    upstroke_duration : int
        Duration in milliseconds of upstroke (Default: 50).
        This does not have to be exact up should at least be
        longer than the upstroke.
    sigmoid_fit : bool
        If True then use a sigmoid function to fit the data
        of the upstroke and report the maximum derivate of
        the sigmoid as the maximum upstroke.

    Notes
    -----
    Brief outline of current algorithm:
    1. Interpolate data to have time resolution of 1 ms.
    2. Find first intersection with ADP50 line
    3. Select 25 ms before and after the point found in 2
    4. Normalize the 50 ms (upstroke_duration) from 3 so that we have
    a max and min value of 1 and 0 respectively.
    If sigmoid fit is True
    Fit the normalize trace to a sigmoid function and compte the

    else:
    5. Compute the successive differences in amplitude (i.e delta y)
    and report the maximum value of these


    """

    # Interpolate to 1ms precision
    t0, y0 = interpolate(t, y, dt=1.0)
    # Find values beloning to upstroke
    upstroke = find_upstroke_values(t0, y0, upstroke_duration=upstroke_duration)
    dt = np.mean(np.diff(t))
    if len(upstroke) == 0:
        # There is no signal
        index = 0
        value = np.nan
        x0 = None
        s = None
        time_APD20_to_APD80 = np.nan
    else:
        t_upstroke = t0[: len(upstroke)]
        t_upstroke -= t_upstroke[0]
        if sigmoid_fit:

            def sigmoid(x, k, x0):
                return 0.5 * (np.tanh(k * (x - x0)) + 1)

            from scipy.optimize import curve_fit

            popt, pcov = curve_fit(sigmoid, t_upstroke, upstroke, method="dogbox")

            k, x0 = popt
            index = None  # type:ignore
            s = sigmoid(t_upstroke, k, x0)
            value = k / 2
            time_APD20_to_APD80 = time_between_APDs(t_upstroke, s, 20, 80)

        else:
            # Find max upstroke
            index = np.argmax(np.diff(upstroke))  # type:ignore
            value = np.max(np.diff(upstroke))
            x0 = None
            s = None
            time_APD20_to_APD80 = time_between_APDs(t_upstroke, upstroke, 20, 80)

    return Upstroke(
        index=index,
        value=value,
        upstroke=upstroke,
        dt=dt,
        x0=x0,
        sigmoid=s,
        time_APD20_to_APD80=time_APD20_to_APD80,
    )


def maximum_upstroke_velocity(y, t=None, use_spline=True, normalize=False):
    r"""
    Compute maximum upstroke velocity

    Arguments
    ---------
    y : array
        The signal
    t : array
        The time points
    use_spline : bool
        Use spline interpolation
        (Default : True)
    normalize : bool
        If true normalize signal first, so that max value is 1.0,
        and min value is zero before performing the computation.

    Returns
    -------
    float
        The maximum upstroke velocity

    Notes
    -----
    If :math:`y` is the voltage, this feature corresponds to the
    maximum upstroke velocity. In the case when :math:`y` is a continuous
    signal, i.e a 5th order spline interpolant of the data, then we can
    simply find the roots of the second derivative, and evaluate the
    derivative at that point, i.e

    .. math::
        \max \frac{\mathrm{d}y}{\mathrm{d}t}
        = \frac{\mathrm{d}y}{\mathrm{d}t} (t^*),

    where

    .. math::
        \frac{\mathrm{d}^2y}{\mathrm{d}^2t}(t^*) = 0.

    If we decide to use the discrete version of the signal then the
    above each derivative is approximated using a forward difference, i.e

    .. math::
        \frac{\mathrm{d}y}{\mathrm{d}t}(t_i) \approx
        \frac{y_{i+1}-y_i}{t_{i+1} - t_i}.

    """

    if t is None:
        t = range(len(y))

    msg = "The signal and time are not of same lenght"
    assert len(t) == len(y), msg

    # Normalize
    if normalize:
        y = normalize_signal(y)

    if use_spline:
        f = UnivariateSpline(t, y, s=0, k=5)
        h = f.derivative()
        max_upstroke_vel = np.max(h(h.derivative().roots()))
    else:
        max_upstroke_vel = np.max(np.divide(np.diff(y), np.diff(t)))

    return max_upstroke_vel


def integrate_apd(y, t=None, percent=0.3, use_spline=True, normalize=False):
    r"""
    Compute the integral of the signals above
    the APD p line

    Arguments
    ---------
    y : array
        The signal
    t : array
        The time points
    use_spline : bool
        Use spline interpolation
        (Default : True)
    normalize : bool
        If true normalize signal first, so that max value is 1.0,
        and min value is zero before performing the computation.


    Returns
    -------
    integral : float
        The integral above the line defined by the APD p line.

    Notes
    -----
    This feature represents the integral of the signal above the
    horizontal line defined by :math:`\mathrm{APD} p`. First let

    .. math::
        y_p = \tilde{y} - \left(1- \frac{p}{100} \right),

    and let :math:`t_1` and :math:`t_2` be the two solutions solving
    :math:`y_p(t_i) = 0, i=1,2` (assume that we have 2 solutions only).
    The integral we are seeking can now be computed as follows:

    .. math::
        \mathrm{Int} \; p = \int_{t_1}^{t_2}
        \left[ y - y(t_1) \right] \mathrm{d}t

    """

    if normalize:
        y = normalize_signal(y)

    if t is None:
        t = range(len(y))

    coords, vals = apd(y, percent, t, return_coords=True, use_spline=use_spline)

    if use_spline:
        Y = y - vals[0]
        f = UnivariateSpline(t, Y, s=0, k=3)
        integral = f.integral(*coords)

    else:
        Y = y - vals[-1]

        t1 = t.tolist().index(coords[0])
        t2 = t.tolist().index(coords[1]) + 1

        integral = np.sum(np.multiply(Y[t1:t2], np.diff(t)[t1:t2]))

    return integral


def normalize_signal(V, v_r=None):
    """
    Normalize signal to have maximum value 1
    and zero being the value equal to v_r (resting value).
    If v_r is not provided the minimum value
    in V will be used as v_r

    Arguments
    ---------
    V : array
        The signal
    v_r : float
        The resting value

    """

    # Maximum valu
    v_max = np.max(V)

    # Baseline or resting value
    if v_r is None:
        v_r = np.min(V)

    return (np.array(V) - v_r) / (v_max - v_r)


def time_unit(time_stamps):

    dt = np.mean(np.diff(time_stamps))
    # Assume dt is larger than 0.5 ms and smallar than 0.5 seconds
    unit = "ms" if dt > 0.5 else "s"
    return unit


def average_intensity(data, **kwargs):
    """
    Compute the average_intensity of the frame stack.
    The available keyword arguments depends on the averaging_type
    chosen

    Arguments
    ---------
    X : :class:`numpy.ndarray`
        The frame stack (size :math:`M \times N \times T`)
    averaging_type : str
        How you want to average the frame stack. Possible values
        are ['all', 'temporal', 'spatial', 'global']

    """
    logger.debug("Compute average intensity")
    from .load import MPS

    if isinstance(data, MPS):
        # Get the frames
        data = data.frames

    mask = kwargs.pop("mask", None)
    alpha = kwargs.pop("alpha", 1.0)
    averaging_type = kwargs.pop("averaging_type", "temporal")

    if alpha == 1.0:
        # Mean of everything
        if mask is None:
            avg = average.get_average_all(data)
        else:
            avg = average.masked_average(data, mask)

    else:

        if averaging_type == "spatial":
            avg = average.get_spatial_average(data, **kwargs)

        elif averaging_type == "temporal":
            avg = average.get_temporal_average(data, **kwargs)
        else:
            msg = (
                "Unknown averaging_type {}. Expected averaging type to "
                'be one of ["all", "spatial", "temporal"]'
            )
            logger.error(msg)
            raise ValueError(msg)

    return avg


def filt(y, kernel_size=None):
    """
    Filer signal using a median filter.
    Default kernel_size is 3
    """
    if kernel_size is None:
        kernel_size = 3

    logger.debug("\nFilter image")
    from scipy.signal import medfilt

    smooth_trace = medfilt(y, kernel_size)
    return smooth_trace


def correct_background(x, y, *args, **kwargs):
    r"""
    Perform at background correction.
    First estimate background :math:`b`, and let
    :math:`F_0 = b(0)`. The corrected background is
    then :math:`\frac{y - b}{F_0}`

    Arguments
    ---------
    x : np.mdarray
        Time points
    y : Fluorecense amplitude

    Returns
    -------
    np.mdarrray
        The corrected trace
    """

    bkg = background(x, y, *args, **kwargs)
    F0 = bkg[0]
    corrected = (1 / F0) * (y - bkg)
    return corrected


def corrected_apd(apd, beat_rate, formula="friderica"):
    """Correct the given APD (or any QT measurement) for the beat rate.
    normally the faster the HR (or the shorter the RR interval),
    the shorter the QT interval, and vice versa

    Friderica formula (default):

    .. math::

        APD (RR)^{-1/3}

    Bazett formula:

    .. math::

        APD (RR)^{-1/2}

    """

    formulas = ["friderica", "bazett"]
    msg = f"Expected formula to be one of {formulas}, got {formula}"
    assert formula in formulas, msg
    RR = np.divide(60, beat_rate)

    if formula == "friderica":
        return np.multiply(apd, pow(RR, -1 / 3))
    else:
        return np.multiply(apd, pow(RR, -1 / 2))


def analyze_apds(
    chopped_data: List[List[float]],
    chopped_times: List[List[float]],
    max_allowed_apd_change: Optional[float] = None,
    fname: str = "",
    plot=True,
) -> APDAnalysis:

    apd_levels = np.sort(np.append(np.arange(0.1, 0.91, 0.2), 0.8))
    apds = {
        int(100 * k): [apd(y, k, t) for y, t in zip(chopped_data, chopped_times)]
        for k in apd_levels
    }

    apd_points = {
        int(100 * k): [
            apd(y, k, t, return_coords=True)
            for y, t in zip(chopped_data, chopped_times)
        ]
        for k in apd_levels
    }

    median_freq = np.median(beating_frequency_modified(chopped_data, chopped_times))
    beat_rates = {}
    apd_dt = {}
    for k, apdx in apd_points.items():
        apd_dt[k] = [v[0][0] for v in apdx]
        if len(apd_dt[k]) > 1:
            # Time between APD80
            # Let last beat have the median beat rate
            beat_rates[k] = [
                60 * (v1 - v0) / 1000.0 for v0, v1 in zip(apd_dt[k][:-1], apd_dt[k][1:])
            ] + [60.0 / median_freq]

    if len(apds[80]) > 1:
        a, const_APD80 = np.polyfit(apd_dt[80], apds[80], deg=1)
        slope_APD80 = a * 1000 * 60  # d APD / min
    else:
        slope_APD80 = np.nan
        const_APD80 = np.nan

    capds = {
        k: corrected_apd(v, b).tolist()
        for (k, v), b in zip(apds.items(), beat_rates.values())
    }

    if len(apds[80]) > 1:
        a, const_cAPD80 = np.polyfit(apd_dt[80], capds[80], deg=1)
        slope_cAPD80 = a * 1000 * 60  # d APD / min
    else:
        slope_cAPD80 = np.nan
        const_cAPD80 = np.nan

    triangulation = []
    for y, t in zip(chopped_data, chopped_times):
        apd30 = apd(y, 0.3, t, return_coords=True)
        apd80 = apd(y, 0.8, t, return_coords=True)

        tri = apd80[0][-1] - apd30[0][-1]
        if tri < 0:
            # Something is wrong
            tri = np.nan
        triangulation.append(tri)

    median_apds = {k: np.median(v) for k, v in apds.items()}

    if max_allowed_apd_change is not None:
        max_diff = {k: float(max_allowed_apd_change) for k in apds.keys()}
    else:
        max_diff = {k: float(np.std(v)) for k, v in apds.items()}

    is_significant = {
        k: np.abs(np.array(v) - median_apds[k]) > max_diff[k] for k, v in apds.items()
    }

    msg = "Found the following number of significant beats based on APDs: \n"
    msg += "\n".join([f"APD{k}: {sum(v)}" for k, v in is_significant.items()])
    logger.info(msg)

    res = APDAnalysis(
        apds=apds,
        capds=capds,
        apd_points=apd_points,
        is_significant=is_significant,
        triangulation=triangulation,
        slope_APD80=slope_APD80,
        slope_cAPD80=slope_cAPD80,
        const_APD80=const_APD80,
        const_cAPD80=const_cAPD80,
        apd_dt=apd_dt,
    )

    if plot:
        plot_apd_analysis(
            chopped_data=chopped_data, chopped_times=chopped_times, res=res, fname=fname
        )

    return res


def plot_apd_analysis(chopped_data, chopped_times, res: APDAnalysis, fname=""):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return

    width = 1 / (len(res.apds) + 1)
    n = len(next(iter(res.apds.values())))

    fig, ax = plt.subplots(3, 1, figsize=(15, 10))
    for i, (label, y) in enumerate(res.apds.items()):
        x = np.arange(n) + i * width - 0.5
        ax[0].bar(x, y, width=width, label=f"{label:.0f}")
        ax[1].plot(np.arange(n) - 0.5, y, marker="o", label=f"{label:.0f}")
    ax[1].plot(
        np.arange(n) - 0.5,
        res.const_APD80 + (res.slope_APD80 / (60 * 1000)) * np.array(res.apd_dt[80]),
        "k--",
        label=f"{res.slope_APD80:.2f} x + C (APD80)",
    )
    ax[1].plot(
        np.arange(n) - 0.5,
        res.const_cAPD80 + (res.slope_cAPD80 / (60 * 1000)) * np.array(res.apd_dt[80]),
        "k:",
        label=f"{res.slope_cAPD80:.2f} x + C (cAPD80)",
    )
    # Mark significatn patches
    for sig, patch in zip(
        np.array(list(res.is_significant.values())).flatten(), ax[0].patches
    ):
        if not sig:
            continue
        ax[0].text(
            patch.get_x() + patch.get_width() / 2.0,
            patch.get_height(),
            "*",
            ha="center",
        )

    for i, axi in enumerate(ax[:2]):
        if i == 1:
            axi.set_xlabel("Beat numbers")
        axi.set_ylabel("APD [ms]")
        axi.set_xticks(np.arange(n) - 0.5)
        axi.set_xticklabels(np.arange(n))
        axi.legend()
        axi.grid()

    for c, t in zip(chopped_data, chopped_times):
        ax[2].plot(t, c, color="k")
    for label, vals in res.apd_points.items():
        x = [v[0][0] for v in vals] + [v[0][1] for v in vals]
        y = [v[1][0] for v in vals] * 2
        ax[2].plot(x, y, linestyle="", marker="o", label=label)

    ax[2].legend()
    ax[2].grid()
    ax[2].set_xlabel("Time [ms]")
    ax[2].set_ylabel(r"$\Delta F / F$")

    if fname != "":
        fig.savefig(fname)
    plt.close()


def background(x, y, order=2, s=0.01, fct="atq"):
    """
    Compute an estimation of the background (aka baseline)
    in chemical spectra


    This is a reimplementation of a MATLAB script that
    can be found `here <https://se.mathworks.com/matlabcentral
    /fileexchange/27429-background-correction>`

    .. rubric:: References

    [1] Mazet, V., Carteret, C., Brie, D., Idier, J. and Humbert, B.,
    2005. Background removal from spectra by designing and minimising
    a non-quadratic cost function. Chemometrics and intelligent
    laboratory systems, 76(2), pp.121-133.

    """

    # Rescaling
    N = len(x)
    x, i = np.sort(x), np.argsort(x)

    y = y[i]

    maxy = np.max(y)
    dely = (maxy - np.min(y)) / 2
    # Normalize time
    x = 2 * np.divide(np.subtract(x, x[N - 1]), x[N - 1] - x[0]) + 1
    y = np.subtract(y, maxy) / dely + 1

    # Vandermonde matrix
    T = np.vander(x, order + 1)
    Tinv = np.linalg.pinv(T.T.dot(T)).dot(T.T)

    # Initialisation (least-squares estimation)
    a = Tinv.dot(y)
    z = T.dot(a)

    #  Other variables
    alpha = 0.99 * 1 / 2  # Scale parameter alpha
    it = 0  # Iteration number
    zp = np.ones(N)  # Previous estimation

    # Iterate
    while np.sum((z - zp) ** 2) / np.sum(zp ** 2) > 1e-9:

        it = it + 1  # Iteration number
        zp = z  # Previous estimation
        res = y - z  # Residual

        d = np.zeros(len(res))

        # Estimate d
        if fct == "sh":
            d[np.abs(res) < s] += res[np.abs(res) < s] * (2 * alpha - 1)
            d[res <= -s] -= alpha * 2 * s + res[res <= -s]
            d[res >= s] -= res[res >= s] - alpha * 2 * s

        elif fct == "ah":
            d[np.abs(res) < s] += res[np.abs(res) < s] * (2 * alpha - 1)
            d[res >= s] -= res[res >= s] - alpha * 2 * s
        elif fct == "stq":
            d[np.abs(res) < s] += res[np.abs(res) < s] * (2 * alpha - 1)
            d[res >= s] -= res[res >= s] - alpha * 2 * s

        elif fct == "atq":
            d[res < s] += res[res < s] * (2 * alpha - 1)
            d[res >= s] -= res[res >= s]

        # Estimate z
        a = Tinv.dot(y + d)  # Polynomial coefficients a
        z = T.dot(a)  # Polynomial

    # Rescaling
    j = np.argsort(i)
    z = (z[j] - 1) * dely + maxy

    a[0] = a[0] - 1
    a = a * dely  # + maxy
    return z


def chop_data(data, time, **kwargs):

    if time is None:
        time = np.arange(len(data))

    pacing = kwargs.pop("pacing", np.zeros(len(time)))

    if all(pacing == 0):
        logger.debug("Chop data without pacing")
        return chop_data_without_pacing(data, time, **kwargs)
    else:
        logger.debug("Chop data with pacing")
        return chop_data_with_pacing(data, time, pacing, **kwargs)


class EmptyChoppingError(ValueError):
    pass


def chop_data_without_pacing(
    data,
    time,
    threshold_factor=None,
    extend_front=None,
    extend_end=None,
    min_window=None,
    max_window=None,
    winlen=None,
    N=None,
    **kwargs,
):

    r"""
    Chop data into beats

    Arguments
    ---------
    data : list or array
        The data that you want to chop
    time : list or array
        The time points. If none is provided time will be
        just the indices.
    threshold_factor : float
        Thresholds for where the signal should be chopped (Default: 0.3)
    extend_front : scalar
        Extend the start of each subsignal this many milliseconds
        before the threshold is detected. Default: 300 ms
    extend_end : scalar
        Extend the end of each subsignal this many milliseconds.
        Default 60 ms.
    min_window : scalar
        Length of minimum window
    max_window : scalar
        Length of maximum window
    N : int
        Length of output signals

    Returns
    -------
    chopped_data : list
        List of chopped data
    chopped_times : list
        List of chopped times
    chopped_pacing : list
        List of chopped pacing amps (which are all zero)



    Notes
    -----

    The signals extracted from the MPS data consist of data from several beats,
    and in order to e.g compute properties from an average beat, we need a way
    to chop the signal into different beats. Suppose we have the
    signal :math:`y(t),
    t \in [0, T]`, where we assume that filtering and background correction
    have allready been applied. Suppose we have :math:`N` sub-signals,
    :math:`z_1, z_2, \cdots, z_N`, each representing one beat. Let
    :math:`\tau_i = [\tau_i^0, \tau_i^1], i = 1, \cdots N` be
    non-empty intervals corresponding to the support of each sub-signal
    ( :math:`\tau_i = \{ t \in [0,T]: z_i(t) \neq 0 \}`), with
    :math:`\tau_i^j < \tau_{i+1}^j \forall i` and :math:`\bigcup_{i=1}^N
    \tau_i = [0, T]`. Note that the intersection :math:` \tau_i \cap
    \tau_{i+1}` can be non-empty. The aim is now to find good candidates
    for :math:`\tau_i^j , i = 1 \cdots N, j = 0,1`. We have two different
    scenarios, namely with or without pacing information.

    If pacing information is available we can e.g set :math:`\tau_i^0`
    to be 30 ms before each stimulus is applied and :math:`\tau_i^1`
    to be 60 ms before the next stimulus is applied.

    If pacing information is not available then we need to estimate the
    beginning and end of each interval. We proceed as follows:


    1.  Choose some threshold value :math:`\eta` (default :math:`\eta = 0.5`),
        and compute :math:`y_{\eta}`

    2. Let :math:`h = y - y_{\eta}`,  and define the set

    .. math::
        \mathcal{T}_0 = \{ t : h(t) = 0 \}

    3. Sort the elements of :math:`\mathcal{T}_0` in increasing order,
       i.e :math:`\mathcal{T}_0 = (t_1, t_2, \cdots, t_M)`,
       with :math:`t_i < t_{i+1}`.

    4. Select a minimum duration of a beat :math:`\nu` (default
       :math:`\nu = 50` ms) and define the set

    .. math::
        \mathcal{T}_1 = \{t_i \in \mathcal{T}_0 : t_{i+1} - t_i > \eta \}

    to be be the subset of :math:`\mathcal{T}_0` where the difference between
    to subsequent time stamps are greater than :math:`\eta`.

    5. Let :math:`\delta > 0` (default :math:`\delta` = 0.1 ms), and set

    .. math::
        \{\tilde{\tau_1}^0, \tilde{\tau_2}^0, \cdots,\tilde{\tau_N}^0\}
        := \{ t \in  \mathcal{T}_1 : h(t + \delta) > 0 \} \\
        \{\tilde{\tau_1}^1, \tilde{\tau_2}^1, \cdots,\tilde{\tau_N}^1\}
        := \{ t \in  \mathcal{T}_1 : h(t + \delta) < 0 \} \\

    Note that we need to also make sure that all beats have a start and an end.

    6. Extend each subsignal at the beginning and end, i,e
       if :math:`a,b \geq 0`, define

    .. math::
        \tau_i^0 = \tilde{\tau_i}^0 - a \\
        \tau_i^1 = \tilde{\tau_i}^1 + b
    """
    logger.debug("Chopping without pacing")
    threshold_factor = threshold_factor if threshold_factor is not None else 0.3
    extend_end = extend_end if extend_end is not None else 300
    extend_front = extend_front if extend_front is not None else 100
    min_window = min_window if min_window is not None else 50
    max_window = max_window if max_window is not None else 2000
    winlen = winlen if winlen is not None else 50

    chop_pars = chopping_parameters(
        use_pacing_info=False, extend_end=extend_end, extend_front=extend_front
    )
    logger.debug(f"Use chopping parameters: {chop_pars}")
    empty = chopped_data(data=[], times=[], pacing=[], parameters=chop_pars)

    # Make a spline interpolation
    data_spline = UnivariateSpline(time, data, s=0)

    try:
        starts, ends, zeros = locate_chop_points(
            time, data, threshold_factor, chop_pars
        )
    except EmptyChoppingError:
        return empty

    if len(zeros) <= 3:
        ## Just return the original data
        return chopped_data(
            data=[data],
            times=[time],
            pacing=[np.zeros_like(data)],
            parameters=chop_pars,
        )
    starts, ends = filter_start_ends_in_chopping(starts, ends, extend_front, extend_end)

    while len(ends) > 0 and ends[-1] > time[-1]:
        ends = ends[:-1]
        starts = starts[:-1]

    if len(ends) == 0:
        return empty

    # Update the ends to be the start of the next trace
    for i, s in enumerate(starts[1:]):
        ends[i] = s
    ends[-1] = min(ends[-1] + extend_end, time[-2])

    # Storage
    cutdata = []
    times = []

    for s, e in zip(starts, ends):

        if N is None:
            # Find the correct time points
            s_idx = next(i for i, si in enumerate(time) if si > s)
            e_idx = next(i for i, ei in enumerate(time) if ei > e)
            t = time[s_idx - 1 : e_idx + 1]

        else:
            t = np.linspace(s, e, N)

        if len(t) == 0:
            continue

        if t[-1] - t[0] < min_window:
            # Skip this one
            continue

        if t[-1] - t[0] > max_window:

            t_end = next(i for i, ti in enumerate(t) if ti - t[0] > max_window) + 1
            t = t[:t_end]

        sub = data_spline(t)

        cutdata.append(sub)
        times.append(t)

    pacing = [np.zeros(len(ti)) for ti in times]
    return chopped_data(data=cutdata, times=times, pacing=pacing, parameters=chop_pars)


def filter_start_ends_in_chopping(
    starts: Union[List[float], np.ndarray],
    ends: Union[List[float], np.ndarray],
    extend_front: float = 0,
    extend_end: float = 0,
):
    starts = np.array(starts)
    ends = np.array(ends)
    # If there is no starts return nothing
    if len(starts) == 0:
        raise EmptyChoppingError

    # The same with no ends
    if len(ends) == 0:
        raise EmptyChoppingError

    # If the first end is lower than the first start
    # we should drop the first end
    if ends[0] < starts[0]:
        ends = ends[1:]

    # And we should check this one more time
    if len(ends) == 0:
        raise EmptyChoppingError

    # Subtract the extend front
    starts = np.subtract(starts, extend_front)
    # If the trace starts in the middle of an event, that event is thrown out
    if starts[0] < 0:
        starts = starts[1:]

    new_starts = []
    new_ends = []
    for s in np.sort(starts):
        # if len(new_ends) !
        new_starts.append(s)
        for e in np.sort(ends):
            if e > s:
                new_ends.append(e)
                break
        else:
            # If no end was appended
            # we pop the last start
            new_starts.pop()

    return np.array(new_starts), np.array(new_ends)


def locate_chop_points(time, data, threshold_factor, chop_pars, winlen=50, eps=0.1):
    """FIXME"""
    # Some perturbation away from the zeros
    # eps = 0.1  # ms

    threshold = ((np.max(data) - np.min(data)) * threshold_factor) + np.min(data)

    # Data with zeros at the threshold
    data_spline_thresh = UnivariateSpline(time, data - threshold, s=0)
    # Localization of the zeros
    zeros_threshold_ = data_spline_thresh.roots()

    # Remove close indices
    inds = winlen < np.diff(zeros_threshold_)
    zeros_threshold = np.append(zeros_threshold_[0], zeros_threshold_[1:][inds])

    # Find the starts
    starts = zeros_threshold[data_spline_thresh(zeros_threshold + eps) > 0]

    if len(starts) == 0:
        logger.info("Chould not chop data. Try to reduce threshold factor")
        # Then something is wrong try agains with a lower threshold factor
        threshold_factor /= 2

        threshold = ((np.max(data) - np.min(data)) * threshold_factor) + np.min(data)

        # Data with zeros at the threshold
        data_spline_thresh = UnivariateSpline(time, data - threshold, s=0)
        # Localization of the zeros
        zeros_threshold_ = data_spline_thresh.roots()

        # Remove close indices
        inds = winlen < np.diff(zeros_threshold_)
        if len(zeros_threshold_) > 1:
            zeros_threshold = np.append(zeros_threshold_[0], zeros_threshold_[1:][inds])
            # Find the starts
            starts = zeros_threshold[data_spline_thresh(zeros_threshold + eps) > 0]

    # Find the endpoint where we hit the threshold
    ends = zeros_threshold[data_spline_thresh(zeros_threshold + eps) < 0]

    return starts, ends, zeros_threshold


def chop_data_with_pacing(
    data,
    time,
    pacing,
    extend_front=0,
    extend_end=0,
    min_window=300,
    max_window=2000,
    **kwargs,
):
    """
    Chop data based on pacing

    Arguments
    ---------
    data : array
        The data to be chopped
    time : array
        The time stamps for the data
    pacing : array
        The pacing amplitude
    extend_front : scalar
        Extend the start of each subsignal this many milliseconds
        before the threshold is detected. Default: 300 ms
    extend_end : scalar
        Extend the end of each subsignal this many milliseconds.
        Default 60 ms.
    min_window : int
        Minimum size of chopped signal

    Returns
    -------
    chopped_data : list
        List of chopped data
    chopped_times : list
        List of chopped times
    chopped_pacing : list
        List of chopped pacing amps


    Notes
    -----
    We first find where the pacing for each beat starts by finding
    the indices where the pacing amplitude goes form zero to
    something positive. Then we iterate over these indices and let
    the start of each beat be the start of the pacing (minus `extend_front`,
    which has a default value of zero) and the end will be the time of the
    beginning of the next beat (plus `extend_end` which is also has as
    default value of zero).
    """
    extend_end = extend_end if extend_end is not None else 0
    extend_front = extend_front if extend_front is not None else 0
    min_window = min_window if min_window is not None else 300
    max_window = max_window if max_window is not None else 2000

    # Find indices for start of each pacing
    start_pace_idx = np.where(np.diff(np.array(pacing, dtype=float)) > 0)[0].tolist()
    N = len(data)
    start_pace_idx.append(N)

    idx_freq = 1.0 / np.mean(np.diff(time))
    shift_start = int(idx_freq * extend_front)
    shift_end = int(idx_freq * extend_end)

    chopped_y = []
    chopped_pacing = []
    chopped_times = []

    for i in range(len(start_pace_idx) - 1):
        # Time
        end = min(start_pace_idx[i + 1] + shift_end, N)
        start = max(start_pace_idx[i] - shift_start, 0)

        t = time[start:end]
        if t[-1] - t[0] < min_window:
            continue

        if t[-1] - t[0] > max_window:
            t_end = next(i for i, ti in enumerate(t) if ti - t[0] > max_window) + 1
            end = start + t_end
            t = time[start:end]

        chopped_times.append(t)

        # Data
        c = data[start:end]
        chopped_y.append(c)

        # Pacing
        p = pacing[start:end]
        chopped_pacing.append(p)

    chop_pars = chopping_parameters(
        use_pacing_info=True, extend_end=extend_end, extend_front=extend_front
    )

    return chopped_data(
        data=chopped_y, times=chopped_times, pacing=chopped_pacing, parameters=chop_pars
    )


def compute_features(chopped_data, use_spline=True, normalize=False):
    r"""
    Analyze signals. Compute all features and
    include only the relevant beats


    Arguments
    ---------
    chopped_data : namedtuple
        The chopped data, time and pacing as a named tuple
    use_spline : bool
        Use spline interpolation
        (Default : True)
    normalize : bool
        If true normalize signal first, so that max value is 1.0,
        and min value is zero before performing the computation.


    Returns
    -------
    data : dict
        A dictionary with the following structure


    Notes
    -----
    In some cases, all the beats are not necessary representative
    for the underlying signal, for example if there is a lot of noise
    present. In this case it would be more robust compute the features by
    only including those beats that are closest to the average signals.
    Suppose we have :math:`N` sub-signals, :math:`z_1, z_2, \cdots, z_N`,
    each representing one beat. Further let :math:`f` denote a function
    which takes a sub-signal as input and output some of the properties,
    e.g :math:`f= \mathrm{APD}30`. Define

    .. math::
        \bar{f} = \frac{1}{N}\sum_{i = 1}^N f(z_i),

    to be the average value of a given feature :math:`f` of all
    sub-signals, and

    .. math::
        \sigma (f) = \sqrt{\frac{1}{N-1}\sum_{i = 1}^N
        \left(f(z_i) - \bar{f} \right)^2},

    be the standard deviation. Now, let :math:`\mathcal{D}`
    be the set of all sub-signals that are within 1 standard deviation
    of the mean, i.e

    .. math::
        \mathcal{D} = \{ z : | f(z) - \bar{f} | < \sigma(f) \}.

    Then a more robust estimate of the average value of :math:`f`
    (than :math:`\bar{f}`) is

    .. math::
        f_{\mathcal{D}} = \frac{1}{|\mathcal{D}|}\sum_{z \in
        \mathcal{D}}^N f(z).

    If the set :math:`\mathcal{D}` is empty, there can be two reasons
    for this, namely the signal is very noisy, or the standard deviation
    is very small. We assume the latter, and will in these cases only
    return :math:`\bar{f}`.

    """
    num_events = len(chopped_data.data)

    apd30 = np.zeros(num_events)
    apd50 = np.zeros(num_events)
    apd80 = np.zeros(num_events)
    apd90 = np.zeros(num_events)
    dFdt_max = np.zeros(num_events)
    int30 = np.zeros(num_events)
    tau75 = np.zeros(num_events)
    upstroke80 = np.zeros(num_events)

    for i, (data, time) in enumerate(zip(chopped_data.data, chopped_data.times)):
        unit = time_unit(time)
        unitfactor = 1000.0 if unit == "s" else 1.0

        time = np.multiply(unitfactor, time)

        apd90[i] = apd(data, 0.9, time, use_spline=use_spline) or np.nan
        apd80[i] = apd(data, 0.8, time, use_spline=use_spline) or np.nan
        apd50[i] = apd(data, 0.5, time, use_spline=use_spline) or np.nan
        apd30[i] = apd(data, 0.3, time, use_spline=use_spline) or np.nan
        tau75[i] = tau(time, data, 0.75)
        upstroke80 = upstroke(time, data, 0.8)

        dFdt_max[i] = (
            maximum_upstroke_velocity(
                data, time, use_spline=use_spline, normalize=normalize
            )
            or np.nan
        )
        int30[i] = (
            integrate_apd(data, time, 0.3, use_spline=use_spline, normalize=normalize)
            or np.nan
        )
    f = dict(
        apd30=apd30,
        apd50=apd50,
        apd80=apd80,
        apd90=apd90,
        dFdt_max=dFdt_max,
        int30=int30,
        tau75=tau75,
        upstroke80=upstroke80,
    )

    for k, v in f.items():
        f[k] = v[~np.isnan(v)]

    return f


def exclude_x_std(data, x=None, skips=None):
    """
    Given a list of values return a list of all
    values that are within x factor of the mean.

    Arguments
    ---------
    data : dict
        A dictionary of lists of values, e.g a list of apd30
    x : float
        The number of standard deviations to be
        included, ex x = 1.0. If none is provided
        everything will be included

    Returns
    -------
    new_data : dict
        A dictionary with new lists containing only those elements
        that are within :math:`x` std of the mean
    included_indice : dict
        Indices included for each (key, value) pair.
    means : array
        The mean of the elements in ``new_data``


    """
    if skips is None:
        skips = []
    # List with new data
    new_data = {k: [] for k in data.keys()}
    included_indices = {k: [] for k in data.keys()}

    for k, v in data.items():
        if x is not None:
            for j, s in enumerate(v):

                # Include only signals within factor *
                # standard deviation from the mean
                if -np.std(v) * x < (s - np.mean(v)) < np.std(v) * x:

                    new_data[k].append(s)
                    included_indices[k].append(j)

        # Check if list is empty
        if not new_data[k]:
            if x is not None and k not in skips:
                msg = (
                    "No data were within {} std for key {}" ". Include all data"
                ).format(k, x)
                logger.info(msg)
            # If it is empty lets include everything
            new_data[k] = v
            included_indices[k] = range(len(v))

    excluded_data = namedtuple(
        "excluded_data", ("new_data, included_indices, " "all_included_indices")
    )
    intsect = utils.get_intersection(included_indices)
    all_included_indices = list(intsect)
    return excluded_data(
        new_data=new_data,
        included_indices=included_indices,
        all_included_indices=all_included_indices,
    )


def analyze_mps_func(mps_data, mask=None, **kwargs):
    avg = average_intensity(mps_data.frames, mask=mask)
    time_stamps = mps_data.time_stamps
    pacing = mps_data.pacing
    info = mps_data.info
    metadata = mps_data.metadata

    analyzer = AnalyzeMPS(
        avg=avg,
        time_stamps=time_stamps,
        pacing=pacing,
        info=info,
        metadata=metadata,
        **kwargs,
    )
    analyzer.analyze_all()
    return analyzer.data


def analyze_frequencies(
    chopped_data: List[List[float]],
    chopped_times: List[List[float]],
    time_unit: str = "ms",
    fname: str = "",
    plot=True,
) -> np.ndarray:

    freqs = beating_frequency_modified(chopped_data, chopped_times, time_unit)

    mean_freq = np.median(freqs)
    std_freq = np.std(freqs)

    is_significant = np.abs(freqs - mean_freq) > std_freq
    logger.info(
        f"Found {sum(is_significant)} significant beats with regard to beat frequency"
    )

    if plot:

        try:
            import matplotlib.pyplot as plt
        except ImportError:
            return freqs

        fig, ax = plt.subplots(2, 1, figsize=(15, 10))
        x = np.arange(len(freqs)) + 1
        ax[0].bar(x, freqs)
        for sig, patch in zip(is_significant, ax[0].patches):
            if not sig:
                continue
            ax[0].text(
                patch.get_x() + patch.get_width() / 2.0,
                patch.get_height(),
                "*",
                ha="center",
            )

        # ax.plot(x, freqs, "k-o")
        ax[0].set_ylabel("Frequency [Hz]")
        ax[0].set_xlabel("Beat number")

        points_x = [t[int(np.argmax(c))] for c, t in zip(chopped_data, chopped_times)]
        points_y = [np.max(c) for c in chopped_data]

        for c, t in zip(chopped_data, chopped_times):
            ax[1].plot(t, c)

        ax[1].plot(points_x, points_y, linestyle="", marker="o", color="r")
        ax[1].set_xlabel("Time [ms]")
        ax[1].set_ylabel(r"$\Delta F / F$")
        if fname != "":
            fig.savefig(fname)
        plt.close()

    return freqs


def detect_ead(
    y: Union[List[float], np.ndarray], sigma: float = 3, prominence_level: float = 0.04
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
        has_ead, peak_index = detect_ead(c)
        if has_ead and peak_index is not None:
            num_eads += 1
            peak_inds[i] = peak_index

    if plot:
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            return num_eads

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

    return num_eads


def poincare_plot(
    chopped_data: List[List[float]],
    chopped_times: List[List[float]],
    apds: List[int],
    fname: str = "",
) -> Dict[int, List[float]]:
    """
    Create poincare plots for given APDs

    Arguments
    ---------
    chopped_data : list
        List of the amplitude of each beat
    chopped_data : list
        List of the time stamps of each beat
    apds : list
        List of APDs to be used, e.g [30, 50, 80]
        will plot the APS30, APD50 and APD80.
    fname : str
        Path to filename to save the figure.

    Notes
    -----
    For more info see <http://doi.org/10.1371/journal.pcbi.1003202>

    Returns
    -------
    dict:
        A dictionary with the key being the :math:`x` in APD:math:`x`
        and the value begin the points being plotted.

    """
    apds_points = {
        k: [apd(y, k / 100, t) for y, t in zip(chopped_data, chopped_times)]
        for k in apds
    }

    if len(chopped_data) <= 1:
        # No possible to plot poincare plot with 1 or zero elements
        return apds_points

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return apds_points

    fig, ax = plt.subplots()
    for k, v in apds_points.items():
        ax.plot(v[:-1], v[1:], label=f"APD{k}", marker=".")
    ax.legend()
    ax.grid()
    ax.set_xlabel("APD(n-1)[ms]")
    ax.set_ylabel("APD(n) [ms]")
    if fname != "":
        fig.savefig(fname)
    else:
        plt.show()
    plt.close()
    return apds_points


class AnalyzeMPS:
    def __init__(
        self,
        avg: List[float],
        time_stamps: List[float],
        pacing: Optional[List[float]] = None,
        **kwargs,
    ):

        self.avg = avg
        self.time_stamps = time_stamps
        self.pacing = pacing
        self._handle_arguments(kwargs)

        self.unchopped_data: Dict[str, Any] = {}
        self.chopped_data: Dict[str, Any] = {}
        self.features: Dict[str, Any] = {}
        self._all_features: Dict[str, Any] = {}
        self.features_computed = False

        self.unchopped_data["original_trace"] = np.copy(avg)
        self.unchopped_data["original_times"] = np.copy(time_stamps)
        self.unchopped_data["original_pacing"] = np.copy(pacing)

    @property
    def chopping_parameters(self):
        return dict(
            threshold_factor=self.parameters["chopping_threshold_factor"],
            extend_front=self.parameters["chopping_extend_front"],
            extend_end=self.parameters["chopping_extend_end"],
            min_window=self.parameters["chopping_min_window"],
            max_window=self.parameters["chopping_max_window"],
        )

    @property
    def data(self):
        return {
            "chopped_data": self.chopped_data,
            "unchopped_data": self.unchopped_data,
            "features": self.features,
            "chopping_parameters": self.chopping_parameters,
            "attributes": self.info,
        }

    def dump_data(self, dump_all=False):
        if self.outdir is None:
            return

        utils.dump_data(self.data, self.outdir.joinpath("data"), "npy")
        if dump_all:
            chopped_data_padded = bin_utils.padding(self.chopped_data)
            unchopped_data_padded = bin_utils.padding(self.unchopped_data)
            utils.to_csv(chopped_data_padded, self.outdir.joinpath("chopped_data"))
            utils.to_csv(unchopped_data_padded, self.outdir.joinpath("unchopped_data"))

            with open(self.outdir.joinpath("data.txt"), "w") as f:
                f.write(self.info_txt)

            # utils.dump_data(data, self.outdir.joinpath("data"), "mat")
            with open(self.outdir.joinpath("metadata.json"), "w") as f:
                json.dump(self.metadata, f, indent=4, default=bin_utils.json_serial)

            about_str = bin_utils.about()
            with open(self.outdir.joinpath("ABOUT.md"), "w") as f:
                f.write(about_str)

    def analyze_all(self):
        self.analyze_unchopped_data()
        self.analyze_chopped_data()
        self.dump_data(dump_all=True)

    def analyze_unchopped_data(self):

        self.dump_data()

        if self.parameters["plot"]:
            plotter.plot_single_trace(
                self.time_stamps, self.avg, self.outdir.joinpath("original_trace")
            )

        self.remove_spikes()
        self.filter()

        if self.parameters["plot"]:
            plotter.plot_twin_trace(
                self.time_stamps,
                self.time_stamps,
                self.avg,
                self.pacing,
                self.outdir.joinpath("pacing"),
            )

        self.ignore_pacing()
        self.remove_points()
        self.background_correction()
        self.dump_data()

    def _set_choopped_data(self):
        self._chopped_data = chop_data(
            self.avg, self.time_stamps, pacing=self.pacing, **self.chopping_parameters
        )
        if len(self._chopped_data.data) == 0:
            self._chopped_data.data.append(self.avg)
            self._chopped_data.times.append(self.time_stamps)
            self._chopped_data.pacing.append(self.pacing)

        for i, (t, d, p) in enumerate(
            zip(
                self._chopped_data.times,
                self._chopped_data.data,
                self._chopped_data.pacing,
            )
        ):
            self.chopped_data["time_{}".format(i)] = t
            self.chopped_data["trace_{}".format(i)] = d
            self.chopped_data["pacing_{}".format(i)] = p

    def analyze_chopped_data(self):
        self._set_choopped_data()

        if self.parameters["plot"]:
            poincare_plot(
                self._chopped_data.data,
                self._chopped_data.times,
                apds=[30, 50, 70, 80, 90],
                fname=self.outdir.joinpath("poincare_plot"),
            )

        self.compute_features()
        self.run_ead_and_apd_analysis()
        self.filter_chopped_data()
        self.get_pacing_avg()
        self.get_average_all()
        self.get_average_1std()
        self.plot_chopped_averages()
        self.dump_data()

    def run_ead_and_apd_analysis(self):
        if not hasattr(self, "_chopped_data"):
            self._set_choopped_data()

        apd_analysis = analyze_apds(
            self._chopped_data.data,
            self._chopped_data.times,
            plot=self.parameters["plot"],
            max_allowed_apd_change=self.parameters["max_allowed_apd_change"],
            fname="" if self.outdir is None else self.outdir.joinpath("apd_analysis"),
        )
        self.features["slope_APD80"] = apd_analysis.slope_APD80
        self.features["slope_cAPD80"] = apd_analysis.slope_cAPD80
        self._all_features["triangulation"] = apd_analysis.triangulation

        analyze_frequencies(
            self._chopped_data.data,
            self._chopped_data.times,
            plot=self.parameters["plot"],
            time_unit=self.info.get("time_unit", "ms"),
            fname=""
            if self.outdir is None
            else self.outdir.joinpath("frequency_analysis"),
        )

        self.features["num_eads"] = analyze_eads(
            self._chopped_data.data,
            self._chopped_data.times,
            sigma=self.parameters["ead_sigma"],
            prominence_threshold=self.parameters["ead_prominence_threshold"],
            fname="" if self.outdir is None else self.outdir.joinpath("EAD_analysis"),
        )

    @property
    def info_txt(self):
        if not hasattr(self, "_features_all"):
            self.filter_chopped_data()

        info_txt = ["{0:20}\t{1:>20}".format(k, v) for k, v in self.info.items()]

        return "\n".join(
            np.concatenate(
                (
                    ["All features"],
                    self._features_all,
                    ["\n"],
                    ["Features within 1std"],
                    self._features_1std,
                    ["\n"],
                    ["Info"],
                    info_txt,
                )
            )
        )

    def compute_features(self):
        if not hasattr(self, "_chopped_data"):
            self._set_choopped_data()

        freqs = beating_frequency_modified(
            self._chopped_data.data,
            self._chopped_data.times,
            self.info.get("time_unit", "ms"),
        )
        self.features["beating_frequency"] = np.mean(freqs)

        self._all_features.update(
            compute_features(
                self._chopped_data,
                use_spline=self.parameters["use_spline"],
                normalize=self.parameters["normalize"],
            )
        )
        for k in [30, 50, 80, 90]:
            self._all_features[f"capd{k}"] = corrected_apd(
                self._all_features[f"apd{k}"], 60 / self.features["beating_frequency"]
            )

    def compute_mean_features(self):
        if not hasattr(self, "_all_features"):
            self.compute_features()

        self._excluded = exclude_x_std(
            self._all_features, self.parameters["std_ex"], skips=["upstroke80"]
        )

        mean_features = {}
        self._features_1std = []
        self._features_all = []
        template = "{0:20}\t{1:10.1f}  +/-  {2:10.3f}"
        for k in self._all_features.keys():

            v = self._excluded.new_data[k]
            self._features_1std.append(template.format(k, np.mean(v), np.std(v)))
            mean_features[k] = np.mean(v)

            v = self._all_features[k]
            self._features_all.append(template.format(k, np.mean(v), np.std(v)))
        self.features.update(mean_features)

        self.features_computed = True

    def filter_chopped_data(self):
        if not self.features_computed:
            self.compute_mean_features()

        logger.debug("Plot chopped data")
        titles = ["apd30", "apd50", "apd80", "dFdt_max", "int30", "tau75"]
        xs = []
        ys = []
        for k in titles:
            logger.debug(k)
            x = [
                np.subtract(xi, xi[0])
                for xi in [
                    t
                    for j, t in enumerate(self._chopped_data.times)
                    if j in self._excluded.included_indices[k]
                ]
            ]
            y = [
                d
                for j, d in enumerate(self._chopped_data.data)
                if j in self._excluded.included_indices[k]
            ]

            xs.append(x)
            ys.append(y)

        if self.parameters["plot"]:
            plotter.plot_multiple_traces(
                xs,
                ys,
                self.outdir.joinpath("chopped_data_features"),
                titles,
                ylabels=[r"$\Delta F / F$" for _ in xs],
                deep=True,
            )

        self._xs_all = [np.subtract(xi, xi[0]) for xi in self._chopped_data.times]
        self._ys_all = deepcopy(self._chopped_data.data)
        self.features["num_beats"] = len(self._xs_all)

    def get_average_all(self):
        if not hasattr(self, "_xs_all"):
            self.filter_chopped_data()

        N = self.parameters["N"]
        logger.debug("Get average of everything")
        # Average of everything
        if self.parameters["spike_duration"] > 0:
            self.chopped_data["trace_all"] = average.get_subsignal_average(
                self._chopped_data.data
            )
            idx = np.argmax([len(xi) for xi in self._xs_all])
            self.chopped_data["time_all"] = self._xs_all[idx]
        else:
            (
                self.chopped_data["trace_all"],
                self.chopped_data["time_all"],
            ) = average.get_subsignal_average_interpolate(
                self._chopped_data.data, self._xs_all, N
            )

    def get_pacing_avg(self):
        if not hasattr(self, "_xs_all"):
            self.filter_chopped_data()

        N = self.parameters["N"]

        if self.parameters["spike_duration"] > 0:
            idx = np.argmax([len(xi) for xi in self._xs_all])
            pacing_avg = self._chopped_data.pacing[idx]

        else:
            pacing_avg = np.interp(
                np.linspace(
                    np.min(self._chopped_data.times[0]),
                    np.max(self._chopped_data.times[0]),
                    N,
                ),
                self._chopped_data.times[0],
                self._chopped_data.pacing[0],
            )
        pacing_avg[pacing_avg <= 2.5] = 0.0
        pacing_avg[pacing_avg > 2.5] = 5.0
        self.chopped_data["pacing_all"] = pacing_avg
        self.chopped_data["pacing_1std"] = pacing_avg

    def get_average_1std(self):
        if not hasattr(self, "_xs_all"):
            self.filter_chopped_data()

        N = self.parameters["N"]
        # Average of everyting within 1 std
        logger.debug("Get average of everything within 1 std")

        if len(self._excluded.all_included_indices) == 0:
            self._xs_1std, self._ys_1std = self._xs_all, self._ys_all
        else:
            xs_1std_ = np.array(self._chopped_data.times, dtype=object)[
                self._excluded.all_included_indices
            ]
            self._xs_1std = [np.subtract(xi, xi[0]) for xi in xs_1std_]
            self._ys_1std = np.array(self._chopped_data.data, dtype=object)[
                self._excluded.all_included_indices
            ]

        if self.parameters["spike_duration"] > 0:
            self.chopped_data["trace_1std"] = average.get_subsignal_average(
                self._ys_1std
            )
            self.chopped_data["time_1std"] = self._xs_1std[
                np.argmax([len(xi) for xi in self._xs_1std])
            ]
        else:
            (
                self.chopped_data["trace_1std"],
                self.chopped_data["time_1std"],
            ) = average.get_subsignal_average_interpolate(
                self._ys_1std, self._xs_1std, N
            )

    def plot_chopped_averages(self):

        if not self.parameters["plot"]:
            return
        if not hasattr(self, "_xs_1st"):
            self.get_average_1std()

        logger.debug("Plot global chopped average")

        plotter.plot_multiple_traces(
            [self._xs_all, self._xs_1std],
            [self._ys_all, self._ys_1std],
            self.outdir.joinpath("chopped_data"),
            ["all", "1 std"],
            ylabels=[r"$\Delta F / F$" for _ in range(len(self._xs_all))],
            deep=True,
        )

        plotter.plot_multiple_traces(
            [self.chopped_data["time_all"], self.chopped_data["time_1std"]],
            [self.chopped_data["trace_all"], self.chopped_data["trace_1std"]],
            self.outdir.joinpath("average"),
            ["all", "1 std"],
            ylabels=[r"$\Delta F / F$" for _ in range(len(self._xs_all))],
        )

        plotter.plot_twin_trace(
            self.chopped_data["time_1std"],
            self.chopped_data["time_1std"],
            self.chopped_data["trace_1std"],
            self.chopped_data["pacing_1std"],
            self.outdir.joinpath("average_pacing"),
        )

    def background_correction(self):
        if self.parameters["background_correction"]:

            logger.debug("Apply background correction")

            background_ = background(self.time_stamps, self.avg)
            bkgplot_y = np.transpose([np.copy(self.avg), background_])
            bkgplot_x = np.transpose([self.time_stamps, self.time_stamps])
            self.avg = correct_background(self.time_stamps, self.avg)

            logger.debug("Plot corrected trace trace")
            if self.parameters["plot"]:
                plotter.plot_multiple_traces(
                    [bkgplot_x, self.time_stamps],
                    [bkgplot_y, self.avg],
                    self.outdir.joinpath("corrected_trace"),
                    ylabels=["Pixel intensity", r"$\Delta F / F$"],
                    titles=["background", "corrected"],
                )
        else:
            background_ = np.zeros_like(self.avg)

        self.unchopped_data["trace"] = np.copy(self.avg)
        self.unchopped_data["time"] = np.copy(self.time_stamps)
        self.unchopped_data["pacing"] = np.copy(self.pacing)
        self.unchopped_data["background"] = background_

    def remove_points(self):

        for p in self.parameters["remove_points_list"]:
            logger.debug(f"Remove points: '{p}'")
            t_ = np.copy(self.time_stamps)
            self.time_stamps, self.avg = remove_points(self.time_stamps, self.avg, *p)
            _, self.pacing = remove_points(t_, self.pacing, *p)

    def ignore_pacing(self):
        if self.parameters["ignore_pacing"]:
            logger.debug("Ignore pacing")
            self.pacing = np.multiply(self.pacing, 0.0)

    def filter(self):
        if self.parameters["filter"]:
            logger.debug("Filter signal")
            self.avg = filt(self.avg)

    def remove_spikes(self):
        if self.parameters["spike_duration"] == 0:
            return

        sd = self.parameters["spike_duration"]
        logger.debug(f"Remove spikes with spike duration {sd}")
        self.avg = remove_spikes(self.avg, self.pacing, sd)
        self.time_stamps = remove_spikes(self.time_stamps, self.pacing, sd)
        self.pacing = remove_spikes(self.pacing, self.pacing, sd)

        self.unchopped_data["original_trace_no_spikes"] = np.copy(self.avg)
        self.unchopped_data["original_times_no_spikes"] = np.copy(self.time_stamps)
        self.unchopped_data["original_pacing_no_spikes"] = np.copy(self.pacing)

    def _handle_arguments(self, kwargs):

        parameters = AnalyzeMPS.default_parameters()
        parameters.update(kwargs)
        self.info = parameters.pop("info")
        self.metadata = parameters.pop("metadata")
        self.outdir = parameters["outdir"]
        if self.outdir is not None:
            self.outdir = Path(self.outdir)
            self.outdir.mkdir(exist_ok=True, parents=True)
        else:
            parameters["plot"] = False
        self.parameters = parameters

    @staticmethod
    def default_parameters():
        """Default parametert for the analysis script.
        Current default values are (with types)

        .. code::

            spike_duration: int = 0
            filter: bool = True
            ignore_pacing: bool = False
            remove_points_list: tuple = ()
            chopping_threshold_factor: float = 0.3
            chopping_extend_front: Optional[float] = None
            chopping_extend_end: Optional[float] = None
            chopping_min_window: float = 200
            chopping_max_window: float = 2000
            std_ex: float = 1.0
            use_spline: bool = True
            normalize: bool = False
            outdir: str = ""
            plot: bool = False
            bar=None
            info=None
            metadata=None
            background_correction=True
            max_allowed_apd_change=None
            ead_sigma: float = 3
            ead_prominence_threshold: float = 0.04
            N=200

        Returns
        -------
        dict
            Dictionary with default parameters
        """

        return dict(
            spike_duration=0,
            filter=True,
            ignore_pacing=False,
            remove_points_list=(),
            chopping_threshold_factor=0.3,
            chopping_extend_front=None,
            chopping_extend_end=None,
            chopping_min_window=200,
            chopping_max_window=2000,
            std_ex=1.0,
            use_spline=True,
            normalize=False,
            outdir=None,
            plot=False,
            bar=None,
            info={},
            metadata={},
            background_correction=True,
            max_allowed_apd_change=None,
            ead_sigma=3,
            ead_prominence_threshold=0.04,
            N=200,
        )


def frame2average(frame, times=None, normalize=True, background_correction=True):
    """
    Compute average pixel intensity of the frames

    Arguments
    ---------
    frames : np.ndarray
        The frames
    times : np.ndarray
        The time stamps
    normalize : bool
        In true normalize average so that it takes value
        between zero and one. Default: True
    background_correction : bool
        If true, apply backround correction to the average.
        You don't want to do this if you only have a single beat.
        Default: True.

    Returns
    -------
    trace : np.ndarray
        The average trace

    """

    avg_ = average.get_average_all(frame)
    if background_correction:
        assert (
            times is not None
        ), "Please provide time stamps for background correection"
        avg = correct_background(times, avg_)
    else:
        avg = avg_

    if normalize:
        trace = normalize_signal(avg)
    else:
        trace = avg

    return trace


def local_averages(
    frames,
    times=None,
    N=10,
    x_start=0,
    y_start=0,
    x_end=None,
    y_end=None,
    background_correction=True,
    normalize=False,
    loglevel=logging.INFO,
    **kwargs,
):
    """
    Compute the local averages

    Arguments
    ---------
    frames : np.ndarray
        The image sequence on the form (nx, ny, T) where nx
        and ny are repspectively the number of pixels in the
        x and y direction and T are the number of frames.
    times : np.ndarray or list
        An array of time stamps.
    N : int
        Maximum number of grid points. The axis with the greatest number
        of pixels will be partitioned into a coarser grid of size n. The
        other axis will be scaled so that each grid point is approximately
        square.
    x_start : int
        Index where to start in x-direction
    x_end : int
        Index where to end in x-direction
    y_start : int
        Index where to start in y-direction
    y_end : int
        Index where to end in y-direction
    backround_correction : bool
        If you want to apply background correction. You typically want
        to allways do this except in the case when you only have
        a single beat. Default: True.
    normalize : bool
        If True, normalize all averages to have values between 0 and 1, if
        False keep the original values. Default: False
    loglevel : int
        Verbosity. Default: INFO (=20). For more info see the logging library.
    """

    logger = utils.get_logger(__name__, loglevel)
    logger.debug("Compute local averages")
    grid = utils.get_grid_settings(
        N, x_start=x_start, y_start=y_start, x_end=x_end, y_end=y_end, frames=frames
    )

    futures = np.empty((grid.nx, grid.ny), dtype=object)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for i in range(grid.nx):
            for j in range(grid.ny):

                x0 = x_start + i * grid.dx
                x1 = min(x0 + grid.dx, grid.x_end)

                y0 = y_start + j * grid.dy
                y1 = min(y0 + grid.dy, grid.y_end)

                if y0 >= y1:
                    continue
                if x0 >= x1:
                    continue

                logger.debug(f"x0 = {x0}, x1 = {x1}, y0 = {y0}, y1 = {y1}")
                im = frames[x0:x1, y0:y1, :]
                kwargs = dict(
                    frame=im,
                    times=times,
                    normalize=normalize,
                    background_correction=background_correction,
                )
                futures[i, j] = executor.submit(_frames2average, kwargs)

    local_averages = np.zeros((grid.nx, grid.ny, len(times)))
    for i in range(grid.nx):
        for j in range(grid.ny):
            local_averages[i, j, :] = futures[i, j].result()

    return local_averages


def _frames2average(kwargs):
    return frame2average(**kwargs)


def analyze_local_average(frames, times=None, mask=None, N=10):

    loc = local_averages(frames, times, N=N)
    avg = frame2average(
        frame=frames, times=times, normalize=False, background_correction=True
    )
    avg_std = np.std(avg)
    if mask is None:
        # We include everything
        mask = np.ones((loc.shape[0], loc.shape[1]), dtype=bool)
    new_mask = np.zeros((loc.shape[0], loc.shape[1]), dtype=bool)
    # import matplotlib as mpl
    # import matplotlib.pyplot as plt
    # fig, ax = plt.subplots()
    for i in range(loc.shape[0]):
        for j in range(loc.shape[1]):
            local = loc[i, j, :]

            if mask[i, j] and np.std(local) > avg_std:
                new_mask[i, j] = True

                # ax.plot(times, loc[i, j, :])
    # ax.plot(times, avg, color="k")
    # plt.show()
    # grid = utils.get_grid_settings(N, frames=frames)

    # fig, ax = plt.subplots()
    # ax.imshow(frames.T[0])

    # for i in range(grid.nx):
    #     for j in range(grid.ny):

    #         facecolor = "red" if mask[i, j] else "yellow"
    #         p = mpl.patches.Rectangle(
    #             (i * grid.dx, j * grid.dy),
    #             grid.dx,
    #             grid.dy,
    #             linewidth=1,
    #             edgecolor="b",
    #             facecolor=facecolor,
    #             alpha=0.2,
    #         )
    #         ax.add_patch(p)
    # plt.show()
    # from IPython import embed

    # embed()
    # exit()
    return loc, new_mask


def baseline_intensity(
    frames,
    times,
    N=10,
    x_start=0,
    y_start=0,
    x_end=None,
    y_end=None,
    normalize=False,
    loglevel=logging.INFO,
    **kwargs,
):
    """
    Compute the baseline intensity in local windows

    Arguments
    ---------
    frames : np.ndarray
        The image sequence on the form (nx, ny, T) where nx
        and ny are repspectively the number of pixels in the
        x and y direction and T are the number of frames.
    times : np.ndarray or list
        An array of time stamps.
    N : int
        Maximum number of grid points. The axis with the greatest number
        of pixels will be partitioned into a coarser grid of size n. The
        other axis will be scaled so that each grid point is approximately
        square.
    x_start : int
        Index where to start in x-direction
    x_end : int
        Index where to end in x-direction
    y_start : int
        Index where to start in y-direction
    y_end : int
        Index where to end in y-direction
    backround_correction : bool
        If you want to apply background correction. You typically want
        to allways do this except in the case when you only have
        a single beat. Default: True.
    normalize : bool
        If True, normalize all averages to have values between 0 and 1, if
        False keep the original values. Default: False
    loglevel : int
        Verbosity. Default: INFO (=20). For more info see the logging library.
    """

    loc = local_averages(
        frames=frames,
        times=times,
        N=N,
        x_start=x_start,
        y_start=y_start,
        x_end=x_end,
        y_end=y_end,
        background_correction=False,
        normalize=normalize,
        loglevel=loglevel,
    )
    shape = loc.shape[:2]
    baseline_intensity = np.zeros((shape[0], shape[1], len(times)))
    for i in range(shape[0]):
        for j in range(shape[1]):
            loc_ij = loc[i, j, :]
            baseline_intensity[i, j, :] = background(times, loc_ij)

    return baseline_intensity


def prevalence(
    mps_data,
    snr_factor=1.5,
    N=50,
    frequency_threshold=0.2,
    baseline_threshold=0.1,
    **kwargs,
):
    """
    Compute the prevalence, i.e the percentage of living
    cells in the recording

    Arguments
    ---------
    mps_data : mps.load.MPS
        The data
    snr_factor : float
        Factor multiplied with the global signal to noise ratio (snr),
        to determine wether a region is noise or not. If a local region
        has a larger values than snr_factor * global snr it will be
        classied as noise. Default: 1.5
    N : int
        Size of grid along the major axis (minor axis will
        be scaled proportionally). Default: 50
    frequency_threshold : float
        Percentage (between 0 and 1) of how far from the global
        frequency a local signal can be before it is classified
        as non-beating.
        Default: 0.2
    baseline_threshold : float
        Percentage (between 0 and 1) of how far from the min and
        max values of the baseline intensities of the beating regions
        that should define the rage for baseline intensity of tissue
        Default: 0.1

    Returns
    -------
    mps_prevalence : namedtuple
        A named tuple with the following fields
    prevalence: float
        percentage of tissue that is beating
    tissue_covered_area : float
        percentage of area in the image that is
        classified as tissue
    is_beating : array
        boolean array who's true values are classied
        as beating regions
    is_tissue : array
        boolean array who's true values are classied
        as regions with tissue

    """

    # Get the local traces
    loc = local_averages(
        frames=mps_data.frames, times=mps_data.time_stamps, N=N, **kwargs
    )

    avg = average.get_average_all(mps_data.frames)
    chopped_data = chop_data_without_pacing(avg, mps_data.time_stamps)
    global_freq = beating_frequency(
        chopped_data.times, unit=mps_data.info.get("time_unit", "ms")
    )
    logger.info(f"Global frequency: {global_freq}")
    global_snr = snr(correct_background(mps_data.time_stamps, avg))
    logger.info(f"Global SNR: {global_snr}")

    # Find the beating frequency for each local average
    shape = loc.shape[:2]
    freq = np.zeros(shape)
    local_snr = np.zeros(shape)
    is_tissue = np.ones(shape, dtype=bool)
    for i, j in it.product(range(shape[0]), range(shape[1])):
        loc_ij = loc[i, j, :]
        chopped_data = chop_data_without_pacing(loc_ij, mps_data.time_stamps)
        freq[i, j] = beating_frequency(
            chopped_data.times, unit=mps_data.info.get("time_unit", "ms")
        )
        local_snr[i, j] = snr(loc[i, j, :])

    # Set non-beating regions where the local signal to
    # noise ratio is high
    is_beating = local_snr <= snr_factor * global_snr

    # Frequencies that deviates too much from the
    # global frequency is most likely noise
    freq_dev = np.where(
        ~np.isclose(freq, global_freq, frequency_threshold * global_freq)
    )
    is_beating[freq_dev] = False

    # The the baselines that are normally subtracted
    baselines = baseline_intensity(
        frames=mps_data.frames, times=mps_data.time_stamps, N=N, **kwargs
    )

    baseline_min = baselines.min(-1)[is_beating]
    baseline_max = baselines.max(-1)[is_beating]
    # Let the range be the the between the 10% lowest
    # and 90% highest just to remove outliers
    baseline_range = (
        sorted(baseline_min)[int(len(baseline_min) * baseline_threshold)],
        sorted(baseline_max)[int(len(baseline_min) * (1 - baseline_threshold))],
    )

    # Check the baseline values in the non-beating tissue
    for i, j in zip(*np.where(~is_beating)):

        # If the baselines values are outside the range
        # we will classify it at non-tissue
        if (
            baselines[i, j, :].min() < baseline_range[0]
            or baselines[i, j, :].max() > baseline_range[1]
        ):
            is_tissue[i, j] = False

    return mps_prevalence(
        prevalence=is_beating.sum() / is_tissue.sum(),
        tissue_covered_area=is_tissue.sum() / np.multiply(*is_tissue.shape),
        is_beating=is_beating,
        is_tissue=is_tissue,
    )
