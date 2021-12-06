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
from textwrap import dedent
from typing import Any, Dict, List, Optional

import ap_features as apf
import numpy as np
from scipy.interpolate import UnivariateSpline

from . import average, ead, plotter, utils

logger = utils.get_logger(__name__)

mps_prevalence = namedtuple(
    "mps_prevalence",
    "prevalence, tissue_covered_area, is_beating, is_tissue",
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
        frames,
        alpha=cutoff_factor,
        return_indices=True,
    )

    mask = np.ones_like(frames.T[0], dtype=bool).reshape(-1)
    mask[inds] = False
    mask = mask.reshape(frames.T[0].T.shape)
    return mask


def average_intensity(data, mask=None, alpha=1.0, averaging_type="spatial"):
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

    if alpha == 1.0:
        # Mean of everything
        if mask is None:
            avg = average.get_average_all(data)
        else:
            avg = average.masked_average(data, mask)

    else:

        if averaging_type == "spatial":
            avg = average.get_spatial_average(data, alpha=alpha)

        elif averaging_type == "temporal":
            avg = average.get_temporal_average(data, alpha=alpha)
        else:
            msg = (
                "Unknown averaging_type {}. Expected averaging type to "
                'be one of ["all", "spatial", "temporal"]'
            )
            logger.error(msg)
            raise ValueError(msg)

    return avg


def analyze_apds(
    chopped_data: List[List[float]],
    chopped_times: List[List[float]],
    max_allowed_apd_change: Optional[float] = None,
    fname: str = "",
    plot=True,
) -> APDAnalysis:

    apd_levels = (100 * np.sort(np.append(np.arange(0.1, 0.91, 0.2), 0.8))).astype(int)
    apds = {
        k: [apf.features.apd(k, y, t) for y, t in zip(chopped_data, chopped_times)]
        for k in apd_levels
    }

    apd_points = {
        k: [apf.features._apd(k, y, t) for y, t in zip(chopped_data, chopped_times)]
        for k in apd_levels
    }

    median_freq = np.median(
        apf.features.beating_frequency_from_peaks(chopped_data, chopped_times),
    )
    beat_rates = {}
    apd_dt = {}
    for k, apdx in apd_points.items():
        apd_dt[k] = [v[0] for v in apdx]
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
        k: apf.features.corrected_apd(v, b).tolist()
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
        apd30 = apf.features._apd(30, y, t)
        apd80 = apf.features._apd(80, y, t)

        tri = apd80[-1] - apd30[-1]
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
            chopped_data=chopped_data,
            chopped_times=chopped_times,
            res=res,
            fname=fname,
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
        np.array(list(res.is_significant.values())).flatten(),
        ax[0].patches,
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
    for i, (label, vals) in enumerate(res.apd_points.items()):
        x = [v[0] for v in vals] + [v[1] for v in vals]

        y = [
            UnivariateSpline(chopped_times[j], chopped_data[j])(v[0])
            for j, v in enumerate(vals)
        ] + [
            UnivariateSpline(chopped_times[j], chopped_data[j])(v[1])
            for j, v in enumerate(vals)
        ]
        ax[2].plot(x, y, linestyle="", marker="o", label=label)

    ax[2].legend()
    ax[2].grid()
    ax[2].set_xlabel("Time [ms]")
    ax[2].set_ylabel(r"$\Delta F / F$")

    if fname != "":
        fig.savefig(fname)
    plt.close()


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
        unit = apf.utils.time_unit(time)
        unitfactor = 1000.0 if unit == "s" else 1.0

        time = np.multiply(unitfactor, time)

        apd90[i] = apf.features.apd(90, data, time, use_spline=use_spline) or np.nan
        apd80[i] = apf.features.apd(80, data, time, use_spline=use_spline) or np.nan
        apd50[i] = apf.features.apd(50, data, time, use_spline=use_spline) or np.nan
        apd30[i] = apf.features.apd(30, data, time, use_spline=use_spline) or np.nan
        tau75[i] = apf.features.tau(time, data, 0.75)
        upstroke80 = apf.features.upstroke(time, data, 0.8)

        dFdt_max[i] = (
            apf.features.maximum_upstroke_velocity(
                data,
                time,
                use_spline=use_spline,
                normalize=normalize,
            )
            or np.nan
        )
        int30[i] = (
            apf.features.integrate_apd(
                data,
                time,
                0.3,
                use_spline=use_spline,
                normalize=normalize,
            )
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
        "excluded_data",
        ("new_data, included_indices, " "all_included_indices"),
    )
    intsect = utils.get_intersection(included_indices)
    all_included_indices = list(intsect)
    return excluded_data(
        new_data=new_data,
        included_indices=included_indices,
        all_included_indices=all_included_indices,
    )


def analyze_frequencies(
    chopped_data: List[List[float]],
    chopped_times: List[List[float]],
    time_unit: str = "ms",
    fname: str = "",
    plot=True,
) -> np.ndarray:

    freqs = apf.features.beating_frequency_from_peaks(
        chopped_data,
        chopped_times,
        time_unit,
    )

    mean_freq = np.median(freqs)
    std_freq = np.std(freqs)

    is_significant = np.abs(freqs - mean_freq) > std_freq
    logger.info(
        f"Found {sum(is_significant)} significant beats with regard to beat frequency",
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
        k: [apf.features.apd(k, y, t) for y, t in zip(chopped_data, chopped_times)]
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


def analyze_mps_func(
    mps_data, mask=None, analysis_window_start=0, analysis_window_end=-1, **kwargs
):
    avg = average_intensity(mps_data.frames, mask=mask)
    time_stamps = mps_data.time_stamps
    pacing = mps_data.pacing

    info = mps_data.info
    metadata = mps_data.metadata

    start_index = 0
    end_index = len(time_stamps)
    if analysis_window_start > 0:
        try:
            start_index = next(
                i for i, v in enumerate(time_stamps) if v >= analysis_window_start
            )
        except StopIteration:
            pass
    if analysis_window_end != -1:
        try:
            end_index = next(
                i for i, v in enumerate(time_stamps) if v >= analysis_window_end
            )
        except StopIteration:
            pass

    avg = avg[start_index:end_index]
    time_stamps = time_stamps[start_index:end_index]
    pacing = pacing[start_index:end_index]

    analyzer = AnalyzeMPS(
        avg=avg,
        time_stamps=time_stamps,
        pacing=pacing,
        info=info,
        metadata=metadata,
        **kwargs,
    )
    analyzer.analyze_all()
    try:
        import matplotlib.pyplot as plt

        plt.close("all")
    finally:
        return analyzer.data


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
        self.intervals: List[apf.chopping.Intervals] = []
        self.features: Dict[str, Any] = {}
        self._all_features: Dict[str, Any] = {}
        self.features_computed = False
        self.upstroke_times: List[float] = []

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
            "intervals": self.intervals,
            "upstroke_times": self.upstroke_times,
        }

    def dump_data(self, dump_all=False):
        if self.outdir is None:
            return

        utils.dump_data(self.data, self.outdir.joinpath("data"), "npy")
        if dump_all:
            chopped_data_padded = utils.padding(self.chopped_data)
            unchopped_data_padded = utils.padding(self.unchopped_data)
            utils.to_csv(chopped_data_padded, self.outdir.joinpath("chopped_data"))
            utils.to_csv(unchopped_data_padded, self.outdir.joinpath("unchopped_data"))

            with open(self.outdir.joinpath("data.txt"), "w") as f:
                f.write(self.info_txt)

            # utils.dump_data(data, self.outdir.joinpath("data"), "mat")
            with open(self.outdir.joinpath("metadata.json"), "w") as f:
                json.dump(self.metadata, f, indent=4, default=utils.json_serial)

            about_str = AnalyzeMPS.about()
            with open(self.outdir.joinpath("ABOUT.md"), "w") as f:
                f.write(about_str)

    def analyze_all(self):
        self.analyze_unchopped_data()
        self.analyze_chopped_data()
        self.dump_data(dump_all=True)

    @staticmethod
    def about():
        return dedent(
            r"""
            # About

            The data in the folder was generated using the following settings

            This folder contains the following files

            * **original_trace**
                - This is the the raw trace obtained after averaging the frames
                in the stack without any background correction or filtering
            * **corrected_trace**
                - Left panel: This plots the original trace and the backgrund that we
                subtract in the corrected trace.
                - Right panel: The is the corrected version of the original trace where we
                have performed a background correction and filtering.

            * **chopped_data_features**
                - Different panels with chopped data
                - For e.g ADP30, we compute APD30 of all beats plotted in chopped_data:all
                the panel with title all. Then we copmuted the mean and standard deviation (std)
                of all those. Next we exclude those beats that are outside 1 std (default x = 1.0)
                of the mean. This panel shows the beats that are within 1 std of the mean.
            * **chopped_data**
                - Left panel: This if all the beats that we were able to extract from
                the corrected trace
                - Right panel: This is shows the intersection of all beats plotted in chopped_data_z
                described above.
            * **average**
                - These are the average of the traces in chopped_data
            * **data.txt**
                - This contains a short summary of analysis.
            * **data.x where x is either mat of npy**
                - This contains all the output data that can be loaded in python (data.npy)
                of Matlab (data.mat)
            * **unchopped_data.csv**
                - This contains the unchopped traces, i.e the original trace and the corrected
                trace ns a structured formated that are easy to view
            * **chopped_data.csv**
                - This contains the chopped traces, i.e the trace of each beat, the average trace etc,
                in a structured formated that are easy to view
            * **settings.js**
                - Settings used to to perform the analysis. These settings
                can be parsed to the analyze_mps script
            * **metadata.js**
                - Metadata stored within the mps file.

            """,
        )

    def analyze_unchopped_data(self):

        self.dump_data()

        if self.parameters["plot"]:
            plotter.plot_single_trace(
                self.time_stamps,
                self.avg,
                self.outdir.joinpath("original_trace"),
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

        self._chopped_data = apf.chopping.chop_data(
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
            ),
        ):
            self.chopped_data["time_{}".format(i)] = t
            self.chopped_data["trace_{}".format(i)] = d
            self.chopped_data["pacing_{}".format(i)] = p
        self.intervals = self._chopped_data.intervals
        self.upstroke_times = self._chopped_data.upstroke_times

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

        self.features["num_eads"] = ead.analyze_eads(
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
                ),
            ),
        )

    def compute_features(self):
        if not hasattr(self, "_chopped_data"):
            self._set_choopped_data()

        freqs = apf.features.beating_frequency_from_peaks(
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
            ),
        )
        for k in [30, 50, 80, 90]:
            self._all_features[f"capd{k}"] = apf.features.corrected_apd(
                self._all_features[f"apd{k}"],
                60 / self.features["beating_frequency"],
            )

    def compute_mean_features(self):
        if not hasattr(self, "_all_features"):
            self.compute_features()

        self._excluded = exclude_x_std(
            self._all_features,
            self.parameters["std_ex"],
            skips=["upstroke80"],
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
                self._chopped_data.data,
            )
            idx = np.argmax([len(xi) for xi in self._xs_all])
            self.chopped_data["time_all"] = self._xs_all[idx]
        else:
            (
                self.chopped_data["trace_all"],
                self.chopped_data["time_all"],
            ) = average.get_subsignal_average_interpolate(
                self._chopped_data.data,
                self._xs_all,
                N,
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
                self._ys_1std,
            )
            self.chopped_data["time_1std"] = self._xs_1std[
                np.argmax([len(xi) for xi in self._xs_1std])
            ]
        else:
            (
                self.chopped_data["trace_1std"],
                self.chopped_data["time_1std"],
            ) = average.get_subsignal_average_interpolate(
                self._ys_1std,
                self._xs_1std,
                N,
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

            bkg = apf.background.correct_background(self.time_stamps, self.avg, "full")
            background_ = bkg.background
            bkgplot_y = np.transpose([np.copy(self.avg), background_])
            bkgplot_x = np.transpose([self.time_stamps, self.time_stamps])
            self.avg = bkg.corrected

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
            self.time_stamps, self.avg = apf.filters.remove_points(
                self.time_stamps, self.avg, *p
            )
            _, self.pacing = apf.filters.remove_points(t_, self.pacing, *p)

    def ignore_pacing(self):
        if self.parameters["ignore_pacing"]:
            logger.debug("Ignore pacing")
            self.pacing = np.multiply(self.pacing, 0.0)

    def filter(self):
        if self.parameters["filter"]:
            logger.debug("Filter signal")
            self.avg = apf.filters.filt(self.avg)

    def remove_spikes(self):
        if self.parameters["spike_duration"] == 0:
            return

        sd = self.parameters["spike_duration"]
        logger.debug(f"Remove spikes with spike duration {sd}")
        self.avg = apf.filters.remove_spikes(self.avg, self.pacing, sd)
        self.time_stamps = apf.filters.remove_spikes(self.time_stamps, self.pacing, sd)
        self.pacing = apf.filters.remove_spikes(self.pacing, self.pacing, sd)

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
            filter=False,
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


def frame2average(frame, times=None, normalize=False, background_correction=True):
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
        avg = apf.background.correct_background(times, avg_, "full").corrected
    else:
        avg = avg_

    if normalize:
        trace = apf.utils.normalize_signal(avg)
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
        N,
        x_start=x_start,
        y_start=y_start,
        x_end=x_end,
        y_end=y_end,
        frames=frames,
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

    return np.fliplr(local_averages)


def _frames2average(kwargs):
    return frame2average(**kwargs)


def analyze_local_average(frames, times=None, mask=None, N=10):

    loc = local_averages(frames, times, N=N)
    avg = frame2average(
        frame=frames,
        times=times,
        normalize=False,
        background_correction=True,
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
            baseline_intensity[i, j, :] = apf.background.background(
                times,
                loc_ij,
                "full",
            ).corrected

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
    chopped_data = apf.chopping.chop_data_without_pacing(avg, mps_data.time_stamps)
    global_freq = apf.features.beating_frequency(
        chopped_data.times,
        unit=mps_data.info.get("time_unit", "ms"),
    )
    logger.info(f"Global frequency: {global_freq}")
    global_snr = snr(
        apf.background.correct_background(mps_data.time_stamps, avg, "full").corrected,
    )
    logger.info(f"Global SNR: {global_snr}")

    # Find the beating frequency for each local average
    shape = loc.shape[:2]
    freq = np.zeros(shape)
    local_snr = np.zeros(shape)
    is_tissue = np.ones(shape, dtype=bool)
    for i, j in it.product(range(shape[0]), range(shape[1])):
        loc_ij = loc[i, j, :]
        chopped_data = apf.chopping.chop_data_without_pacing(
            loc_ij,
            mps_data.time_stamps,
        )
        freq[i, j] = apf.features.beating_frequency(
            chopped_data.times,
            unit=mps_data.info.get("time_unit", "ms"),
        )
        local_snr[i, j] = snr(loc[i, j, :])

    # Set non-beating regions where the local signal to
    # noise ratio is high
    is_beating = local_snr <= snr_factor * global_snr

    # Frequencies that deviates too much from the
    # global frequency is most likely noise
    freq_dev = np.where(
        ~np.isclose(freq, global_freq, frequency_threshold * global_freq),
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
