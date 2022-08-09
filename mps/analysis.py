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
import operator as op
import json
import logging
from collections import namedtuple, defaultdict
from collections.abc import Iterable
from pathlib import Path
from textwrap import dedent
from typing import Any, Dict, List, Optional

import ap_features as apf
import numpy as np

from . import average, plotter, utils

logger = logging.getLogger(__name__)

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

Features = namedtuple(
    "Features",
    [
        "features",  # Mean features
        "features_beats",  # Features for all beats
        "features_included_beats",  # Features for the included beats
        "features_1std_str",  # Same as features but as strings
        "features_all_str",  # Sane as 'features_1std_str' but for all beats
        "included_indices",  # List of common indices included
        "included_indices_beats",  # Included indices for all beats
    ],
)

ExcludedData = namedtuple(
    "ExcludedData",
    ("new_data, included_indices, " "all_included_indices"),
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
    beats: List[apf.Beat],
    max_allowed_apd_change: Optional[float] = None,
    fname: str = "",
    plot=True,
) -> APDAnalysis:

    apd_levels = (100 * np.sort(np.append(np.arange(0.1, 0.91, 0.2), 0.8))).astype(int)
    apds = {k: [b.apd(k) for b in beats] for k in apd_levels}
    apd_points = {k: [b.apd_point(k) for b in beats] for k in apd_levels}
    apd_dt = {k: [v[0] for v in value] for k, value in apd_points.items()}
    capds = {k: [b.capd(k) for b in beats] for k in apd_levels}

    slope_APD80, const_APD80 = apf.beat.apd_slope(beats=beats, factor=80)
    slope_cAPD80, const_cAPD80 = apf.beat.apd_slope(
        beats=beats,
        factor=80,
        corrected_apd=True,
    )

    triangulation = [b.triangulation(low=30, high=80) for b in beats]

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
    logger.debug(msg)

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
            beats=beats,
            res=res,
            fname=fname,
        )

    return res


def plot_apd_analysis(beats: List[apf.Beat], res: APDAnalysis, fname=""):
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

    for beat in beats:
        ax[2].plot(beat.t, beat.y, color="k")
    for i, (label, vals) in enumerate(res.apd_points.items()):
        x = [v[0] for v in vals] + [v[1] for v in vals]  # type: ignore

        y = [beats[j].as_spline(k=3, s=0)(v[0]) for j, v in enumerate(vals)] + [
            beats[j].as_spline(k=3, s=0)(v[1]) for j, v in enumerate(vals)
        ]
        ax[2].plot(x, y, linestyle="", marker="o", label=label)

    ax[2].legend()
    ax[2].grid()
    ax[2].set_xlabel("Time [ms]")
    ax[2].set_ylabel(r"$\Delta F / F$")

    if fname != "":
        fig.savefig(fname)
    plt.close()


def compute_features(beats: List[apf.Beat], use_spline=True, normalize=False):
    r"""
    Analyze signals. Compute all features and
    include only the relevant beats


    Arguments
    ---------
    beats : List[apf.Beat]
        List of beats
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
    num_beats = len(beats)
    if num_beats == 0:
        return {"num_beats": 0}

    assert beats[0].parent is not None, "Beats must have a parent"

    features: Dict[str, Any] = dict(
        apd30=np.zeros(num_beats),
        apd50=np.zeros(num_beats),
        apd80=np.zeros(num_beats),
        apd90=np.zeros(num_beats),
        capd30=np.zeros(num_beats),
        capd50=np.zeros(num_beats),
        capd80=np.zeros(num_beats),
        capd90=np.zeros(num_beats),
        dFdt_max=np.zeros(num_beats),
        int30=np.zeros(num_beats),
        tau75=np.zeros(num_beats),
        upstroke80=np.zeros(num_beats),
        ttp=np.zeros(num_beats),
    )

    for i, beat in enumerate(beats):
        beat.ensure_time_unit("ms")

        features["apd90"][i] = beat.apd(90, use_spline=use_spline) or np.nan
        features["apd80"][i] = beat.apd(80, use_spline=use_spline) or np.nan
        features["apd50"][i] = beat.apd(50, use_spline=use_spline) or np.nan
        features["apd30"][i] = beat.apd(30, use_spline=use_spline) or np.nan
        features["capd90"][i] = beat.capd(90, use_spline=use_spline) or np.nan
        features["capd80"][i] = beat.capd(80, use_spline=use_spline) or np.nan
        features["capd50"][i] = beat.capd(50, use_spline=use_spline) or np.nan
        features["capd30"][i] = beat.capd(30, use_spline=use_spline) or np.nan
        features["ttp"][i] = beat.ttp() or np.nan
        features["tau75"][i] = beat.tau(0.75) or np.nan
        features["upstroke80"][i] = beat.upstroke(0.8) or np.nan
        features["dFdt_max"][i] = (
            beat.maximum_upstroke_velocity(
                use_spline=use_spline,
                normalize=normalize,
            )
            or np.nan
        )
        features["int30"][i] = (
            beat.integrate_apd(
                0.3,
                use_spline=use_spline,
                normalize=normalize,
            )
            or np.nan
        )
    features["beating_frequencies"] = beats[0].parent.beating_frequencies
    for k, v in features.items():
        features[k] = v[~np.isnan(v)]
        num_beats = min(len(features[k]), num_beats)

    features["beating_frequency"] = beats[0].parent.beating_frequency
    features["num_beats"] = num_beats
    return features


def find_included_indices(data, x=None, use=None):
    """Given a list of values return a list of all
    values that are within x factor of the mean.

    Parameters
    ----------
    data : dict
        A dictionary of lists of values, e.g a list of apd30
    x : float
        The number of standard deviations to be
        included, ex x = 1.0. If none is provided
        everything will be included
    use : List[str]
        List of features to use

    Returns
    -------
    Tuple[Dict[str, List[int]], List[int]]
        List of indices or each beat and list of common indices

    """
    if use is None:
        use = ["apd30", "apd80"]

    included_indices = defaultdict(list)  # {k: [] for k in data.keys()}

    # Find minimum number of beats:
    num_beats = 0
    for k, v in data.items():
        if not isinstance(v, Iterable):
            continue
        if num_beats == 0:
            num_beats = len(v)
        num_beats = min(num_beats, len(v))

    for k, v in data.items():
        if not isinstance(v, Iterable) or k not in use:
            continue
        if len(v) == 0:
            continue
        for j, s in enumerate(v):

            # Include only signals within factor *
            # standard deviation from the mean
            if -np.std(v) * x < (s - np.mean(v)) < np.std(v) * x:
                included_indices[k].append(j)

        # Check if list is empty
        if not len(included_indices[k]) == 0:
            included_indices[k] = list(range(num_beats))

    intsect = utils.get_intersection(included_indices)
    return dict(included_indices), [ind for ind in intsect if ind < num_beats]


def exclude_x_std(data, x=None, use=None):
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
    use : List[str]
        List of features to use

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
    included_indices, all_included_indices = find_included_indices(data, x, use)
    new_data = {}
    all_included_indices = np.array(all_included_indices)
    if len(all_included_indices) == 0:
        all_included_indices = np.arange(data["num_beats"])

    for k, v in data.items():
        if isinstance(v, Iterable) and len(v) > 0:
            new_data[k] = np.array(v)[all_included_indices].tolist()
        else:
            new_data[k] = v

    return ExcludedData(
        new_data=new_data,
        included_indices=included_indices,
        all_included_indices=all_included_indices,
    )


def analyze_frequencies(
    beats,
    time_unit: str = "ms",
    fname: str = "",
    plot=True,
) -> np.ndarray:

    freqs = apf.features.beating_frequency_from_peaks(
        [beat.y for beat in beats],
        [beat.t for beat in beats],
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

        points_x = [beat.t[int(np.argmax(beat.y))] for beat in beats]
        points_y = [np.max(beat.y) for beat in beats]

        for beat in beats:
            ax[1].plot(beat.t, beat.y)

        ax[1].plot(points_x, points_y, linestyle="", marker="o", color="r")
        ax[1].set_xlabel("Time [ms]")
        ax[1].set_ylabel(r"$\Delta F / F$")
        if fname != "":
            fig.savefig(fname)
        plt.close()

    return freqs


def analyze_eads(
    beats: List[apf.Beat],
    sigma: float = 1,
    prominence_threshold: float = 0.07,
    plot=True,
    fname: str = "",
) -> int:
    """
    Loop over all beats and check for EADs.

    Arguments
    ---------
    beats : List[apf.Beat]
        List of beats
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
    for i, beat in enumerate(beats):
        has_ead, peak_index = beat.detect_ead(
            sigma=sigma,
            prominence_level=prominence_threshold,
        )
        if has_ead and peak_index is not None:
            num_eads += 1
            peak_inds[i] = peak_index

    if plot:
        plot_ead(fname, beats, peak_inds)

    return num_eads


def plot_ead(fname, beats, peak_inds):

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return

    fig, ax = plt.subplots()
    num_eads = 0
    for i, beat in enumerate(beats):
        ax.plot(beat.t, beat.y)
        if i in peak_inds:
            ax.plot([beat.t[peak_inds[i]]], [beat.y[peak_inds[i]]], "ro")
    ax.set_title(f"EADs are marked with red dots. Found EADs in {num_eads} beats")
    if fname != "":
        fig.savefig(fname)
    plt.close()


class Collector:
    def __init__(self, outdir=None, plot: bool = False, params=None):

        self.unchopped_data: Dict[str, Any] = {}
        self.chopped_data: Dict[str, Any] = {}
        self.intervals: List[apf.chopping.Intervals] = []
        self.features = Features({}, {}, {}, [], [], [], {})
        self._all_features: Dict[str, Any] = {}
        self.upstroke_times: List[float] = []
        self.chopping_parameters: Dict[str, Any] = {}
        self.params = {} if params is None else params

        if outdir is not None:
            outdir = Path(outdir)
            logger.info(f"Create directory {outdir}")
            outdir.mkdir(exist_ok=True, parents=True)
        self.outdir = outdir
        self.plot = plot

    @property
    def info(self) -> Dict[str, Any]:
        return self._info

    @info.setter
    def info(self, info: Optional[Dict[str, Any]]) -> None:
        if info is None:
            info = {}
        self._info = info

    @property
    def metadata(self) -> Dict[str, Any]:
        return self._metadata

    @metadata.setter
    def metadata(self, metadata: Optional[Dict[str, Any]]) -> None:
        if metadata is None:
            metadata = {}
        self._metadata = metadata

    def register_unchopped(self, trace: apf.Beats, label: str = "") -> None:
        if label != "":
            label += "_"

        self.y = np.copy(trace.y)
        self.t = np.copy(trace.t)
        self.pacing = np.copy(trace.pacing)
        self.unchopped_data[label + "trace"] = self.y
        self.unchopped_data[label + "times"] = self.t
        self.unchopped_data[label + "pacing"] = self.pacing
        self.dump()
        if self.plot:
            fname = (
                self.outdir.joinpath("trace")
                if label == ""
                else self.outdir.joinpath(label.strip("_"))
            )
            trace.plot(fname=fname.with_suffix(".png"))

    def register_chopped_data(
        self,
        chopped_data: apf.chopping.ChoppedData,
        aligned_beats=None,
    ) -> None:

        if len(chopped_data.data) == 0 and hasattr(self, "y"):
            chopped_data.data.append(np.copy(self.y))
            chopped_data.times.append(np.copy(self.t))
            chopped_data.pacing.append(np.copy(self.pacing))

        for i, (t, d, p) in enumerate(
            zip(
                chopped_data.times,
                chopped_data.data,
                chopped_data.pacing,
            ),
        ):
            if aligned_beats is not None and len(aligned_beats) > i:
                beat = aligned_beats[i]
                alinged_t = beat.t
                alinged_d = beat.y
            else:
                alinged_t = t
                alinged_d = d

            self.chopped_data[f"time_{i}"] = t
            self.chopped_data[f"trace_{i}"] = d
            self.chopped_data[f"pacing_{i}"] = p
            self.chopped_data[f"aligned_time_{i}"] = alinged_t
            self.chopped_data[f"aligned_trace_{i}"] = alinged_d
        self.intervals = chopped_data.intervals
        self.upstroke_times = chopped_data.upstroke_times

    @property
    def data(self):
        return {
            "chopped_data": self.chopped_data,
            "unchopped_data": self.unchopped_data,
            "features": self.features.features,
            "all_features": self.features._asdict(),
            "chopping_parameters": self.chopping_parameters,
            "attributes": self.info,
            "intervals": self.intervals,
            "upstroke_times": self.upstroke_times,
            "settings": self.params,
        }

    @property
    def info_txt(self) -> str:
        template = "{0:20}\t{1:>20}"
        info_txt = [template.format(k, v) for k, v in self.info.items()]

        return "\n".join(
            np.concatenate(
                (
                    ["All features"],
                    self.features.features_all_str,
                    ["\n"],
                    ["Features within 1std"],
                    [
                        template.format(
                            "Included beats:",
                            ",".join(map(str, self.features.included_indices)),
                        ),
                    ],
                    self.features.features_1std_str,
                    ["\n"],
                    ["Info"],
                    info_txt,
                    ["\nSettings"],
                    [template.format(k, str(v)) for k, v in self.params.items()],
                ),
            ),
        )

    def dump(self):
        if self.outdir is None:
            return
        np.save(self.outdir.joinpath("data"), self.data)

    def dump_all(self):
        if self.outdir is None:
            return
        self.dump()
        chopped_data_padded = utils.padding(self.chopped_data)
        unchopped_data_padded = utils.padding(self.unchopped_data)
        utils.to_csv(chopped_data_padded, self.outdir.joinpath("chopped_data"))
        utils.to_csv(unchopped_data_padded, self.outdir.joinpath("unchopped_data"))

        with open(self.outdir.joinpath("data.txt"), "w") as f:
            f.write(self.info_txt)

        # utils.dump_data(data, self.outdir.joinpath("data"), "mat")
        with open(self.outdir.joinpath("metadata.json"), "w") as f:
            json.dump(self.metadata, f, indent=4, default=utils.json_serial)

        about_str = about()
        with open(self.outdir.joinpath("README.txt"), "w") as f:
            f.write(about_str)


def about():
    import datetime
    from . import __version__
    import scipy
    import sys

    sys_version = sys.version.replace("\n", " ")
    return dedent(
        f"""
        Time
        ----
        These results were generated on
        {datetime.datetime.now()} (local time)
        {datetime.datetime.utcnow()} (UTC)

        Setup
        -----
        The following packages
        mps version {__version__}
        ap_features version {apf.__version__}
        numpy version {np.__version__}
        scipy version {scipy.__version__}

        and python version
        {sys_version}
        and platform {sys.platform}

        Content
        -------
        This folder contains the following files

        - apd_analysis.png
            Figure showing plots of action potential durations
            and corrected actions potential durations (using
            the friderica correction formula). Top panel shows
            bars of the APD for the different beats where the *
            indicated if the APD for that beat is significantly different.
            The middle panel show the APD for each each plotted as a
            line as well as a linear fit of the APD80 and cAPD80 line.
            The intention behind this panel is to see it there is any
            correlation between the APD and the beat number
            The lower panel shows where the cut is made to compute
            the different APDs
        - average_pacing.png
            These are the average trace with pacing
        - average.png
            These are the average of the traces in chopped_data
        - background.png
            This plots the original trace and the background that we
            subtract in the corrected trace.
        - chopped_data_aligned.png
            All beats plotted on with the time starting at zero
        - chopped_data.csv
            A csv file with the chopped data
        - chopped_data.png
            Left panel: All the beats that we were able to extract from the corrected trace
            Right panel: The intersection of all beats where APD80 and APD30 are within 1
            standard deviation of the mean.
        - data.npy
            A file containing all the results. You can load the results
            in python as follows
            >> import numpy as np
            >> data = np.load('data.npy', allow_pickle=True).item()
            data is now a python dictionary.
        - data.txt
            This contains a short summary of analysis.
        - EAD_analysis.png
            A plot of the beats containing EADs
        - frequency_analysis.png
            A plot of the frequencies of the different beats computed
            as the time between peaks.
        - metadata.json
            Metadata stored inside the raw data
        - original_pacing.png
            This is the the raw trace obtained after averaging the frames
            in the stack without any background correction or filtering
            together with the pacing amplitude.
        - original.png
            This is the the raw trace obtained after averaging the frames
            in the stack without any background correction or filtering
        - sliced_filtered.png
            Original trace where we have performed slicing and filtering
        - trace.png
            The is the corrected version of the original trace where we
            have performed a background correction and filtering
        - unchopped_data.csv
            A csv file with all the unchopped data (original and corrected)
        """,
    )


def analyze_unchopped_data(
    mps_data,
    collector,
    mask=None,
    analysis_window_start=0,
    analysis_window_end=-1,
    spike_duration=0,
    filter_signal=False,
    ignore_pacing=False,
    remove_points_list=(),
    background_correction=True,
    **kwargs,
) -> apf.Beats:

    avg = average_intensity(mps_data.frames, mask=mask)
    time_stamps = mps_data.time_stamps
    pacing = mps_data.pacing

    collector.info = mps_data.info
    collector.metadata = mps_data.metadata

    trace = apf.Beats(
        y=avg,
        t=time_stamps,
        pacing=pacing,
    )
    # Original data
    collector.register_unchopped(trace, "original")
    if collector.plot and collector.outdir is not None:
        trace.plot(
            collector.outdir.joinpath("original_pacing.png"),
            include_pacing=True,
        )
    # Remove spikes and filter
    filtered_trace = trace.remove_spikes(spike_duration=spike_duration)

    if filter_signal:
        logger.info("Filter signal")
        filtered_trace = trace.filter()
    # Sliced data
    sliced_filt_trace = filtered_trace.slice(analysis_window_start, analysis_window_end)
    collector.register_unchopped(sliced_filt_trace, "sliced_filtered")

    if ignore_pacing:
        sliced_filt_trace.pacing[:] = 0

    for p in remove_points_list:
        sliced_filt_trace = sliced_filt_trace.remove_points(*p)

    # Corrected data
    background_correction_method = "full" if background_correction else "none"
    corrected_trace: apf.Beats = sliced_filt_trace.correct_background(
        background_correction_method,
    )
    collector.register_unchopped(corrected_trace, "")
    collector.unchopped_data["background"] = np.copy(corrected_trace.background)
    if collector.plot and collector.outdir is not None:
        corrected_trace.plot(
            collector.outdir.joinpath("background.png"),
            include_background=True,
        )
    return corrected_trace


def compute_all_features(
    beats,
    outdir=None,
    plot=True,
    use_spline: bool = True,
    normalize: bool = False,
    max_allowed_apd_change: Optional[float] = None,
    ead_sigma: int = 3,
    ead_prominence_threshold: float = 0.04,
    std_ex: float = 1.0,
    **kwargs,
):
    features = compute_features(beats, use_spline=use_spline, normalize=normalize)

    apd_analysis = analyze_apds(
        beats=beats,
        max_allowed_apd_change=max_allowed_apd_change,
        plot=plot,
        fname="" if outdir is None else outdir.joinpath("apd_analysis.png"),
    )
    features["slope_APD80"] = apd_analysis.slope_APD80
    features["slope_cAPD80"] = apd_analysis.slope_cAPD80
    features["triangulation"] = apd_analysis.triangulation
    features["num_eads"] = analyze_eads(
        beats=beats,
        sigma=ead_sigma,
        prominence_threshold=ead_prominence_threshold,
        fname="" if outdir is None else outdir.joinpath("EAD_analysis.png"),
    )
    excluded = exclude_x_std(
        features,
        std_ex,
        use=["apd30", "apd80"],
    )
    mean_features = {}
    features_1std = []
    features_all = []
    template = "{0:20}\t{1:10.1f}  +/-  {2:10.3f}"
    for k in features.keys():
        v = excluded.new_data[k]
        features_1std.append(template.format(k, mean(v), std(v)))
        mean_features[k] = mean(v)

        v = features[k]
        features_all.append(template.format(k, mean(v), std(v)))

    return Features(
        features=mean_features,
        features_beats=features,
        features_included_beats=excluded.new_data,
        features_1std_str=features_1std,
        features_all_str=features_all,
        included_indices=excluded.all_included_indices,
        included_indices_beats=excluded.included_indices,
    )


def mean(x):
    if np.isscalar(x):
        return x
    if len(x) == 0:
        return np.nan
    return np.mean(x)


def std(x):
    if np.isscalar(x):
        return x
    if len(x) == 0:
        return np.nan
    return np.std(x)


def enlist(x):
    if not isinstance(x, (list, tuple)):
        return [x]
    return x


def analyze_chopped_data(
    trace: apf.Beats,
    collector: Collector,
    threshold_factor=0.3,
    extend_front=None,
    extend_end=None,
    min_window=200,
    max_window=2000,
    use_spline: bool = True,
    normalize: bool = False,
    max_allowed_apd_change: Optional[float] = None,
    ead_sigma: int = 3,
    ead_prominence_threshold: float = 0.04,
    std_ex: float = 1.0,
    N: int = 200,
    **kwargs,
):

    chopping_parameters = dict(
        threshold_factor=threshold_factor,
        extend_front=extend_front,
        extend_end=extend_end,
        min_window=min_window,
        max_window=max_window,
    )

    trace.chopping_options.update(chopping_parameters)
    chopped_data = trace.chopped_data

    beats = apf.beat.chopped_data_to_beats(chopped_data, parent=trace)
    aligned_beats = apf.beat.align_beats(beats, parent=trace)
    collector.register_chopped_data(chopped_data, aligned_beats=aligned_beats)

    if collector.plot and collector.outdir is not None:
        apf.plot.plot_beats(beats, fname=collector.outdir.joinpath("chopped_data.png"))
        apf.plot.plot_beats(
            beats,
            fname=collector.outdir.joinpath("chopped_data_aligned.png"),
            align=True,
        )
    fname = ""
    if collector.plot and collector.outdir is not None:
        collector.outdir.joinpath("frequency_analysis.png")
    analyze_frequencies(
        beats,
        fname=fname,
        plot=collector.plot,
    )
    collector.features = compute_all_features(
        beats,
        outdir=collector.outdir,
        plot=collector.plot,
        use_spline=use_spline,
        normalize=normalize,
        max_allowed_apd_change=max_allowed_apd_change,
        ead_sigma=ead_sigma,
        ead_prominence_threshold=ead_prominence_threshold,
        std_ex=std_ex,
    )
    collector.dump()

    if len(collector.features.included_indices) > 0:
        inds = op.itemgetter(*collector.features.included_indices)
        included_beats = enlist(inds(beats))
    else:
        included_beats = beats
    average_all = apf.beat.average_beat(beats, N=N)
    average_1std = apf.beat.average_beat(included_beats, N=N)
    collector.chopped_data["time_1std"] = average_1std.t
    collector.chopped_data["trace_1std"] = average_1std.y
    collector.chopped_data["pacing_1std"] = average_1std.pacing
    collector.chopped_data["time_all"] = average_all.t
    collector.chopped_data["trace_all"] = average_all.y
    collector.chopped_data["pacing_all"] = average_all.pacing

    if collector.plot and collector.outdir is not None:
        plotter.plot_multiple_traces(
            [[b.t for b in beats], [b.t for b in included_beats]],
            [[b.y for b in beats], [b.y for b in included_beats]],
            collector.outdir.joinpath("chopped_data.png"),
            ["all", "1 std"],
            ylabels=[r"$\Delta F / F$" for _ in range(len(beats))],
            deep=True,
        )

        plotter.plot_multiple_traces(
            [collector.chopped_data["time_all"], collector.chopped_data["time_1std"]],
            [collector.chopped_data["trace_all"], collector.chopped_data["trace_1std"]],
            collector.outdir.joinpath("average.png"),
            ["all", "1 std"],
            ylabels=[r"$\Delta F / F$" for _ in range(len(beats))],
        )

        plotter.plot_twin_trace(
            collector.chopped_data["time_1std"],
            collector.chopped_data["time_1std"],
            collector.chopped_data["trace_1std"],
            collector.chopped_data["pacing_1std"],
            collector.outdir.joinpath("average_pacing.png"),
        )


def analyze_mps_func(
    mps_data,
    mask=None,
    analysis_window_start=0,
    analysis_window_end=-1,
    spike_duration=0,
    filter_signal=False,
    ignore_pacing=False,
    remove_points_list=(),
    threshold_factor=0.3,
    extend_front=None,
    extend_end=None,
    min_window=200,
    max_window=2000,
    use_spline=True,
    normalize=False,
    outdir=None,
    plot=False,
    background_correction=True,
    max_allowed_apd_change=None,
    ead_sigma=3,
    ead_prom=0.04,
    std_ex: float = 1.0,
    N=200,
    **kwargs,
):
    logger.info("Analyze MPS data")
    params = dict(
        analysis_window_start=analysis_window_start,
        analysis_window_end=analysis_window_end,
        spike_duration=spike_duration,
        filter_signal=filter_signal,
        ignore_pacing=ignore_pacing,
        remove_points_list=remove_points_list,
        background_correction=background_correction,
        use_spline=use_spline,
        normalize=normalize,
        max_allowed_apd_change=max_allowed_apd_change,
        ead_sigma=ead_sigma,
        ead_prom=ead_prom,
        std_ex=std_ex,
        N=N,
        threshold_factor=threshold_factor,
        extend_front=extend_front,
        extend_end=extend_end,
        min_window=min_window,
        max_window=max_window,
    )

    collector = Collector(outdir=outdir, plot=plot, params=params)

    corrected_trace = analyze_unchopped_data(
        mps_data=mps_data, collector=collector, mask=mask, **params
    )
    try:
        analyze_chopped_data(trace=corrected_trace, collector=collector, **params)
    except Exception:
        logger.warning("Failed analyzing chopped data")

    collector.dump_all()

    try:
        import matplotlib.pyplot as plt

        plt.close("all")
    finally:
        return collector.data


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
