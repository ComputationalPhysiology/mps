#!/usr/bin/env python3
# c) 2001-2019 Simula Research Laboratory ALL RIGHTS RESERVED
#
# END-USER LICENSE AGREEMENT
# PLEASE READ THIS DOCUMENT CAREFULLY. By installing or using this
# software you agree with the terms and conditions of this license
# agreement. If you do not accept the terms of this license agreement
# you may not install or use this software.

# Permission to use, copy, modify and distribute any part of this
# software for non-profit educational and research purposes, without
# fee, and without a written agreement is hereby granted, provided
# that the above copyright notice, and this license agreement in its
# entirety appear in all copies. Those desiring to use this software
# for commercial purposes should contact Simula Research Laboratory AS:
# post@simula.no
#
# IN NO EVENT SHALL SIMULA RESEARCH LABORATORY BE LIABLE TO ANY PARTY
# FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
# INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE
# "MPS" EVEN IF SIMULA RESEARCH LABORATORY HAS BEEN ADVISED
# OF THE POSSIBILITY OF SUCH DAMAGE. THE SOFTWARE PROVIDED HEREIN IS
# ON AN "AS IS" BASIS, AND SIMULA RESEARCH LABORATORY HAS NO OBLIGATION
# TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
# SIMULA RESEARCH LABORATORY MAKES NO REPRESENTATIONS AND EXTENDS NO
# WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESSED, INCLUDING, BUT
# NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS
__author__ = "Henrik Finsberg (henriknf@simula.no), 2017--2019"
__maintainer__ = "Henrik Finsberg"
__email__ = "henriknf@simula.no"
from collections import OrderedDict
import json
from timeit import default_timer as timer
from typing import List
import datetime
from textwrap import dedent
import os.path
import argparse
import sys
import logging
import shutil
import numpy as np
import scipy.io

import mps

from .plotter import plt


def version_parser(parser):
    parser.add_argument(
        "--version",
        action="version",
        version="MPS version {version}".format(version=get_version()),
    )


brightfield_patterns = ["BF", "bf", "brightfield", "Brightfield"]


class Script:

    valid_extensions = [".nd2", ".czi", ".stk"]
    skip_startswith = ["."]
    skip_endswith = ["~", "#"]
    exclude_patterns: List[str] = []
    logger = mps.utils.get_logger("base")
    _make_outdir = False

    def __init__(self, args):
        self.args = args
        self.check_args()

    def check_args(self):
        pass

    def _check_args(self, fname):

        _, ext = os.path.splitext(fname)
        if self.valid_extensions and ext not in self.valid_extensions:
            logging.debug(
                f"Invalid file extension {ext}, expected on of {self.valid_extensions}"
            )
            return False

        for pattern in self.exclude_patterns:
            if pattern in fname:
                return False

        for sw in self.skip_startswith:
            if os.path.basename(fname).startswith(sw):
                logging.debug(f"Skip {fname}")
                return False

        for ew in self.skip_endswith:
            if os.path.basename(fname).endswith(ew):
                logging.debug(f"Skip {fname}")
                return False

        assert os.path.isfile(fname), f"File {fname} does not exist"
        return True

    def main_func(self, func):

        log_level = logging.DEBUG if self.args.get("verbose", False) else logging.INFO
        self.logger.setLevel(log_level)
        for h in self.logger.handlers:
            h.setLevel(log_level)

        mps.set_log_level(log_level)
        settings = "\n    Settings used:\n" + "\n".join(
            ["\t{}\t{}".format(k, v) for k, v in sorted(self.args.items())]
        )
        self.logger.debug(settings)

        outdir = None
        if os.path.isfile(self.args["file"]):
            if self._check_args(self.args["file"]):
                if self._make_outdir:
                    outdir = make_outdir(self.args)
                func(self.args["file"], outdir)
        elif os.path.isdir(self.args["file"]):
            for root, _, files in os.walk(self.args["file"]):
                for f in files:
                    path = os.path.join(root, f)
                    if self._check_args(path):
                        if self._make_outdir:
                            outdir = make_outdir({"outdir": "", "file": path})
                        try:
                            func(path, outdir)
                        except Exception as ex:
                            self.logger.warning(ex, exc_info=True)
        else:
            print(f"Invalid input self.args['file']")

    @staticmethod
    def get_parser():
        raise NotImplementedError

    @classmethod
    def run(cls, sysargs=None, gooeyargs=None):

        if gooeyargs is not None:
            args = gooeyargs
        else:
            parser = cls.get_parser()
            version_parser(parser)

            if sysargs is not None and len(sysargs) == 0:
                parser.print_help()
                sys.exit(1)

            if len(sys.argv) == 1:
                parser.print_help()
                sys.exit(1)

            args = vars(parser.parse_args(args=sysargs))

        c = cls(args)
        c.main_func()


def make_outdir(args):

    if args["outdir"] == "":
        outdir = os.path.splitext(args["file"])[0]
    else:
        outdir = args["outdir"]

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    return outdir


def get_version():
    return mps.__version__


def padding(data, fill_value=0):
    """
    Make sure all traces are of the same lenght and
    fill in extra values if neccessary
    """
    N = np.max([len(a) for a in data.values()])
    padded_data = OrderedDict()
    for k, v in data.items():
        padded_data[k] = np.concatenate((v, fill_value * np.ones(N - len(v))))

    return padded_data


def json_serial(obj):
    """JSON serializer for objects not serializable by default json code"""

    if isinstance(obj, (datetime.datetime)):
        return obj.isoformat()
    elif isinstance(obj, (np.ndarray)):
        return obj.tolist()
    raise TypeError("Type %s not serializable" % type(obj))


def dump_settings(outdir, args):
    from . import __version__

    args["time"] = datetime.datetime.now()
    args["mps_version"] = __version__
    args["full_path"] = os.path.abspath(args["file"])
    with open(os.path.join(outdir, "settings.json"), "w") as f:
        json.dump(args, f, indent=4, default=json_serial)


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

    """
    )


class analyze_mps(Script):

    logger = mps.utils.get_logger("analyze_mps")
    _make_outdir = True
    exclude_patterns = brightfield_patterns

    def main_func(self):
        def func(path, outdir):
            averaging_parameters = dict(
                alpha=self.args["alpha"],
                averaging_type=self.args["averaging_type"],
                local=self.args["local"],
            )

            start = timer()
            args = self._local_args
            args["outdir"] = outdir
            args["file"] = path
            mps_data = mps.MPS(path, averaging_parameters)
            mps.analysis.analyze_mps_func(mps_data, plot=True, **args)
            dump_settings(outdir, args)
            end = timer()
            self.logger.info(
                (
                    "Finished analyzing MPS data. Data stored in {}. ".format(outdir)
                    + "\nTotal elapsed time: {} seconds".format(end - start)
                )
            )

        super().main_func(func)

    @staticmethod
    def get_parser():

        descr = "Convert czi and nd2 files containing MPS data into averages"
        usage = "Example: \nanalyze_mps voltage.czi"

        parser = argparse.ArgumentParser(description=descr, usage=usage, add_help=True)

        parser.add_argument(
            action="store",
            dest="file",
            type=str,
            help="Path to the .*czi or .*nd2 file",
        )

        parser.add_argument(
            "-o",
            "--outdir",
            action="store",
            dest="outdir",
            default="",
            type=str,
            help=(
                "Output directory for where you want to store the output. "
                "If not provided a folder with the same name as "
                "the basename of the input file in the currect directory will "
                "be used. Example: `analyze_mps voltage.czi -o figures` will store "
                "the output in the folder 'figures'."
            ),
        )

        return _analyze_mps_args(parser)

    def _check_args(self, fname):

        if not super()._check_args(fname):
            return False

        args = self.args.copy()
        outdir = make_outdir({"outdir": "", "file": fname})
        args["outdir"] = outdir
        args["file"] = fname

        if args["reuse_settings"] and os.path.isfile(
            os.path.join(outdir, "settings.json")
        ):
            self.logger.debug("Reuse settings")
            with open(os.path.join(outdir, "settings.json"), "r") as f:
                old_args = json.load(f)

            args.update(**old_args)

        if os.path.exists(outdir):
            if args["overwrite"]:
                try:
                    shutil.rmtree(outdir)
                except Exception as ex:
                    self.logger.warning(ex, exc_info=True)

            else:
                # Check if we can find the version number of the
                # software used to analyze the data
                try:
                    with open(os.path.join(outdir, "settings.json"), "r") as f:
                        settings = json.load(f)
                    version = settings["mps_version"]
                except (FileNotFoundError, KeyError):
                    old_dir = os.path.join(outdir, "old")
                else:
                    old_dir = os.path.join(outdir, "old_version_{}".format(version))

                if not os.path.exists(old_dir):
                    os.makedirs(old_dir)

                    for f in os.listdir(outdir):

                        if f == os.path.basename(old_dir):
                            continue
                        if f.startswith("._"):
                            try:
                                os.remove(os.path.join(outdir, f))
                            except FileNotFoundError:
                                pass

                            continue
                        shutil.move(os.path.join(outdir, f), os.path.join(old_dir, f))
        self._local_args = args
        return True


def _analyze_mps_args(parser):

    parser.add_argument(
        "-d",
        "--delete",
        action="append",
        dest="remove_points_list",
        default=[],
        type=int,
        nargs=2,
        help=(
            "Add the timepoints that should be removed. This is useful if the "
            "know that part of the signal is really bad, and could potentially corrupt "
            "the whole average. You need to provide a start and a stop."
            "\nExample :\n`./analyze_mps my_file.czi -d 100 120`\n will remove the "
            "points with times between 100 and 120. If you want to remove even more points "
            "then you can simple do ./analyze_mps my_file.czi -d 100 120 -d 110 120. "
            "This will remove the points 100 to 120 (i.e 20 points) and 130 (=110+20) to 140."
        ),
    )

    parser.add_argument(
        "-f",
        "--filter",
        action="store_false",
        default=True,
        dest="filter",
        help=(
            "Filter signal using a median filter. Default=True. To turn of filtering "
            "do `analyze_mps voltage.czi -f`"
        ),
    )

    parser.add_argument(
        "-t",
        "--averaging_type",
        action="store",
        default="temporal",
        dest="averaging_type",
        choices=["spatial", "temporal", "global"],
        help=(
            "Type of averainging of the images. Only applicable if alpha < 1.0. "
            "Default='temporal'. To change to e.g 'global' do `analyze_mps voltage.czi -t 'global'`"
        ),
    )

    parser.add_argument(
        "-a",
        "--alpha",
        action="store",
        default=1.0,
        dest="alpha",
        type=float,
        help=(
            "When taking the average over the images include only values larger "
            "than the bound when all values are sorted. If alpha = 1.0, then "
            "all values will be used when taking the average. If alpha = 0.1 then "
            "only 10 percent of the pixels will be included."
        ),
    )

    parser.add_argument(
        "--ead-prom",
        action="store",
        default=0.4,
        dest="ead_prominence_threshold",
        type=float,
        help=(
            "How prominent a peak should be in order to be "
            "characterized as an EAD. This value shold be "
            "between 0 and 1, with a greater value being "
            "more prominent. Defaulta: 0.04 "
        ),
    )

    parser.add_argument(
        "--ead-sigma",
        action="store",
        default=3.0,
        dest="ead_sigma",
        type=float,
        help=(
            "Standard deviation in the gaussian smoothing kernal applied "
            "before estimating EADs. Default: 3.0"
        ),
    )

    parser.add_argument(
        "-l",
        "--local",
        action="append",
        dest="local",
        default=[],
        type=int,
        nargs=4,
        help=(
            "If you want to average over a smaller local window in the image you can "
            "provide the indices here. Example `analyze_mps my_file.czi -l 100 120 200 250`, "
            "will average over the window 100 < x < 120, and 200 < y 250"
        ),
    )

    parser.add_argument(
        "-s",
        "--spline",
        dest="use_spline",
        action="store_false",
        default=True,
        help=(
            "Use spline to approximate signals when estimating features. Default=True. "
            "To turn of spline interpolation do `analyze_mps voltage.czi -s`"
        ),
    )

    parser.add_argument(
        "-n",
        "--normalize",
        dest="normalize",
        action="store_true",
        default=False,
        help=(
            "Normalize signals to have maximum value 1 and minimum value 0. "
            "Note that this will change the upstroke velocity and the integral. "
            "You want in general not to do this. "
        ),
    )

    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Print a lot more. Default : False. Example: `analyze_mps voltage.czi -v`",
    )

    parser.add_argument(
        "-e",
        "--std_ex_factor",
        action="store",
        default=1.0,
        dest="std_ex",
        type=float,
        help=(
            "Exclude values outside this factor times 1 standard deviation. "
            "Default:1.0, meaning exclude values outside of 1 std from the mean"
        ),
    )

    parser.add_argument(
        "--spike_duration",
        action="store",
        default=0,
        dest="spike_duration",
        type=int,
        help=(
            "Remove spikes from signal by setting this value "
            "greater than 0. This will locate the timing of "
            "pacing (if available) and delete the signal this "
            "amount after pacing starts. If 0 or no pacing is "
            "detected, nothing will be deleted. "
        ),
    )

    parser.add_argument(
        "--chopping.threshold_factor",
        action="store",
        dest="chopping_threshold_factor",
        default=0.3,
        type=float,
        help=("Factor of where to synchronize data when chopping. Default = 0.3"),
    )

    parser.add_argument(
        "--chopping.extend_front",
        action="store",
        dest="chopping_extend_front",
        type=int,
        help=(
            "How many milliseconds extra to extend signal at the beginning when chopping"
        ),
    )

    parser.add_argument(
        "--chopping.extend_end",
        action="store",
        dest="chopping_extend_end",
        type=int,
        help=("How many milliseconds extra to extend signal at the end when chopping"),
    )

    parser.add_argument(
        "--chopping.min_window",
        action="store",
        dest="chopping_min_window",
        default=200,
        type=int,
        help=(
            "Smallest allowed length of beat (in milliseconds) to be included in chopping"
        ),
    )

    parser.add_argument(
        "--chopping.max_window",
        action="store",
        dest="chopping_max_window",
        default=2000,
        type=int,
        help=(
            "Largest allowed length of beat (in milliseconds) to be included in chopping"
        ),
    )

    parser.add_argument(
        "--overwrite",
        action="store_false",
        help=(
            "If True, overwrite existing data if outdir allready exist. If False, then the old"
            "data will be copied to a subdirectory with version number of the software. If "
            'version number is not found it will be saved to a folder called "old".'
        ),
    )

    parser.add_argument(
        "--reuse_settings",
        action="store_true",
        help=(
            "If the output folder contains a file called settings.json and this flag "
            "is turned on, then we will use the settings stored in the settings.json file. "
            "This is handy if you e.g want to compare the output of the software between different "
            "versions, or you to reproduce the exact traces from the raw data."
        ),
    )

    parser.add_argument(
        "--ignore_pacing",
        action="store_true",
        help=("Ignore pacing data, for example if the pacing is wrong"),
    )

    return parser


class mps_summary(Script):

    logger = mps.utils.get_logger("mps_summary")
    exclude_patterns = brightfield_patterns
    _make_outdir = False

    def check_args(self):
        msg = "The folder '{}' does not exist".format(self.args["folder"])
        assert os.path.exists(self.args["folder"]), msg

        msg = "Given path must be a valid folder"
        assert os.path.isdir(self.args["folder"]), msg

    @staticmethod
    def get_parser():

        descr = "Create a summary pdf of all files in the a directory."
        usage = (
            "Assume all files are located in a folder called 'fiiles', "
            "then Example: \nmps_summary files\n will create a pdf "
            "with a summary of all the experiments in that folder."
        )

        parser = argparse.ArgumentParser(description=descr, usage=usage, add_help=True)

        parser.add_argument(
            action="store", dest="folder", type=str, help="Path to the folder"
        )

        parser.add_argument(
            "-o",
            "--outdir",
            action="store",
            dest="outdir",
            default="",
            type=str,
            help=(
                "Output directory for where you want to store the output. "
                "If not provided a folder with the same name as "
                "the basename of the input file in the currect directory will "
                "be used. Example: `analyze_mps voltage.czi -o figures` will store "
                "the output in the folder 'figures'."
            ),
        )

        parser.add_argument(
            "-f",
            "--filname",
            action="store",
            dest="filename",
            default="mps_summary",
            type=str,
            help=(
                "Name of the pdf and csv file that is the output from the mps_summary "
                "script. Default_ mps_summary"
            ),
        )

        parser.add_argument(
            "-s",
            "--silent",
            dest="silent",
            action="store_true",
            default=False,
            help="Turn off printing",
        )

        parser.add_argument(
            "-i",
            "--ignore_pacing",
            dest="ignore_pacing",
            action="store_true",
            default=False,
            help=(
                "Ignore pacing when computing features. This is usefule it the cells "
                " didn't catch the pacing, but you want to do the analysis anyway"
            ),
        )

        return parser

    @classmethod
    def get_data(cls, f, ignore_pacing=False):

        try:
            cls.logger.info("Processing file {}".format(f))
            data = mps.MPS(f)
            average = mps.analysis.average_intensity(data.frames)
            time_stamps = data.time_stamps
            pacing = data.pacing
            background = mps.analysis.background(time_stamps, average)
            avg = mps.analysis.filt(average - background, 3)

            pacing_chop = np.zeros_like(pacing) if ignore_pacing else pacing
            chopped_data = mps.analysis.chop_data(avg, time_stamps, pacing=pacing_chop)
            apd_analysis = mps.analysis.analyze_apds(
                chopped_data.data, chopped_data.times, plot=False
            )

            tau75s = [
                mps.analysis.tau(t, c, 0.75)
                for c, t in zip(chopped_data.data, chopped_data.times)
            ]
            upstroke80s = [
                mps.analysis.upstroke(t, c, 0.8)
                for c, t in zip(chopped_data.data, chopped_data.times)
            ]
            ttp = [
                mps.analysis.time_to_peak(t, c, p)
                for c, t, p in zip(
                    chopped_data.data, chopped_data.times, chopped_data.pacing
                )
            ]
            nbeats = len(chopped_data.data)
            freqs = mps.analysis.beating_frequency_modified(
                chopped_data.data, chopped_data.times
            )

            return dict(
                average=average,
                time_stamps=time_stamps,
                pacing=pacing,
                apd30s=apd_analysis.apds[30],
                apd80s=apd_analysis.apds[80],
                apd90s=apd_analysis.apds[90],
                capd30s=apd_analysis.capds[30],
                capd80s=apd_analysis.capds[80],
                capd90s=apd_analysis.capds[90],
                slope_APD80=apd_analysis.slope_APD80,
                slope_cAPD80=apd_analysis.slope_cAPD80,
                triangulation=apd_analysis.triangulation,
                tau75s=tau75s,
                upstroke80s=upstroke80s,
                nbeats=nbeats,
                ttp=ttp,
                freqs=freqs,
            )

        except Exception as ex:
            logging.warning(ex, exc_info=True)
            return None

    def main_func(self):

        level = logging.WARNING if self.args["silent"] else logging.INFO
        self.logger.setLevel(level)

        files = [
            f
            for f in os.listdir(self.args["folder"])
            if ((f.endswith("nd2") or f.endswith("czi")) and ("BF" not in f))
        ]

        fig, ax = plt.subplots(len(files), 2)
        xl_datas = []

        for i, f in enumerate(files):
            ax[i, 0].set_title(f.replace("_", r"\_"))
            path = os.path.join(self.args["folder"], f)
            data = mps_summary.get_data(path, self.args["ignore_pacing"])
            if data is None:
                continue

            ax[i, 0].plot(
                data["time_stamps"],
                mps.analysis.normalize_signal(data["average"]),
                color="b",
            )
            ax[i, 0].plot(
                data["time_stamps"],
                mps.analysis.normalize_signal(data["pacing"]),
                label="Pacing ignored!",
                color="r",
            )
            if self.args["ignore_pacing"]:
                ax[i, 0].legend()

            max_std = 20
            apd_keys = ["apd30s", "apd80s", "apd90s", "capd30s", "capd80s", "capd90s"]
            nremoved = {k: 0 for k in apd_keys}
            for k in nremoved.keys():
                while np.std(data[k]) > max_std:
                    # Find the value that is most distant from the mean and remove it
                    idx = np.argmax(np.abs(np.subtract(data[k], np.mean(data[k]))))
                    nremoved[k] += 1
                    data[k].pop(idx)

            keys = [
                "apd30s",
                "apd80s",
                "apd90s",
                "upstroke80s",
                "tau75s",
                "ttp",
                "freqs",
                "slope_APD80",
                "triangulation",
            ]
            labels = [
                "APD30",
                "APD80",
                "APD90",
                "upstroke80",
                "tau75",
                "time to peak",
                "Frequency",
                "Slope APD80",
                "Triangulation",
            ]
            pos = np.linspace(0, 1, len(keys) + 1)[1:]
            for j, (key, label) in enumerate(zip(keys, labels)):

                if key == "slope_APD80":
                    s1 = f"{data[key]:.2f}"
                else:
                    s1 = "{:.2f} +/- {:.2f} ms".format(
                        np.mean(data[key]), np.std(data[key])
                    )
                s = "{}: {}".format(label, s1)
                ax[i, 1].text(0.0, pos[j], s)

            ax[i, 1].axis("off")

            xl_datas.append(
                [
                    f,
                    np.mean(data["apd30s"]),
                    np.std(data["apd30s"]),
                    nremoved["apd30s"],
                    np.mean(data["capd30s"]),
                    np.std(data["capd30s"]),
                    nremoved["capd30s"],
                    np.mean(data["apd80s"]),
                    np.std(data["apd80s"]),
                    nremoved["apd80s"],
                    np.mean(data["capd80s"]),
                    np.std(data["capd80s"]),
                    nremoved["capd80s"],
                    np.mean(data["apd30s"]) / np.mean(data["apd80s"]),
                    np.mean(data["triangulation"]),
                    np.std(data["triangulation"]),
                    np.mean(data["apd90s"]),
                    np.std(data["apd90s"]),
                    nremoved["apd90s"],
                    np.mean(data["upstroke80s"]),
                    np.std(data["upstroke80s"]),
                    np.mean(data["tau75s"]),
                    np.std(data["tau75s"]),
                    np.mean(data["ttp"]),
                    np.std(data["ttp"]),
                    data["nbeats"],
                    np.mean(data["freqs"]),
                    np.std(data["freqs"]),
                    data["slope_APD80"],
                    data["slope_cAPD80"],
                ]
            )

        header = [
            "Filename",
            "APD30 [ms]",
            "APD30 stdev",
            "APD30 #removed",
            "cAPD30 [ms]",
            "cAPD30 stdev",
            "cAPD30 #removed",
            "APD80 [ms]",
            "APD80 stdev",
            "APD80 #removed",
            "cAPD80 [ms]",
            "cAPD80 stdev",
            "cAPD80 #removed",
            "ratio APD30/APD80" "APD90 [ms]",
            "triangulation [ms]",
            "triangulation stdev",
            "APD90 [ms]",
            "APD90 stdev",
            "APD90 #removed",
            "Upstroke80 [ms]",
            "Upstroke80 stdev",
            "tau75 [ms]",
            "tau75 stdev",
            "Time to peak [ms]",
            "Time to peak stdev",
            "Number of beats",
            "Frequency [Hz]",
            "Frequency stdev",
            "Slope APD80",
            "Slope cAPD80",
        ]

        mps.utils.to_csv(
            xl_datas, os.path.join(self.args["folder"], self.args["filename"]), header
        )
        fig.set_figheight(5 * len(files))
        fig.set_figwidth(15)
        fig.savefig(
            os.path.join(self.args["folder"], "{}.pdf".format(self.args["filename"])),
            dpi=100,
        )


class mps2mp4(Script):

    logger = mps.utils.get_logger("mps2mp4")

    def check_args(self):

        msg = "The file {} does not exist".format(self.args["file"])
        assert os.path.isfile(self.args["file"]), msg

        if self.args["outfile"] == "":
            self.args["outfile"] = os.path.splitext(self.args["file"])[0]

    @staticmethod
    def get_parser():

        descr = "Convert stack of czi files containing MPS data into movie"
        usage = "Example: \n./mps2mp4 1uM 1hz berst.czi -o figures"
        parser = argparse.ArgumentParser(description=descr, usage=usage, add_help=True)

        parser.add_argument(
            action="store", dest="file", type=str, help="Path to the mps file"
        )

        parser.add_argument(
            "-o",
            "--outfile",
            action="store",
            dest="outfile",
            default="",
            type=str,
            help=(
                "Output name for where you want to store the output movie. "
                "If not provided a the same name as the basename of the czi file will be used."
            ),
        )

        parser.add_argument(
            "-s",
            "--synch",
            action="store_true",
            dest="synch",
            help=("Start video at same time as start of pacing"),
        )

        return parser

    def main_func(self):

        mps_data = mps.MPS(self.args["file"])
        if self.args["synch"]:
            idx = next(i for i, p in enumerate(mps_data.pacing) if p > 0)
        else:
            idx = 0

        mps.utils.frames2mp4(
            mps_data.frames.T[idx:, :, :], self.args["outfile"], mps_data.framerate
        )
        self.logger.info("Movie saved to {}".format(self.args["outfile"]))


class collect_mps(Script):

    logger = mps.utils.get_logger("collect_mps")
    _make_outdir = True
    exclude_patterns = brightfield_patterns

    @staticmethod
    def get_parser():

        descr = (
            "Gather Voltage and Calcium data into one file that can be "
            "used as input to the inversion algorithm."
        )
        usage = (
            "Example: \ncollect_mps '1uM 1hz berst' '1uM 1hz gcamp' -o '1uM 1hz'\n"
            + "Note that first argument should be voltage and second should be calcium"
        )
        parser = argparse.ArgumentParser(description=descr, usage=usage, add_help=True)

        parser.add_argument(
            action="store",
            dest="voltage",
            type=str,
            help="Path to the directory containing the voltage output",
        )

        parser.add_argument(
            action="store",
            dest="calcium",
            type=str,
            help="Path to the directory containing the calcium output",
        )

        parser.add_argument(
            "-o",
            "--outdir",
            action="store",
            dest="outdir",
            default="",
            type=str,
            help=(
                "Output directory for where you want to store the output. "
                "If not provided a folder with the same name as "
                "the basename of the czi file in the currect directory will "
                "be used"
            ),
        )

        _collect_mps_args(parser)
        return parser

    def check_args(self):

        if self.args["output_format"] == "yml":
            try:
                import yaml
            except ImportError:
                msg = "Cannot save to yml format. yaml is not installed"
                self.logger.error(msg)
                raise ImportError(msg)

        for name, s in zip(
            ["voltage", "calcium"], [self.args["voltage"], self.args["calcium"]]
        ):

            msg = "The folder {} does not exist".format(s)
            assert os.path.isdir(s), msg

            f = os.path.join(s, "data.npy")
            msg = (
                "The folder {0} does not contain any data files. "
                "Please rerun analyze_mps {0}.czi"
            ).format(s)
            assert os.path.isfile(f), msg

            self.logger.info("Found {} input: {}".format(name, s))

        if self.args["outdir"] == "":
            self.logger.warning("Output directory not specified")
            self.args["outdir"] = os.path.join(
                os.path.dirname(self.args["voltage"]),
                "_".join(
                    [
                        os.path.basename(self.args["voltage"]),
                        os.path.basename(self.args["calcium"]),
                    ]
                ),
            )
        self.logger.info(
            "Output will be saved to {}".format(os.path.abspath(self.args["outdir"]))
        )
        if not os.path.isdir(self.args["outdir"]):
            os.makedirs(self.args["outdir"])

    def main_func(self):

        self.logger.info("Collect voltage and calcium data...")
        voltage_data = np.load(
            os.path.join(self.args["voltage"], "data.npy"), allow_pickle=True
        ).item()
        calcium_data = np.load(
            os.path.join(self.args["calcium"], "data.npy"), allow_pickle=True
        ).item()

        data = mps.utils.collect_data(voltage_data, calcium_data)

        # Plot the data
        fig, ax = plt.subplots(1, 2)
        ax[0].plot(data["t_V"][: len(data["V"])], data["V"])
        ax[0].set_title("Voltage")
        ax[0].set_ylabel("Pixel intensity")
        ax[1].plot(data["t_Ca"][: len(data["Ca"])], data["Ca"])
        ax[1].set_title("Calcium")
        fig.savefig(os.path.join(self.args["outdir"], "input_data.pdf"))
        plt.close()

        fig, ax = plt.subplots(2, 1, sharex=True)
        try:
            v = voltage_data["unchopped_data"]["trace"]
            c = calcium_data["unchopped_data"]["trace"]
        except KeyError:
            v = voltage_data["unchopped_data"]["data"]
            c = calcium_data["unchopped_data"]["data"]

        p_v = voltage_data["unchopped_data"]["pacing"]
        p_c = calcium_data["unchopped_data"]["pacing"]

        ax[0].plot(voltage_data["unchopped_data"]["time"][: len(v)], v, "b-")

        ax[0].set_title("Voltage")
        ax[0].set_ylabel("Pixel intensity")
        ax0 = ax[0].twinx()
        ax0.plot(voltage_data["unchopped_data"]["time"][: len(p_v)], p_v, "r-")
        ax0.set_ylabel("Pacing")
        ax[1].plot(calcium_data["unchopped_data"]["time"][: len(c)], c, "b-")
        ax[1].set_title("Calcium")
        ax[1].set_ylabel("Pixel intensity")
        ax[1].set_xlabel("Time ({})".format(voltage_data["attributes"]["time_unit"]))
        ax1 = ax[1].twinx()
        ax1.plot(calcium_data["unchopped_data"]["time"][: len(p_c)], p_c, "r-")
        ax1.set_ylabel("Pacing")
        fig.savefig(os.path.join(self.args["outdir"], "full_traces.pdf"))
        plt.close()

        # Path to output file (without extension)
        out = os.path.join(self.args["outdir"], self.args["output_name"])
        # Matfiles cannot reade datetime format to remove this
        if self.args["output_format"] == "mat":
            for k in ["attributes_V", "attributes_Ca"]:
                data[k].pop("date", None)
        # Dump data
        mps.utils.dump_data(data, out, self.args["output_format"])


def _collect_mps_args(parser):

    parser.add_argument(
        "-f",
        "--output_format",
        action="store",
        dest="output_format",
        default="npy",
        type=str,
        choices=["mat", "npy", "yml"],
        help="Format for the output data",
    )

    parser.add_argument(
        "-n",
        "--output_name",
        action="store",
        dest="output_name",
        default="output",
        type=str,
        help="Name of outputfile",
    )

    parser.add_argument(
        "-s",
        "--synch",
        action="store_false",
        dest="synch",
        default=True,
        help=(
            "If True such voltage and calcium data based on pacing data, so that "
            + "the pacing protocol is the same. If no pacing data is provided "
            + "then no synching will be performed. Default: True"
        ),
    )

    return parser


class mps_phase_plot(Script):

    logger = mps.utils.get_logger("mps_phase_plot")
    exclude_patterns = brightfield_patterns

    @staticmethod
    def get_parser():

        descr = (
            "Make a phase plot with voltage on the x-axis and " "calcium on the y-axis."
        )

        parser = argparse.ArgumentParser(description=descr, add_help=True)

        parser.add_argument(
            action="store", dest="voltage", type=str, help="Path to the voltage file"
        )

        parser.add_argument(
            action="store", dest="calcium", type=str, help="Path to the calcium file"
        )

        parser.add_argument(
            "-o",
            "--outfile",
            action="store",
            dest="outfile",
            default="",
            type=str,
            help=(
                "Name of the output file. If not specified it will merge "
                "the name of the voltage and calcium file"
            ),
        )

        return parser

    def check_args(self):

        for name, s in zip(
            ["voltage", "calcium"], [self.args["voltage"], self.args["calcium"]]
        ):

            msg = "The file {} does not exist".format(s)
            assert os.path.isfile(s), msg

        if self.args["outfile"] == "":
            self.logger.warning("Output file not specified")
            self.args["outfile"] = os.path.join(
                os.path.dirname(self.args["voltage"]),
                "_".join(
                    [
                        os.path.splitext(os.path.basename(self.args["voltage"]))[0],
                        os.path.splitext(os.path.basename(self.args["calcium"]))[0],
                    ]
                ),
            )
        outfile = os.path.abspath(self.args["outfile"])
        self.logger.info("Output will be saved to {}".format(outfile))
        outdir = os.path.dirname(outfile)
        if not os.path.isdir(outdir):
            os.makedirs(outdir)

    def main_func(self):

        self.logger.info("Plot phase plot of voltage and calcium")
        voltage = mps.MPS(self.args["voltage"])
        calcium = mps.MPS(self.args["calcium"])

        voltage_data = mps.analysis.analyze_mps_func(voltage, plot=False)
        calcium_data = mps.analysis.analyze_mps_func(calcium, plot=False)

        from . import plotter
        from . import analysis

        voltage_trace = voltage_data["chopped_data"]["trace_1std"]
        calcium_trace = calcium_data["chopped_data"]["trace_1std"]

        plotter.phase_plots(
            analysis.normalize_signal(voltage_trace),
            analysis.normalize_signal(calcium_trace),
            self.args["outfile"],
        )


class mps_prevalence(Script):

    logger = mps.utils.get_logger("mps_prevalence")
    _make_outdir = True
    exclude_patterns = brightfield_patterns

    @staticmethod
    def get_parser():

        descr = (
            "Estimate the percentage of tissue in the chips, "
            "and the percentage of beating tissue vs non-beating\n\n"
        )

        parser = argparse.ArgumentParser(description=descr, add_help=True)

        parser.add_argument(
            action="store", dest="file", type=str, help="Path to the file"
        )

        parser.add_argument(
            "-o",
            "--outdir",
            action="store",
            dest="outdir",
            default="",
            type=str,
            help=(
                "Output directory for where you want to store the output. "
                "If not provided a folder with the same name as "
                "the basename of the input file in the currect directory will "
                "be used. Example: `analyze_mps voltage.czi -o figures` will store "
                "the output in the folder 'figures'."
            ),
        )

        return _mps_prevalence_args(parser)

    def check_args(self):

        msg = "The file {} does not exist".format(self.args["file"])
        if not os.path.isfile(self.args["file"]):
            raise IOError(msg)

    def main_func(self):
        outdir = make_outdir(self.args)
        mps_data = mps.MPS(self.args["file"])
        mps_prevalence = mps.analysis.prevalence(mps_data, **self.args)

        prevalence = np.zeros(mps_prevalence.is_tissue.shape, dtype=int)
        prevalence[np.where(mps_prevalence.is_tissue & ~mps_prevalence.is_beating)] = 1
        prevalence[np.where(mps_prevalence.is_tissue & mps_prevalence.is_beating)] = 2

        out_string = (
            f"Prevalence: {mps_prevalence.prevalence * 100:.1f}\\%\n"
            f"Tissue covered area: {mps_prevalence.tissue_covered_area * 100:.1f}\\%"
        )
        self.logger.info("\n" + out_string)
        from matplotlib.colors import LinearSegmentedColormap

        colors = [(0, 0, 0), (0.8, 0.8, 0.8), (1, 0, 0)]
        cm = LinearSegmentedColormap.from_list("prevalence", colors, N=3)

        fig, ax = plt.subplots(1, 2)
        ax[0].imshow(mps_data.frames.T[0].T)
        cax = ax[1].imshow(prevalence, cmap=cm)
        ax[1].set_title(out_string)

        cbar = fig.colorbar(cax, ticks=[0, 1, 2])
        cbar.ax.set_yticklabels(["Not tissue", "Non beating", "Beating"])
        fig.savefig(os.path.join(outdir, "prevalence"))


def _mps_prevalence_args(parser):

    parser.add_argument(
        "-N",
        action="store",
        default=50,
        dest="N",
        type=int,
        help=(
            """Size of grid along the major axis (minor axis will
            be scaled proportionally). Default: 50"""
        ),
    )

    parser.add_argument(
        "-bt",
        "--baseline-threshold",
        action="store",
        default=0.1,
        dest="baseline_threshold",
        type=float,
        help=(
            """Percentage (between 0 and 1) of how far from the min and
            max values of the baseline intensities of the beating regions
            that should define the rage for baseline intensity of tissue
            Default: 0.1"""
        ),
    )

    parser.add_argument(
        "-ft",
        "--frequency-threshold",
        action="store",
        default=0.2,
        dest="frequency_threshold",
        type=float,
        help=(
            """Percentage (between 0 and 1) of how far from the global
            frequency a local signal can be before it is classified
            as non-beating. Default: 0.2"""
        ),
    )

    parser.add_argument(
        "-s",
        "--snr-factor",
        action="store",
        default=1.5,
        dest="snr_factor",
        type=float,
        help=(
            """Factor multiplied with the global signal to noise ratio (snr),
            to determine wether a region is noise or not. If a local region
            has a larger values than snr_factor * global snr it will be
            classied as noise. Default: 1.5"""
        ),
    )
    return parser


class analyze_med64(Script):

    logger = mps.utils.get_logger("analyze_med64")
    valid_extensions = [".modat"]
    skip_startswith = ["."]
    skip_endswith = ["~", "#"]
    _make_outdir = True

    @staticmethod
    def get_parser():

        descr = "Analyze MED64 data.\n\n"

        parser = argparse.ArgumentParser(description=descr, add_help=True)

        parser.add_argument(
            action="store",
            dest="file",
            type=str,
            help="Path to the file MED64 file (.modat or converted file)",
        )

        parser.add_argument(
            "-o",
            "--outdir",
            action="store",
            dest="outdir",
            default="",
            type=str,
            help=(
                "Output directory for where you want to store the output. "
                "If not provided a folder with the same name as "
                "the basename of the input file in the currect directory will "
                "be used. Example: `analyze_mps voltage.czi -o figures` will store "
                "the output in the folder 'figures'."
            ),
        )

        return _analyze_med64_args(parser)

    def _check_args(self, path):

        if self.args["skip"]:
            outdir = os.path.splitext(path)[0]
            data_path = os.path.join(outdir, "data." + self.args["output_format"])
            if os.path.isfile(data_path):
                return False
        return super()._check_args(path)

    def main_func(self):
        def func(path, outdir):
            if not self._check_args(path):
                return
            self.logger.info(f"Analyzing file {path}")

            datapath = os.path.join(outdir, "data.") + self.args["output_format"]
            if self.args["reset"] or not os.path.isfile(datapath):

                with mps.MED64(path) as med64:
                    data = dict(
                        positions=med64.electrode_positions(),
                        electrode_arrangement=med64.electrode_arrangement,
                        data=med64._data,
                        times=med64.times,
                        factor=med64.metadata["factor"],
                    )
                if self.args["output_format"] == "mat":
                    scipy.io.savemat(datapath, data)
                else:
                    np.save(datapath, data, allow_pickle=True)

            if self.args["output_format"] == "mat":
                data = mps.utils.loadmat(datapath)
            else:
                data = np.load(datapath, allow_pickle=True).item()

            if not self.args["noplot"]:
                figpath = os.path.join(outdir, "data")
                mps.med64_utils.plot_data(figpath, **data)

            if self.args["with_cv"]:
                filename = os.path.join(outdir, "conduction_velocity")
                cv = mps.med64_utils.analyze_med64(
                    **data, **self.args, filename=filename
                )
                self.logger.info(f"Estimated conduction velocity is {cv:.2f} cm / s")

        super().main_func(func)


def _analyze_med64_args(parser):

    parser.add_argument(
        "-f",
        "--output_format",
        action="store",
        dest="output_format",
        default="mat",
        type=str,
        choices=["mat", "npy"],
        help="Format for the output data",
    )

    parser.add_argument(
        "--noplot",
        action="store_true",
        help=("Turn off plotting of contraction data and mean contraction."),
    )

    parser.add_argument(
        "--skip",
        action="store_true",
        help=("Skip files that have been analyzed before"),
    )

    parser.add_argument(
        "--reset",
        action="store_true",
        help=(
            "Load modat file again, even if it is allready saved in 'mat' or 'npy' format."
        ),
    )

    parser.add_argument(
        "--with_cv", action="store_true", help=("Compute conduction velocity as well")
    )

    parser.add_argument(
        "--cv_alg",
        dest="cv_algorithm",
        default="gradient",
        type=str,
        choices=mps.med64_utils.CV_ALGORITHMS,
        help=(
            "Choice of algorithm to compute conduction velocity. "
            "'gradient' (Default) will compute the speed of the wave "
            "in the direction of the gradient at each "
            "electrode and return the median value. "
            "'optimized' will first try to find the sources of the wave by selecting "
            "the points with the lowest time of activation, and then find the "
            "conduction velocity that minimizes the difference betweed model and data."
        ),
    )

    return parser


class mps_deltaF_summary(Script):

    logger = mps.utils.get_logger("mps_deltaF_summary")
    exclude_patterns = brightfield_patterns
    _make_outdir = False

    def check_args(self):
        msg = "The folder '{}' does not exist".format(self.args["folder"])
        assert os.path.exists(self.args["folder"]), msg

        msg = "Given path must be a valid folder"
        assert os.path.isdir(self.args["folder"]), msg

        if self.args["outdir"] == "":
            self.args["outdir"] = self.args["folder"]

    @staticmethod
    def get_parser():

        descr = "Create a summary csv of all folders' deltaF max in the a directory."
        usage = ""

        parser = argparse.ArgumentParser(description=descr, usage=usage, add_help=True)

        parser.add_argument(
            action="store", dest="folder", type=str, help="Path to the folder"
        )

        parser.add_argument(
            "-o",
            "--outdir",
            action="store",
            dest="outdir",
            default="",
            type=str,
            help=(
                "Output directory for where you want to store the output. "
                "If not provided a folder with the same name as "
                "the basename of the input folder."
            ),
        )

        parser.add_argument(
            "-f",
            "--filname",
            action="store",
            dest="filename",
            default="delta_F_summary",
            type=str,
            help=(
                "Name of the pdf and csv file that is the output from the mps_summary "
                "script. Default: delta_F_summary"
            ),
        )

        return parser

    def main_func(self):
        xl_datas = []
        for root, subfolders, files in os.walk(self.args["folder"]):
            for numpyFile in files:
                if numpyFile.endswith(".npy"):
                    path = os.path.join(
                        root, numpyFile
                    )  # join the path of current folder, and the file in subfolder
                    data = np.load(path, allow_pickle=True).item()  # load the .npy file

                    if "chopped_data" not in data:
                        continue
                    if "trace_1std" not in data["chopped_data"]:
                        continue
                    xl_datas.append(
                        [
                            path,
                            numpyFile,
                            np.max(
                                data["chopped_data"]["trace_1std"]
                            ),  # add the deltaF max
                            np.std(
                                data["chopped_data"]["trace_1std"]
                            ),  # add the deltaF std
                        ]
                    )

        header = ["path", "Filename", "deltaF max", "deltaF stdev"]

        path = os.path.join(self.args["outdir"], self.args["filename"])
        mps.utils.to_csv(xl_datas, path, header)
