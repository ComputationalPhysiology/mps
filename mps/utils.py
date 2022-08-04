#!/usr/bin/env python3
__author__ = "Henrik Finsberg (henriknf@simula.no), 2017--2020"
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

import datetime
import logging
import os
import sys

# Suppress scipy warning
import warnings
from collections import OrderedDict, namedtuple
from pathlib import Path

import numpy as np
import scipy.io as sio

warnings.filterwarnings(action="ignore", module="scipy", message="^internal gelsd")

grid_settings = namedtuple(
    "grid_settings",
    ("nx, ny, dx, dy, x_start," " y_start, x_end, y_end"),
)

logger = logging.getLogger(__name__)


def get_space_step(N, x_start, x_end, y_start, y_end, **kwargs):

    frames = kwargs.get("frames", None)
    if x_end is None:
        if frames is None:
            raise ValueError("Need to provide x_end or frames")
        x_end = np.shape(frames)[0]
    if y_end is None:
        if frames is None:
            raise ValueError("Need to provide y_end or frames")
        y_end = np.shape(frames)[1]

    max_length = max(x_end - x_start, y_end - y_start)
    return int(max_length / N)


def merge_pdfs(lst, fname, cleanup=True):
    """
    Given a list of paths to pdfs, merge these together and save
    the results in a new file

    lst : list
        List with paths to existing pdfs
    fname : str
        Name of the output file where you want to save the
        merged pdf.
    cleanup : bool
        If True, delete the files in the list provided.
        Default: True
    """

    from PyPDF2 import PdfFileReader, PdfFileWriter

    pdf_writer = PdfFileWriter()

    for f in lst:

        logger.debug(f"Read {f}")
        assert Path(f).is_file(), f"File {f} does not exist"
        pdf_reader = PdfFileReader(f)

        for page in range(pdf_reader.getNumPages()):
            pdf_writer.addPage(pdf_reader.getPage(page))

    fname = Path(fname)
    if fname.suffix != ".pdf":
        fname = Path(str(fname) + ".pdf")
    fname.parent.mkdir(parents=True, exist_ok=True)

    with open(fname, "wb") as fh:
        pdf_writer.write(fh)

    if cleanup:
        for f in lst:
            logger.debug(f"Delete {f}")
            Path(f).unlink()


def get_grid_settings(
    N, x_start=0, y_start=0, x_end=None, y_end=None, frames=None, **kwargs
):
    if x_end is None:
        if frames is None:
            raise ValueError("Need to provide x_end or frames")
        x_end = np.shape(frames)[0]
    if y_end is None:
        if frames is None:
            raise ValueError("Need to provide y_end or frames")
        y_end = np.shape(frames)[1]

    dx = dy = get_space_step(N, x_start, x_end, y_start, y_end, frames=frames)

    nx = max(int((x_end - x_start) / dx), 1)
    ny = max(int((y_end - y_start) / dy), 1)

    return grid_settings(
        nx=nx,
        ny=ny,
        dx=dx,
        dy=dy,
        x_start=x_start,
        y_start=y_start,
        x_end=x_end,
        y_end=y_end,
    )


def namedtuple2dict(x):
    keys = [
        attr
        for attr in dir(x)
        if not attr.startswith("_") and attr not in ["index", "count"]
    ]
    return {k: getattr(x, k) for k in keys}


def get_data_from_dict(set_data, set_key, get_data, get_key, get_name):
    """Copy data from one dictionary to another"""

    try:
        set_data[set_key] = get_data[get_key]
        return 0, ""
    except KeyError:
        msg = "Could not found key {} in {} data".format(get_key, get_name)
        return 1, msg


def find_num_beats(data):

    n = 0
    while True:
        if f"trace_{n}" in data["chopped_data"]:
            n += 1
        else:
            break
    return n


def collect_data(voltage_data, calcium_data):
    """
    Collect voltage and calcium data togther so that it can be used
    in the inversion.
    """

    if (
        "pacing" in voltage_data["unchopped_data"].keys()
        and "pacing" in calcium_data["unchopped_data"].keys()
    ):

        synch_data = True
    else:
        synch_data = False

    data = {}

    # Keys for the features
    keys = ["apd30", "apd50", "apd80", "dFdt_max", "int30"]
    voltage_keys = ["APD_V_30", "APD_V_50", "APD_V_80", "dVdt_max", "int30V"]
    calcium_keys = ["APD_Ca_30", "APD_Ca_50", "APD_Ca_80", "dCadt_max", "int30Ca"]

    if synch_data:
        # First we check that the pacing info is the same
        # for voltage and calcium

        pacing_voltage = voltage_data["chopped_data"]["pacing_1std"]
        pacing_calcium = calcium_data["chopped_data"]["pacing_1std"]

        min_len = min(len(pacing_calcium), len(pacing_voltage))
        pacing_calcium = pacing_calcium[-min_len:]
        pacing_voltage = pacing_voltage[-min_len:]
        # Check that the pacing traces are alinged
        if not (np.array(pacing_voltage) == np.array(pacing_calcium)).all():
            msg = "Pacing traces are not aligned."
            logger.error(msg)
            data["pacing"] = np.zeros(len(voltage_data["chopped_data"]["time_1std"]))

        else:
            data["pacing"] = pacing_voltage

    data["attributes_V"] = voltage_data["attributes"]
    data["attributes_Ca"] = calcium_data["attributes"]
    data["chopping_parameters_V"] = voltage_data["chopping_parameters"]
    data["chopping_parameters_Ca"] = calcium_data["chopping_parameters"]
    data["t_V"] = voltage_data["chopped_data"]["time_1std"]
    data["t_Ca"] = calcium_data["chopped_data"]["time_1std"]

    n_V = find_num_beats(voltage_data)
    data["chopped_data_V"] = dict(
        chopped_data=[voltage_data["chopped_data"][f"trace_{i}"] for i in range(n_V)],
        chopped_times=[voltage_data["chopped_data"][f"time_{i}"] for i in range(n_V)],
    )
    n_Ca = find_num_beats(calcium_data)
    data["chopped_data_Ca"] = dict(
        chopped_data=[calcium_data["chopped_data"][f"trace_{i}"] for i in range(n_Ca)],
        chopped_times=[calcium_data["chopped_data"][f"time_{i}"] for i in range(n_Ca)],
    )
    try:
        data["V"] = voltage_data["chopped_data"]["trace_1std"]
        data["Ca"] = calcium_data["chopped_data"]["trace_1std"]
    except KeyError:
        # Make backwards compatible
        data["V"] = voltage_data["chopped_data"]["avg_1std"]
        data["Ca"] = calcium_data["chopped_data"]["avg_1std"]

    # Collect features
    logger.info("Collect features")
    for ck, vk, k in zip(calcium_keys, voltage_keys, keys):
        logger.info("Key {} - Voltage: {}, Calcium: {}".format(k, vk, ck))
        # Voltage
        exit_code, msg = get_data_from_dict(
            data,
            vk,
            voltage_data["features"],
            k,
            "voltage",
        )
        if exit_code:
            # Backward compatible
            exit_code, _ = get_data_from_dict(
                data,
                vk,
                voltage_data["features"],
                k.upper(),
                "voltage",
            )
            if exit_code:
                logger.warning(msg)

        # Calcium
        exit_code, msg = get_data_from_dict(
            data,
            ck,
            calcium_data["features"],
            k,
            "calcium",
        )
        if exit_code:
            # Backward compatible
            exit_code, _ = get_data_from_dict(
                data,
                ck,
                calcium_data["features"],
                k.upper(),
                "calcium",
            )
            if exit_code:
                logger.warning(msg)

    return data


def to_csv(data, path, header=None):
    """
    FIXME
    """

    if isinstance(data, dict):
        header = data.keys()
        data = list(zip(*data.values()))

    import csv

    if sys.platform == "win32":
        kwargs = {}
    else:
        kwargs = dict(delimiter=";")

    with open(Path(path).with_suffix(".csv"), "w", newline="") as f:
        w = csv.writer(f, dialect="excel", **kwargs)
        if header is not None:
            w.writerow(header)
        for i in data:
            w.writerow(i)
    logger.debug("Saved to {}.csv".format(path))


def json_serial(obj):
    """JSON serializer for objects not serializable by default json code"""

    if isinstance(obj, (datetime.datetime)):
        return obj.isoformat()
    elif isinstance(obj, (np.ndarray)):
        return obj.tolist()
    else:
        try:
            return str(obj)
        except Exception:
            raise TypeError("Type %s not serializable" % type(obj))


def to_txt(data, path, header_list=None, delimiter=";", fmt="%10.6g"):

    if header_list is not None:
        header = delimiter.join(
            [("{:>" + fmt[1:].split(".")[0] + "s}").format(h) for h in header_list],
        )

    else:
        header = ""

    np.savetxt(path, data, delimiter=delimiter, header=header, fmt=fmt, comments="")


def dump_data(data, path, fileformat="npy"):
    """
    Dump data to a file with given file format

    Arguments
    ---------
    data : dict or array
        The data you want to dum
    path : str
        Path to the data (without extension)
    fileformat : str
        File format. Either 'npy', 'yml', 'csv' or 'mat'.
        Deafult: 'npy'

    """

    err = False
    if fileformat == "mat":
        from scipy.io import savemat

        savemat(path, data)

    elif fileformat == "npy":
        np.save(path, data)

    elif fileformat == "csv":

        import csv

        try:
            with open(path + ".csv", "w") as f:

                w = csv.writer(f, delimiter=";")

                def write(k, v):

                    if isinstance(v, dict):
                        for k1, v1 in v.items():
                            write(k1, v1)
                    elif isinstance(v, np.ndarray):
                        w.writerow([k])
                        w.writerow(v.tolist())
                    elif isinstance(v, list):

                        if len(v) > 0 and isinstance(v[0], np.ndarray):
                            w.writerow([k])
                            w.writerow(np.array(v).tolist())

                        else:
                            w.writerow([k, v])
                    else:
                        w.writerow([k, v])

                for k, v in data.items():
                    write(k, v)

        except IOError as exc:
            logger.error("I/O error: {}".format(exc))
            err = True

    else:
        err = True
        logger.error("Unknown fileformat {}".format(fileformat))

    if not err:
        logger.info("Data dumped to {}.{}".format(path, fileformat))


def get_intersection(data):
    """
    Get intersection of all values in
    a dictionary
    """
    if isinstance(data, dict):
        vals = list(data.values())
    else:
        vals = data
    if len(vals) > 0:
        return list(set(vals[0]).intersection(*list(vals)))
    return vals


def normalize_frames(X, max_val=255, min_val=0, dtype=np.uint8):
    x_min = np.min(X)
    arr = ((X - x_min) / (np.max(X) - x_min)) * (max_val - min_val) + min_val
    return np.array(arr, dtype=dtype)


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


def frames2mp4(frames, path, framerate=None):
    try:
        import imageio
    except ImportError:
        logger.warning(
            '"imageio" not found. Please install imageio in order to'
            'save frames to mp4 file: "pip install imageio"',
        )
    else:

        if framerate is None:
            # Assume one second in total
            framerate = frames.shape[-1]

        mp4path = "{}.mp4".format(os.path.splitext(path)[0])
        try:
            imageio.mimwrite(mp4path, normalize_frames(frames), fps=framerate)
        except RuntimeError:
            imageio.plugins.ffmpeg.download()
            imageio.mimwrite(mp4path, normalize_frames(frames), fps=framerate)

        logger.info("Video saved to {}".format(mp4path))


def loadmat(filename):
    """Short summary.

    Parameters
    ----------
    filename : type
        Description of parameter `filename`.

    Returns
    -------
    type
        Description of returned object.

    """
    """
    this function should be called instead of direct sio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    """

    def _check_keys(d):
        """
        checks if entries in dictionary are mat-objects. If yes
        todict is called to change them to nested dictionaries
        """
        for key in d:
            if isinstance(d[key], sio.matlab.mio5_params.mat_struct):
                d[key] = _todict(d[key])
        return d

    def _todict(matobj):
        """
        A recursive function which constructs from matobjects nested dictionaries
        """
        d = {}
        for strg in matobj._fieldnames:
            elem = matobj.__dict__[strg]
            if isinstance(elem, sio.matlab.mio5_params.mat_struct):
                d[strg] = _todict(elem)
            elif isinstance(elem, np.ndarray):
                d[strg] = _tolist(elem)
            else:
                d[strg] = elem
        return d

    def _tolist(ndarray):
        """
        A recursive function which constructs lists from cellarrays
        (which are loaded as numpy ndarrays), recursing into the elements
        if they contain matobjects.
        """
        elem_list = []
        for sub_elem in ndarray:
            if isinstance(sub_elem, sio.matlab.mio5_params.mat_struct):
                elem_list.append(_todict(sub_elem))
            elif isinstance(sub_elem, np.ndarray):
                elem_list.append(_tolist(sub_elem))
            else:
                elem_list.append(sub_elem)
        return elem_list

    data = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)
