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

import io
import logging
import multiprocessing
import os
import time
import zipfile
from collections import namedtuple
from concurrent.futures import ThreadPoolExecutor
from copy import deepcopy

import numpy as np

from . import analysis, czifile, utils
from .nd2file import ND2File

try:
    import tifffile
except ImportError:
    has_tifffile = False
else:
    has_tifffile = True


logger = utils.get_logger(__name__)


mps_data_fields = ["frames", "time_stamps", "info", "metadata", "pacing"]
mps_data = namedtuple(  # type: ignore
    "mps_data", mps_data_fields, defaults=(None,) * len(mps_data_fields)  # type: ignore
)


def info_dictionary(time_stamps):
    dt = np.mean(np.diff(time_stamps))
    time_unit = analysis.time_unit(time_stamps)

    info = dict(num_frames=len(time_stamps), dt=dt, time_unit=time_unit)
    return info


def get_single_frame(path, index=0):
    """
    Get a single frame from a czi or nd2 file.

    Arguments
    ---------
    path : str
        Path to the file
    index : int
        The index of the frame (Default: 0)

    Returns
    -------
    frame : numpy.ndarray
        The frame
    """

    ext = os.path.splitext(path)[1]
    if ext == ".nd2":
        with ND2File(path) as frames:
            frame = frames.image(index).squeeze()

    elif ext == ".czi":
        with czifile.CziFile(path) as f:
            frames = f.asarray().squeeze()
        frame = frames[:, :, index]

    else:
        raise ValueError("Unkown extension {}".format(ext))

    return frame


def load_nd2(fname):
    """
    Load ND2 file (Nikon image file)

    Arguments
    ---------
    nd2file : str
        Path to the file

    Returns
    -------
    I : array
        The images
    attributes : dict
        Dictionary with metadata

    """
    logger.info("Load nd2 file {}".format(fname))
    max_workers = multiprocessing.cpu_count() // 2
    t0 = time.time()
    with ND2File(fname) as frames:

        num_frames = frames.imagecount
        metadata = frames.metadata

        images = np.zeros((frames.height, frames.width, num_frames), dtype=np.uint16)
        time_stamps = np.zeros(num_frames)

        def func(t):
            # print(t)
            images[:, :, t] = frames.image(t).squeeze()
            time_stamps[t] = frames.get_time(t)

        with ThreadPoolExecutor(max_workers) as executor:
            executor.map(func, np.arange(num_frames))

    t1 = time.time()
    logger.info("Loaded {} frames in {} seconds".format(num_frames, t1 - t0))
    # Let time start at zero
    if np.all(time_stamps == 0):
        # Then there is something wrong with the file
        # Assume standard framrate
        if "BF" in fname:  # BrightField
            dt = 25.0
        else:
            dt = 10.0
        time_stamps = np.arange(dt * num_frames, step=dt)

    time_stamps -= time_stamps[0]

    # Make sure time is in milliseconds
    time_unit = analysis.time_unit(time_stamps)
    if time_unit == "s":
        time_stamps *= 1000

    info = info_dictionary(time_stamps)
    info["um_per_pixel"] = metadata["ImageMetadataSeqLV|0"]["SLxPictureMetadata"][
        "dCalibration"
    ]
    if info["um_per_pixel"] == 0.0:
        # We need to get this info in another way

        # This is not correct
        dic = metadata["ImageMetadataSeqLV|0"]["SLxPictureMetadata"]["sPicturePlanes"][
            "sSampleSetting"
        ]["a0"]["pCameraSetting"]["FormatQuality"]["fmtDesc"]
        # pixels = [dic['sizeSensorPixels'][i] for i in ["cx", "cy"]]
        # microns = [dic['sizeSensorMicrons'][i] for i in ["cx", "cy"]]
        logger.warning(
            (
                "Could not find calcibration factor to convert "
                "from pixels to micrometers. Use 0.325 "
                "micrometers per pixel"
            )
        )
        # Scale with the binning
        binning = dic["dBinningX"]

        # microns_per_pixel = np.divide(microns, pixels)
        info["um_per_pixel"] = 0.325 * binning

    info["size_x"] = images.shape[0]
    info["size_y"] = images.shape[1]

    if "CustomData|MyoPacer" in metadata:
        pacing = metadata["CustomData|MyoPacer"]
    elif "CustomData|PXI1Slot2/port0/line6" in metadata:
        pacing = metadata["CustomData|PXI1Slot2/port0/line6"]
    else:
        pacing = np.zeros(len(time_stamps))

    return mps_data(
        frames=images,
        info=info,
        time_stamps=time_stamps,
        metadata=metadata,
        pacing=pacing,
    )


def load_czi(fname):

    with czifile.CziFile(fname) as f:
        images = f.asarray().squeeze()
        metadata = czifile.elem2dict(f.metadata)

        time_stamps = None
        pacing_triggers = []
        for a in f.attachment_directory:

            d = a.data_segment().data()

            if isinstance(d, czifile.TimeStamps):
                time_stamps = np.array(d.time_stamps)

            if isinstance(d, czifile.EventList):
                for event in d.events:
                    pacing_triggers.append(event.time)

    time_stamps -= time_stamps[0]

    # Make sure time is in milliseconds
    time_unit = analysis.time_unit(time_stamps)
    if time_unit == "s":
        time_stamps *= 1000

    info = info_dictionary(time_stamps)
    info["um_per_pixel"] = (
        float(metadata["Metadata"]["Scaling"]["Items"]["Distance"]["Value"]) * 1e6
    )
    info["size_x"] = images.shape[0]
    info["size_y"] = images.shape[1]

    return mps_data(
        frames=np.swapaxes(images, 0, -1),
        time_stamps=time_stamps,
        pacing=np.zeros(len(time_stamps)),
        info=info,
        metadata=metadata,
    )


def load_zip(fname):

    try:
        import xmltodict
    except ImportError:
        msg = "Please install xmltodict - pip install xmltodict"
        raise ImportError(msg)

    with zipfile.ZipFile(fname, "r") as f:

        metaname = [g for g in f.namelist() if os.path.basename(g) == "meta.xml"][0]
        with f.open(metaname) as fid:
            metadata = xmltodict.parse(fid.read())

        dtype = np.dtype(metadata["OME"]["Image"]["Pixels"]["@Type"])
        size_x = int(metadata["OME"]["Image"]["Pixels"]["@SizeX"])
        size_y = int(metadata["OME"]["Image"]["Pixels"]["@SizeY"])

        time_stamps = np.array(
            [float(p["@DeltaT"]) for p in metadata["OME"]["Image"]["Pixels"]["Plane"]]
        )
        # Normalize so that they start at zero
        time_stamps -= time_stamps[0]
        info = info_dictionary(time_stamps)

        info.update(
            **dict(
                size_x=size_x,
                size_y=size_y,
                um_per_pixel=float(
                    metadata["OME"]["Image"]["Pixels"]["@PhysicalSizeY"]
                ),
            )
        )

        frames = np.zeros((size_x, size_y, info["num_frames"]), dtype=dtype)
        for filename in f.namelist():
            if filename == metaname:
                continue

            idx = (
                int(os.path.splitext(os.path.basename(filename))[0].split("_")[-1]) - 1
            )
            try:
                data = utils.loadmat(f.open(name=filename))
            except io.UnsupportedOperation:
                full_path = f.extract(filename)
                data = utils.loadmat(full_path)
                os.remove(full_path)

            frames[:, :, idx] = data["image"].T

        return mps_data(
            frames=frames,
            time_stamps=time_stamps,
            pacing=np.zeros(len(time_stamps)),
            info=info,
            metadata=metadata,
        )


def load_stk(fname):

    if not has_tifffile:
        raise ImportError(
            (
                "tifffile is not installed. Please install "
                "that if you want to load stk files. python -m"
                " pip install tiffile"
            )
        )

    with tifffile.TiffFile(fname) as f:
        metadata = f.stk_metadata
        frames = f.asarray()

    time_stamps = metadata["TimeCreated"]
    time_stamps -= time_stamps[0]

    info = info_dictionary(time_stamps)

    info.update(
        **dict(
            size_x=frames.shape[0] * metadata["XCalibration"],
            size_y=frames.shape[1] * metadata["YCalibration"],
            um_per_pixel=float(
                metadata["YCalibration"]  # Different in x- and y direction
            ),
        )
    )

    return mps_data(
        frames=np.swapaxes(frames, 0, -1),
        time_stamps=time_stamps,
        pacing=np.zeros(len(time_stamps)),
        info=info,
        metadata=metadata,
    )


def time2isoformat(s):
    import datetime

    date, time = s.split(" ")
    return datetime.datetime.fromisoformat(f"{date[:4]}-{date[4:6]}-{date[6:]}T{time}")


def load_tiff_timestamps(f):
    import tifffile

    num_pages = len(f.pages.pages)
    get_page_timstamp = lambda i: time2isoformat(
        tifffile.tifffile.metaseries_description_metadata(f.pages.get(i).description)[
            "PlaneInfo"
        ]["acquisition-time-local"]
    )
    return list(
        map(
            lambda x: (x - get_page_timstamp(0)).total_seconds(),
            map(get_page_timstamp, range(num_pages)),
        )
    )


def load_tiff(fname):

    if not has_tifffile:
        raise ImportError(
            (
                "tifffile is not installed. Please install "
                "that if you want to load stk files. python -m"
                " pip install tiffile"
            )
        )

    with tifffile.TiffFile(fname) as f:
        metadata = f.metaseries_metadata
        frames = f.asarray()
        time_stamps = np.multiply(load_tiff_timestamps(f), 1000)

    info = info_dictionary(time_stamps)

    info.update(
        **dict(
            size_x=metadata["PlaneInfo"]["pixel-size-x"],
            size_y=metadata["PlaneInfo"]["pixel-size-y"],
            um_per_pixel=float(
                metadata["PlaneInfo"][
                    "spatial-calibration-y"
                ]  # Different in x- and y direction
            ),
        )
    )

    return mps_data(
        frames=np.swapaxes(frames, 0, -1),
        time_stamps=time_stamps,
        pacing=np.zeros(len(time_stamps)),
        info=info,
        metadata=metadata,
    )


def load_file(fname, ext):

    if ext == ".czi":
        data = load_czi(fname)

    elif ext == ".nd2":
        data = load_nd2(fname)

    elif ext == ".zip":
        data = load_zip(fname)

    elif ext == ".stk":
        data = load_stk(fname)

    elif ext in [".tiff", ".tif"]:
        data = load_tiff(fname)

    elif ext == ".npy":
        data = mps_data(**np.load(fname, allow_pickle=True).item())

    else:
        raise IOError("Unknown file extension {}".format(ext))

    return data


class MPS(object):
    """
    Create instant of class for analysing MPS data

    Arguments
    --------
    fname : str
        Path to the file
    parameters : dict
        Additional parameters, see
        `MPS.default_parameters`



    """

    def __init__(self, fname="", parameters=None, data=None):

        self.parameters = MPS.default_parameters()
        if parameters is not None:
            self.parameters.update(**parameters)
        logger.setLevel(self.parameters["log_level"])

        if fname == "":
            self._fname = ""
            self._ext = ""
            if data is None:
                raise ValueError("Please provide data of filename")
        else:
            # Check extension
            _, ext = os.path.splitext(fname)

            folder = os.path.abspath(os.path.dirname(fname))
            # Find the correct extension
            for f in os.listdir(folder):
                if os.path.isfile(os.path.join(folder, f)):
                    name_, ext_ = os.path.splitext(f)
                    if name_ == fname:
                        ext = ext_

            self._fname = os.path.abspath(fname)
            self._ext = ext

            data = load_file(self._fname, self._ext)
        self._unpack(data)

    @classmethod
    def from_dict(cls, **kwargs):
        required_fields = ["frames", "time_stamps", "info"]
        data = {}
        for field in required_fields:
            assert field in kwargs, f"Field {field} is not provided"
            data[field] = kwargs[field]
        data["pacing"] = kwargs.get("pacing", np.zeros_like(data["time_stamps"]))
        data["metadata"] = kwargs.get("metadata", {})
        return cls(data=mps_data(**data))

    def __getstate__(self):
        # Copy the object's state from self.__dict__ which contains
        # all our instance attributes. Always use the dict.copy()
        # method to avoid modifying the original state.
        state = deepcopy(self.__dict__)
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)

    @staticmethod
    def default_parameters():
        """
        Returns
        -------
        parameters: dict
            A dictionary with the follings key-value pairs:

            background_polynomial_order : int
                Order of polynomial applying in the background
                correction algoritm
            averaging_type : str
                Type of averaging. Possible inputs are
                ["spatial", "temporal", "global"] (Default: "global")
            alpha : str
                Take mean over value larger than this percentage value
                (Default : 1.0, i.e take mean of all values)
            local : array
                List of indice to take take average over
                [startx, endx, starty, endy]
            log_level : int
                Level of logging: Default 20 (INFO)
        """
        return dict(
            background_polynomial_order=2,
            log_level=logging.INFO,
            alpha=1.0,
            averaging_type="temporal",
            local=[],
        )

    def __repr__(self):
        return ("{self.__class__.__name__}" "({self._fname})").format(self=self)

    @property
    def name(self):
        return os.path.basename(self._fname)

    @property
    def frames(self):
        return self.data.frames

    @property
    def info(self):
        return self.data.info

    def _unpack(self, data):
        self.data = data
        self.time_stamps = data.time_stamps

        # Set dt, num_frames and time_unit
        for k, v in data.info.items():
            setattr(self, k, v)

    @property
    def framerate(self):
        """
        Return number of frames per second
        """

        factor = 1e-3 if self.info["time_unit"] == "ms" else 1.0
        time_increment = self.info["dt"] * factor
        fps = round(1.0 / time_increment)
        return fps

    @property
    def pacing(self):
        return self.data.pacing

    @property
    def metadata(self):
        return self.data.metadata

    @property
    def original_time_stamps(self):
        return self.data.time_stamps
