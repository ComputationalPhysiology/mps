# -*- coding: utf-8 -*-
# flake8: noqa
# Copyright (c) 2013-2016 Christian C. Sachs.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
This is an Nikon ND2 (NIS Elements file format) reader. Information about the ND2 file format
was gained from various other readers, and from ND2 files.
As it uses memory-mapping to open the file, it is most often only useful
for 64-bit Python (otherwise the file size is limited to 2 GiB).
There are some heuristics for parsing files, so that even broken files (i.e. due to a crashed NIS instance)
have a possibility to open.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import mmap
import warnings
from functools import partial
from xml.etree import ElementTree

import numpy as np


class ND2File(object):
    CHUNK_MAP_SIGNATURE = b"ND2 CHUNK MAP SIGNATURE 0000001"
    CHUNK_FILE_SIGNATURE = b"ND2 FILE SIGNATURE CHUNK NAME01"
    CHUNK_MAGIC = 0xABECEDA

    CHUNK_MAGIC_BYTES = np.array([CHUNK_MAGIC], dtype=np.uint32).tobytes()

    allow_file_scanning = True

    def __init__(self, filename):
        self.handle = open(filename, "rb")
        self.mem = mmap.mmap(self.handle.fileno(), 0, access=mmap.ACCESS_READ)
        self.size = len(self.mem)

        self.chunkmap = {}
        self.images = {}

        self.metadata = {}

        self.depths = {8: np.uint8, 16: np.uint16}

        process_map = [
            # (self.xml, [
            # b'ImageTextInfo', b'ImageMetadataSeq|0', b'ImageMetadata', b'ImageCalibration|0', b'ImageAttributes'
            # ]),
            (
                partial(self.binaryblock, dtype=np.double),
                [
                    "CustomData|X",
                    "CustomData|Y",
                    "CustomData|Z",
                    "CustomData|AcqTimesCache",
                    "CustomData|PXI1Slot2/ao0",
                    "CustomData|Camera_ExposureTime1",
                    "CustomData|Camera_ExposureTime2",
                    "CustomData|Camera_Temp1",
                    "CustomData|Camera_Temp2",
                ],
            ),
            (
                partial(self.binaryblock, dtype=np.uint32),
                [
                    "CustomData|MyoPacer",
                    "CustomData|PFS_STATUS",
                    "CustomData|PFS_OFFSET",
                    "CustomData|PXI1Slot2/port0/line6",
                    "CustomData|AcqFramesCache",
                ],
            ),
            (
                self.lv,
                [
                    "ImageTextInfoLV",
                    "ImageMetadataSeqLV|0",
                    "ImageMetadataLV",
                    "CustomData|Z1",
                    "CustomData|LimApi_CustomData",
                    # 'CustomData|JobDefinitionV1_0'
                    "CustomData|JobRunGUIDV1_0",
                    "CustomData|JobDefinitionDescV1_0",
                    "CustomData|TimeMeasurementData_V1_0",
                    "CustomData|RoiMetadata_v1",
                    "CustomData|GUIDStore",
                    "ImageMetadataLV",
                    "ImageCalibrationLV|0",
                    "ImageAttributesLV",
                    "ImageEventsLV",
                ],
            ),
        ]

        process_map = {key: fun for fun, keys in process_map for key in keys}

        self.opportunistic = False

        if self.chunk_name(0) != self.CHUNK_FILE_SIGNATURE:
            raise RuntimeError("This is probably not an ND2 file.")

        if self.mem[-40:-9] != self.CHUNK_MAP_SIGNATURE:
            if not self.allow_file_scanning:
                raise RuntimeError(
                    "Did not find chunk map signature at the end of file. "
                    "This file might not be an ND2 file, or "
                    "the file might be damaged (software crash?) and "
                    "cannot directly be opened with this library.",
                )
            else:
                warnings.warn(
                    "The file is damaged and has no chunk map, will use "
                    "heuristics to calculate image positions.",
                )
                readahead, readback = 25, 25
                self.opportunistic = True
                mini_chunk_map = self.opportunistic_chunk_scanner(count=readahead)
                prefix = b"ImageDataSeq|"

                def parse_images(inner_mini_chunk_map):
                    prefix_len = len(prefix)
                    return np.array(
                        [
                            [
                                int(name_[prefix_len:]),
                                chunk_begin_,
                                chunk_end_,
                                chunk_position_,
                            ]
                            for name_, chunk_begin_, chunk_end_, chunk_position_ in inner_mini_chunk_map
                            if name_.startswith(prefix)
                        ],
                    )

                images = parse_images(mini_chunk_map)
                slope = int((np.diff(images[:, 3]) / np.diff(images[:, 0])).mean())
                intercept = images[0, 3]

                self.opp_slope = slope
                self.opp_intercept = intercept

                self.opp_shifts = []

                # now find the end
                max_images = self.size // self.opp_slope

                def opportunistic_image_position(num):
                    return self.opp_intercept + num * self.opp_slope

                mini_chunk_map += self.opportunistic_chunk_scanner(
                    opportunistic_image_position(max_images - readback),
                    2 * readback,
                )

                images = parse_images(mini_chunk_map)

                for i in range(images[:, 0].min(), images[:, 0].max()):
                    self.images[i] = opportunistic_image_position(i)

                for (name, chunk_begin, chunk_end, chunk_position) in mini_chunk_map:
                    if name.startswith(prefix):
                        continue
                    name = name.decode("ascii")

                    self.chunkmap[name] = chunk_position
                    if name in process_map:
                        self.metadata[name] = process_map[name](chunk_begin, chunk_end)

        if not self.opportunistic:

            map_position = self._sfp(-8, np.uint64)
            pos, epos = self.chunk_position(map_position)

            while pos < epos:
                nameendpos = self.mem.find(b"!", pos, epos)
                name = self.mem[pos:nameendpos].decode("ascii")

                if nameendpos + 1 + 16 > epos:
                    break
                position = self._sfp(nameendpos + 1, np.uint64)

                nsplit = name.split("|")

                if nsplit[0] == "ImageDataSeq":
                    self.images[int(nsplit[1])] = position
                else:

                    self.chunkmap[name] = position
                    if name in process_map:
                        self.metadata[name] = process_map[name](
                            *self.chunk_position(position)
                        )

                pos = nameendpos + 1 + 16

        self.imagecount = len(self.images)

        if "ImageAttributes" in self.metadata:
            attributes = self.metadata["ImageAttributes"]
        elif "ImageAttributesLV" in self.metadata:
            attributes = self.metadata["ImageAttributesLV"]["SLxImageAttributes"]
        else:
            raise RuntimeError("No ImageAttributes")

        if "ImageAttributes" in self.metadata:
            calib = self.metadata["ImageCalibration|0"]
        elif "ImageAttributesLV" in self.metadata:
            calib = self.metadata["ImageCalibrationLV|0"]["SLxCalibration"]
        else:
            raise RuntimeError("No calibration")

        self.width = attributes["uiWidth"]
        self.width_bytes = attributes["uiWidthBytes"]
        self.height = attributes["uiHeight"]
        self.calibration = calib["dCalibration"]
        self.bpc = attributes["uiBpcInMemory"]
        self.bypc = self.bpc // 8
        self.bpcsig = attributes["uiBpcSignificant"]
        self.channels = attributes["uiComp"]
        self.array_stride = attributes["uiWidthBytes"] - (
            self.bypc * attributes["uiWidth"] * self.channels
        )

    def chunk_position(self, begin_position):
        magic, shift = self._dfp(begin_position, np.uint32, count=2)
        shift = int(shift)
        thelen = self._sfp(begin_position + 8, np.uint64)
        if magic != self.CHUNK_MAGIC:
            raise RuntimeError("Error")
        real_bpos = begin_position + 16 + shift
        return real_bpos, real_bpos + thelen

    def chunk_name(self, begin_position):
        return self._readdelimitedstring(begin_position + 16, b"!")

    def opportunistic_chunk_scanner(self, begin_position=0, count=-1):
        position = begin_position
        result = []
        while position < self.size:
            new_pos = self.mem.find(self.CHUNK_MAGIC_BYTES, position)
            chunk_begin, chunk_end = self.chunk_position(new_pos)
            result.append((self.chunk_name(new_pos), chunk_begin, chunk_end, new_pos))
            position = chunk_end
            count -= 1
            if count == 0:
                break
        return result

    def _readdelimitedstring(self, position, delimiter):
        return self.mem[position : self.mem.find(delimiter, position)]

    def _readcstring(self, position):
        return self._readdelimitedstring(position, b"\0")

    def _dfp(self, position, dtype, count=1):
        if position < 0:
            position += self.size
        return np.ndarray(shape=(count,), dtype=dtype, buffer=self.mem, offset=position)

    def _sdfp(self, position, dtype, count=1):
        return map(np.asscalar, self._dfp(position, dtype, count))

    def _sfp(self, position, dtype):
        return np.asscalar(self._dfp(position, dtype, count=1)[0])

    def _readstring(self, position, length):
        return self.mem[position : position + 2 * length].decode(
            "utf-16",
        )  # .encode('utf-8')

    def xml(self, bpos, epos):
        root = ElementTree.fromstring(self.mem[bpos:epos])
        types = {
            "double": float,
            "lx_uint32": int,
            "bool": bool,
            "CLxStringW": lambda s: s.encode("utf-8"),
        }

        def recurse(node):
            if node.attrib["runtype"] in types:
                return types[node.attrib["runtype"]](node.attrib["value"])
            else:
                return dict(map(lambda c: (c.tag, recurse(c)), node))

        return recurse(root[0])

    # refactor this beast!
    def lv(self, pos, epos, num=1):
        tmap = [None, np.uint8, np.uint32, np.uint32, None, np.uint64, np.double]
        result = {}
        while pos < epos and num > 0:
            thetype, thecount = self._sdfp(pos, np.uint8, count=2)
            pos += 2

            thename = self._readstring(pos, thecount - 1)
            pos += 2 * thecount
            value = None
            action = tmap[thetype] if thetype < len(tmap) else None
            if type(action) == type:
                value = self._sfp(pos, action)
                pos += np.dtype(action).itemsize
            elif thetype == 8:
                e = self.mem.find(b"\0\0", pos, epos)
                if (e - pos) % 2:
                    e += 1
                value = self._readstring(pos, (e - pos) // 2)
                pos = e + 2
            elif thetype == 9:
                num = self._sfp(pos, np.uint64)
                pos += 8
                value = self._dfp(pos, np.uint8, count=num).tolist()
                pos += num
            elif thetype == 11:
                t = self._sfp(pos, np.uint32)
                pos += 4
                l_ = self._sfp(pos, np.uint64) - ((thecount + 1) * 2 + 12)
                pos += 8
                value = self.lv(pos, pos + l_, t)
                pos += l_ + t * 8
            else:
                pass
                # print("WARNING", thetype)

            if thename in result:
                if type(result[thename]) == list:
                    result[thename] += [value]
                else:
                    result[thename] = [result[thename], value]
            else:
                result[thename] = value
            num -= 1
        return result

    def binaryblock(self, bpos, epos, dtype=None):
        return np.ndarray(
            shape=((epos - bpos) // np.dtype(dtype).itemsize),
            dtype=dtype,
            buffer=self.mem,
            offset=bpos,
        )

    def image(self, num):
        if self.opportunistic:
            found = False
            bpos, epos = 0, 0

            for shift in [0] + self.opp_shifts:
                try:
                    bpos, epos = self.chunk_position(self.images[num] + shift)
                    found = True
                    if shift != 0:
                        self.images[num] += shift
                except RuntimeError:
                    continue

            if not found:
                warnings.warn(
                    "Searching for the image, this will lead to degraded performance.",
                )
                look_before = 1 * 1024 * 1024
                pos = max(0, self.images[num] - look_before)
                search = ("ImageDataSeq|" + str(num)).encode("ascii")
                pos = self.mem.find(search, pos) - 16
                warnings.warn(
                    "Was looking at %d found it at %d ... %d bytes away."
                    % (self.images[num], pos, self.images[num] - pos),
                )
                self.opp_shifts.append(pos - self.images[num])
                self.images[num] = pos
                bpos, epos = self.chunk_position(pos)
        else:
            bpos, epos = self.chunk_position(self.images[num])

        if self.array_stride == 0:
            return np.ndarray(
                shape=(self.height, self.width, self.channels),
                dtype=self.depths[self.bpc],
                buffer=self.mem,
                offset=bpos + 8,
            )
        else:
            return np.ndarray(
                shape=(self.height, self.width, self.channels),
                dtype=self.depths[self.bpc],
                buffer=self.mem,
                offset=bpos + 8,
                strides=(
                    self.array_stride + self.width * self.bypc * self.channels,
                    self.channels * self.bypc,
                    self.bypc,
                ),
            )

    def image_singlechannel(self, num, channel=0):
        return self.image(num)[:, :, channel]

    def get_time(self, num):
        return self.metadata["CustomData|AcqTimesCache"][num] / 1000.0

    def close(self):
        self.handle.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


class ND2MultiFile(ND2File):
    def __init__(self, filename):
        super(ND2MultiFile, self).__init__(filename)

        # so many heuristics ... maybe just search for uLoopPars recursively?
        # ... yep.

        def find_by_key(tree, needle):
            class FoundIt(Exception):
                def __init__(self, result):
                    self.result = result

            try:

                def _recurse(item):
                    if isinstance(item, dict):
                        for k, v in sorted(item.items()):
                            if k == needle:
                                raise FoundIt(v)
                            if isinstance(item, dict) or isinstance(item, list):
                                _recurse(v)
                    elif isinstance(item, list):
                        for v in item:
                            _recurse(v)
                    else:
                        pass

                _recurse(tree)
            except FoundIt as f:
                return f.result
            return None

        def find_all_by_key(tree, needle):
            results = []

            def _recurse(item):
                if isinstance(item, dict):
                    for k, v in sorted(item.items()):
                        if k == needle:
                            results.append(v)
                        if isinstance(item, dict) or isinstance(item, list):
                            _recurse(v)
                elif isinstance(item, list):
                    for v in item:
                        _recurse(v)
                else:
                    pass

            _recurse(tree)
            return results

        self.experiment = find_by_key(self.metadata, "SLxExperiment")

        if find_by_key(self.experiment, "ppNextLevelEx"):
            self.experiment = find_by_key(self.experiment, "ppNextLevelEx")[""]

        self._points = next(
            iter(
                intmd
                for intmd in find_all_by_key(self.experiment, "Points")
                if find_by_key(intmd, "dPosName") is not None
            ),
        )[""]

        try:
            self.multipoints_of_experiment = [
                {
                    "x": p["dPosX"],
                    "y": p["dPosY"],
                    "z": p["dPosZ"],
                    "pfs": p["dPFSOffset"],
                    "name": p["dPosName"],
                }
                for p in self._points
            ]

            valid_points = find_by_key(self.experiment, "pItemValid")

            if valid_points:
                for i, n in enumerate(valid_points):
                    self.multipoints_of_experiment[i]["valid"] = n == 1
            else:
                for i, n in enumerate(self.multipoints_of_experiment):
                    self.multipoints_of_experiment[i]["valid"] = True

        except KeyError:
            self.multipoints_of_experiment = [{"valid": True}]

        self.multipoints = [p for p in self.multipoints_of_experiment if p["valid"]]
        self.multipointcount = len(self.multipoints)

        self.timepointcount = self.imagecount // self.multipointcount

        # Channel 'colors'
        try:
            tmp = self.metadata["ImageMetadataSeqLV|0"]["SLxPictureMetadata"][
                "sPicturePlanes"
            ]["sPlane"]
        except KeyError:
            # new file
            tmp = self.metadata["ImageMetadataSeqLV|0"]["SLxPictureMetadata"][
                "sPicturePlanes"
            ]["sPlaneNew"]

        tmp = [tmp[n] for n in sorted(tmp.keys())]
        tmp = [c["uiColor"] for c in tmp]

        self.channelcolors = [
            [(c & 0xFF0000) >> 16, (c & 0xFF00) >> 8, (c & 0xFF)][::-1] for c in tmp
        ]

        # probably easier to associate with certain channel parameters

        self.heuristic_pcm = (
            [n for n, c in enumerate(self.channelcolors) if c == [255, 255, 255]]
            or [None]
        )[0]
        self.heuristic_fluorescence = (
            [n for n, c in enumerate(self.channelcolors) if c != [255, 255, 255]]
            or [None]
        )[0]
        self.heuristic_fluorescences = [
            n for n, c in enumerate(self.channelcolors) if c != [255, 255, 255]
        ]

    def calc_num(self, multipoint=0, timepoint=0):
        return timepoint * self.multipointcount + multipoint

    def image(self, multipoint=0, timepoint=0):
        return super(ND2MultiFile, self).image(
            self.calc_num(multipoint=multipoint, timepoint=timepoint),
        )

    def image_singlechannel(self, multipoint=0, timepoint=0, channel=0):
        return self.image(multipoint=multipoint, timepoint=timepoint)[:, :, channel]
