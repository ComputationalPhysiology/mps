"""
Create a summary pdf of all files in the a directory.
"""
import logging
import shutil
from pathlib import Path
from typing import List

try:
    import matplotlib.pyplot as plt

    has_mpl = True
except ImportError:
    has_mpl = False


import ap_features as apf
import numpy as np

from .. import analysis, utils
from ..load import MPS, valid_extensions

logger = logging.getLogger(__name__)


def get_data(path: Path, ignore_pacing: bool = False):
    mps_data = MPS(path)
    try:
        data = analysis.analyze_mps_func(
            mps_data=mps_data,
            ignore_pacing=ignore_pacing,
            plot=False,
            outdir=None,
        )
    except Exception as ex:
        logger.warning(ex, exc_info=True)
        return None
    else:
        return dict(
            average=data["unchopped_data"]["original_trace"],
            time_stamps=data["unchopped_data"]["original_times"],
            pacing=data["unchopped_data"]["original_pacing"],
            apd30s=data["all_features"]["features_beats"]["apd30"],
            apd80s=data["all_features"]["features_beats"]["apd80"],
            apd90s=data["all_features"]["features_beats"]["apd90"],
            capd30s=data["all_features"]["features_beats"]["capd30"],
            capd80s=data["all_features"]["features_beats"]["capd80"],
            capd90s=data["all_features"]["features_beats"]["capd90"],
            slope_APD80=data["all_features"]["features_beats"]["slope_APD80"],
            slope_cAPD80=data["all_features"]["features_beats"]["slope_cAPD80"],
            triangulation=data["all_features"]["features_beats"]["triangulation"],
            tau75s=data["all_features"]["features_beats"]["tau75"],
            upstroke80s=data["all_features"]["features_beats"]["upstroke80"],
            nbeats=data["all_features"]["features_beats"]["num_beats"],
            ttp=data["all_features"]["features_beats"]["ttp"],
            freqs=data["all_features"]["features_beats"]["beating_frequencies"],
        )


def plot(
    files: List[Path],
    folder: Path,
    filename: str = "mps_summary",
    ignore_pacing: bool = False,
):

    fig, ax = plt.subplots(len(files), 2, constrained_layout=True)
    if len(files) == 1:
        ax = np.expand_dims(ax, axis=0)

    xl_datas = []

    figdir = folder.joinpath("tmpfigures")
    figdir.mkdir(exist_ok=True, parents=True)
    image_paths = []

    for i, path in enumerate(files):
        fig, ax = plt.subplots(1, 2)
        ax[0].set_title(path.name.replace("_", r"\_"))
        data = get_data(path, ignore_pacing)
        if data is None:
            continue

        ax[0].plot(
            data["time_stamps"],
            apf.utils.normalize_signal(data["average"]),
            color="b",
        )
        ax[0].plot(
            data["time_stamps"],
            apf.utils.normalize_signal(data["pacing"]),
            label="Pacing ignored!",
            color="r",
        )
        if ignore_pacing:
            ax[0].legend()

        max_std = 20
        apd_keys = ["apd30s", "apd80s", "apd90s", "capd30s", "capd80s", "capd90s"]
        nremoved = {k: 0 for k in apd_keys}
        for k in nremoved.keys():
            while np.std(data[k]) > max_std:
                # Find the value that is most distant from the mean and remove it
                idx = np.argmax(np.abs(np.subtract(data[k], np.mean(data[k]))))
                nremoved[k] += 1
                data[k] = np.delete(data[k], idx)

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
                    float(np.mean(data[key])),
                    float(np.std(data[key])),
                )
            s = "{}: {}".format(label, s1)
            ax[1].text(0.0, pos[j], s)

        ax[1].axis("off")

        figpath = figdir.joinpath(f"image_{i}.png")
        fig.savefig(figpath)
        image_paths.append(figpath)

        xl_datas.append(
            [
                path.name,
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
            ],
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

    utils.to_csv(xl_datas, folder.joinpath(filename), header)

    from PIL import Image

    images = [Image.open(x) for x in image_paths]
    widths, heights = zip(*(i.size for i in images))

    max_width = max(widths)
    total_height = sum(heights)

    new_im = Image.new("RGB", (max_width, total_height))

    y_offset = 0
    for im in images:
        new_im.paste(im, (0, y_offset))
        y_offset += im.size[1]

    new_im.save(folder.joinpath(f"{filename}.pdf"))
    shutil.rmtree(figdir)


def main(
    folder: str,
    filename: str = "mps_summary",
    ignore_pacing: bool = False,
    silent: bool = False,
    include_npy: bool = False,
):

    level = logging.WARNING if silent else logging.INFO
    logger.setLevel(level)

    if not has_mpl:
        logger.error(
            "Cannot run script without matplotlib. Please install that first - 'pip install matplotlib'",
        )
        return

    path = Path(folder)
    if not path.is_dir():
        raise IOError(f"Folder {path} does not exist")

    valid_ext = valid_extensions.copy()
    if include_npy:
        valid_ext.append(".npy")

    exclude_patterns = ["BF"]

    files = [
        f
        for f in path.iterdir()
        if (f.suffix in valid_ext and not any(p in f.name for p in exclude_patterns))
    ]

    if len(files) == 0:
        logger.warning(f"Could not find any files with in folder {path}")
        logger.info(f"Valid extensions {valid_ext}")
        logger.info(f"Exclude patterns {exclude_patterns}")
        return

    plot(files=files, folder=path, filename=filename, ignore_pacing=ignore_pacing)
