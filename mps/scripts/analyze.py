"""
Analyze flourecense data
"""
import datetime
import json
import logging
import os
import shutil
import time
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, Optional

from ..analysis import analyze_mps_func
from ..load import MPS, valid_extensions
from ..utils import json_serial

logger = logging.getLogger(__name__)


def dump_settings(outdir, kwargs):
    from .. import __version__

    kwargs["time"] = datetime.datetime.now()
    kwargs["mps_version"] = __version__
    kwargs["full_path"] = os.path.abspath(kwargs["path"])
    with open(os.path.join(outdir, "settings.json"), "w") as f:
        json.dump(kwargs, f, indent=4, default=json_serial)


def run_folder(**kwargs):

    path = Path(kwargs.get("path"))
    if not path.exists():
        raise IOError(f"Folder {path} does not exist")

    for root, dirs, files in os.walk(path):
        for f in files:

            if Path(f).suffix not in valid_extensions:
                continue

            # Exclude Brightfield
            if "BF" in f:
                continue

            fkwargs = deepcopy(kwargs)
            fkwargs["path"] = Path(root).joinpath(f)
            fkwargs["outdir"] = Path(root).joinpath(Path(f).stem)
            run_file(**fkwargs)


def update_kwargs(kwargs: Dict[str, Any], outdir: Path):
    if not outdir.is_dir():
        return

    if (
        kwargs.get("reuse_settings", False)
        and outdir.joinpath("settings.json").is_file()
    ):
        logger.debug("Reuse settings")
        with open(outdir.joinpath("settings.json"), "r") as f:
            old_kwargs = json.load(f)

        kwargs.update(**old_kwargs)


def check_overwrite(kwargs: Dict[str, Any], outdir: Path):

    if not outdir.is_dir():
        return

    if kwargs.get("overwrite", True):
        try:
            shutil.rmtree(outdir)
        except Exception as ex:
            logger.warning(ex, exc_info=True)

    else:
        # Check if we can find the version number of the
        # software used to analyze the data
        try:
            with open(outdir.joinpath("settings.json"), "r") as f:
                settings = json.load(f)
            version = settings["mps_version"]
        except (FileNotFoundError, KeyError):
            old_dir = outdir.joinpath("old")
        else:
            old_dir = outdir.joinpath(f"old_version_{version}")

        old_dir.mkdir(exist_ok=True, parents=True)
        for p in outdir.iterdir():

            if Path(p).name == old_dir.name:
                # We cannot move the old_dir into itself
                continue
            if p.name.startswith("._"):
                try:
                    p.unlink()
                except FileNotFoundError:
                    pass

                continue
            shutil.move(p.as_posix(), old_dir.joinpath(p.name).as_posix())


def run_file(**kwargs):

    path = Path(kwargs.get("path"))
    if not path.exists():
        raise IOError(f"File {path} does not exist")

    outdir = kwargs.get("outdir")
    if outdir is None:
        outdir = path.parent.joinpath(path.stem)
    outdir = Path(outdir)
    kwargs["outdir"] = outdir.absolute().as_posix()

    if outdir.is_dir():
        update_kwargs(kwargs, outdir)
        check_overwrite(kwargs, outdir)

    start = time.time()

    mps_data = MPS(path)
    analyze_mps_func(mps_data, **kwargs)
    dump_settings(outdir, kwargs)

    end = time.time()
    logger.info(
        (
            f"Finished analyzing MPS data. Data stored in {outdir}. "
            f"\nTotal elapsed time: {end - start} seconds"
        )
    )


def main(
    path: str,
    outdir: Optional[str] = None,
    plot: bool = True,
    filter_signal: bool = True,
    alpha: float = 1.0,
    ead_prom: float = 0.04,
    ead_sigma: float = 3.0,
    std_ex_factor: float = 1.0,
    spike_duration: float = 0.0,
    chopping_threshold_factor: float = 0.3,
    chopping_extend_front: Optional[int] = None,
    chopping_extend_end: Optional[int] = None,
    chopping_min_window: Optional[int] = None,
    chopping_max_window: Optional[int] = None,
    ignore_pacing: bool = False,
    reuse_settings: bool = False,
    overwrite: bool = True,
    verbose: bool = False,
):

    level = logging.DEBUG if verbose else logging.INFO
    logger.setLevel(level)

    logger.info("Run analysis script")

    filepath = Path(path)

    if not filepath.exists():
        raise IOError(f"Path {filepath} does not exist")

    kwargs = dict(
        path=path,
        outdir=outdir,
        plot=plot,
        filter_signal=filter_signal,
        alpha=alpha,
        ead_prom=ead_prom,
        ead_sigma=ead_sigma,
        std_ex_factor=std_ex_factor,
        spike_duration=spike_duration,
        chopping_threshold_factor=chopping_threshold_factor,
        chopping_extend_front=chopping_extend_front,
        chopping_extend_end=chopping_extend_end,
        chopping_min_window=chopping_min_window,
        chopping_max_window=chopping_max_window,
        ignore_pacing=ignore_pacing,
        reuse_settings=reuse_settings,
        overwrite=overwrite,
    )

    if filepath.is_dir():
        run_folder(**kwargs)
    else:
        run_file(**kwargs)
