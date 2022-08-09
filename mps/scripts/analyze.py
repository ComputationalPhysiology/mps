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
from typing import Any
from typing import Dict
from typing import Optional

from ..analysis import analyze_mps_func
from ..load import MPS
from ..load import valid_extensions
from ..utils import json_serial

logger = logging.getLogger(__name__)


def dump_settings(outdir, kwargs):
    from .. import __version__

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True, parents=True)

    kwargs["time"] = datetime.datetime.now()
    kwargs["mps_version"] = __version__
    kwargs["full_path"] = Path(kwargs["path"]).absolute().as_posix()
    with open(outdir.joinpath("settings.json"), "w") as f:
        json.dump(kwargs, f, indent=4, default=json_serial)


def run_folder(**kwargs):

    path = Path(kwargs.get("path"))
    if not path.is_dir():
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
    if not path.is_file():
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
            f"Finished analyzing MPS data. Data stored folder '{outdir}'. "
            f"\nTotal elapsed time: {end - start:.2f} seconds"
        ),
    )


def main(
    path: str,
    outdir: Optional[str] = None,
    plot: bool = True,
    filter_signal: bool = False,
    ead_prom: float = 0.04,
    ead_sigma: float = 3.0,
    std_ex_factor: float = 1.0,
    spike_duration: float = 0.0,
    threshold_factor: float = 0.3,
    extend_front: Optional[int] = None,
    extend_end: Optional[int] = None,
    min_window: float = 50,
    max_window: float = 2000,
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
        ead_prom=ead_prom,
        ead_sigma=ead_sigma,
        std_ex_factor=std_ex_factor,
        spike_duration=spike_duration,
        threshold_factor=threshold_factor,
        extend_front=extend_front,
        extend_end=extend_end,
        min_window=min_window,
        max_window=max_window,
        ignore_pacing=ignore_pacing,
        reuse_settings=reuse_settings,
        overwrite=overwrite,
    )

    if filepath.is_dir():
        run_folder(**kwargs)
    else:
        run_file(**kwargs)
