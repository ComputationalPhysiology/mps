"""
Run script on a folder with files and this will copy the files into
folders with the same pacing frequency
"""
import logging
import shutil
from pathlib import Path

from ..load import MPS
from ..load import valid_extensions

logger = logging.getLogger(__name__)


def move(src, dst):
    shutil.move(src, dst)


def copy(src, dst):
    shutil.copy2(src, dst)


def recursive_version(folder: Path, keep_original: bool):

    # func = copy if keep_original else move
    raise NotImplementedError("Recurive function not implemented yet")


def normal_version(folder: Path, keep_original: bool):

    func = copy if keep_original else move

    for f in folder.iterdir():
        logger.debug(f"Found file {f}")
        if f.suffix not in valid_extensions:
            continue

        logger.debug("Open file")
        data = MPS(f)
        # Round to 1 decimal point
        freq = round(data.pacing_frequency, 1)
        logger.debug(f"Frequency {freq}")

        dst = folder.joinpath(str(freq)).joinpath(f.name)
        dst.parent.mkdir(exist_ok=True, parents=True)

        func(f, dst)


def main(
    folder: str,
    recursive: bool = False,
    verbose: bool = False,
    keep_original: bool = True,
):

    level = logging.DEBUG if verbose else logging.INFO
    logger.setLevel(level)

    logger.info("Split files in to pacing folders")

    folder_path = Path(folder)
    if not folder_path.is_dir():
        raise ValueError(f"Path {folder_path} is not a directory")

    if recursive:
        recursive_version(folder_path, keep_original)
    else:
        normal_version(folder_path, keep_original)

    logger.info("Done spliting files in to pacing folders")
