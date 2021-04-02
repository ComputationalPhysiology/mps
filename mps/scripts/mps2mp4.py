"""
Create movie of data file
"""
import logging
from pathlib import Path
from typing import Optional

from .. import utils
from ..load import MPS

logger = logging.getLogger(__name__)


def main(
    path: str,
    outfile: Optional[str] = None,
    synch: bool = False,
):

    logger.setLevel(logging.INFO)
    file_path = Path(path)
    if not file_path.is_file():
        raise ValueError(f"Path {file_path} is not a file")

    logger.info("Convert file to video")

    if outfile is None:
        outfile = file_path.with_suffix(".mp4").as_posix()

    mps_data = MPS(file_path)

    if synch:
        idx = next(i for i, p in enumerate(mps_data.pacing) if p > 0)
    else:
        idx = 0

    utils.frames2mp4(mps_data.frames.T[idx:, :, :], outfile, mps_data.framerate)
    logger.info(f"Movie saved to {outfile}")
