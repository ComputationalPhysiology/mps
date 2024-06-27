import sys

if getattr(sys, "frozen", False):
    # We are running within pyinstaller
    import matplotlib

    matplotlib.use("Agg")

from importlib.metadata import metadata

import logging as _logging

from . import analysis, average, load, plotter, scripts, utils
from .load import MPS

log_level = _logging.INFO
_logging.basicConfig(format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
meta = metadata("cardiac-mps")
__version__ = meta["Version"]
__author__ = meta["Author-email"]
__license__ = meta["License"]
__email__ = meta["Author-email"]
__program_name__ = meta["Name"]


def set_log_level(level, logger=None):
    loggers = [logger]
    if logger is None:
        loggers = [
            utils.logger,
            load.logger,
            analysis.logger,
            average.logger,
            plotter.logger,
        ] + scripts._loggers

    for logger in loggers:
        logger.setLevel(level)
        for h in logger.handlers:
            h.setLevel(level)
    global log_level
    log_level = level


__all__ = [
    "load",
    "MPS",
    "utils",
    "analysis",
    "plotter",
    "average",
    "log_level",
    "scripts",
]
