"""
Make a phase plot with voltage on the x-axis
and calcium on the y-axis."""
import logging
from pathlib import Path
from typing import Optional

import ap_features as apf

from .. import analysis
from ..load import MPS
from ..plotter import phase_plots

logger = logging.getLogger(__name__)


def main(voltage: str, calcium: str, outfile: Optional[str] = None):

    logger.setLevel(logging.INFO)
    calcium_path = Path(calcium)
    voltage_path = Path(voltage)

    if not calcium_path.is_file():
        raise IOError(f"Invald path for calcium file: {calcium}")

    if not voltage_path.is_file():
        raise IOError(f"Invald path for voltage file: {voltage}")

    if outfile is None:
        outpath = voltage_path.parent.joinpath(
            f"{voltage_path.stem}_{calcium_path.stem}",
        )
    else:
        outpath = Path(outfile)

    logger.info(
        f"Create phase plot of voltage at {voltage_path} and calcium at {calcium_path}.",
    )

    outpath.parent.mkdir(exist_ok=True, parents=True)

    voltage_file = MPS(voltage_path)
    calcium_file = MPS(calcium_path)

    voltage_data = analysis.analyze_mps_func(voltage_file, plot=False)
    calcium_data = analysis.analyze_mps_func(calcium_file, plot=False)

    voltage_trace = voltage_data["chopped_data"]["trace_1std"]
    calcium_trace = calcium_data["chopped_data"]["trace_1std"]

    phase_plots(
        apf.utils.normalize_signal(voltage_trace),
        apf.utils.normalize_signal(calcium_trace),
        outpath,
    )

    logger.info(f"Saved to {outpath}")
