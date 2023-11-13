from . import analyze
from . import mps2mp4
from . import phase_plot
from . import split_pacing
from . import summary

_loggers = [
    split_pacing.logger,  # type:ignore
    analyze.logger,  # type:ignore
    summary.logger,  # type:ignore
    mps2mp4.logger,  # type:ignore
    phase_plot.logger,  # type:ignore
]

__all__ = ["split_pacing", "analyze", "summary", "mps2mp4", "phase_plot"]
