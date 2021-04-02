#!/usr/bin/env python3
__author__ = "Henrik Finsberg (henriknf@simula.no), 2017--2021"
__maintainer__ = "Henrik Finsberg"
__email__ = "henriknf@simula.no"
__program_name__ = "MPS"
__license__ = """
c) 2001-2021 Simula Research Laboratory ALL RIGHTS RESERVED

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
from textwrap import dedent
from typing import Optional

import typer

from mps import __version__, scripts

app = typer.Typer()


def version_callback(show_version: bool):
    """Prints version information."""
    if show_version:
        typer.echo(f"{__program_name__} {__version__}")
        raise typer.Exit()


def license_callback(show_license: bool):
    """Prints license information."""
    if show_license:
        typer.echo(f"{__license__}")
        raise typer.Exit()


@app.callback()
def main(
    version: bool = typer.Option(
        None, "--version", callback=version_callback, is_eager=True, help="Show version"
    ),
    license: bool = typer.Option(
        None, "--license", callback=license_callback, is_eager=True, help="Show license"
    ),
):
    # Do other global stuff, handle other global options here
    return


@app.command(help=scripts.split_pacing.__doc__)
def split_pacing(
    folder: str = typer.Argument(..., help="The folder to be analyzed"),
    recursive: bool = typer.Option(False, help="Recursively go through all sufolders"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="More verbose"),
    keep_original: bool = typer.Option(
        True, help="If True, copy the files, otherwise move them."
    ),
):
    scripts.split_pacing.main(
        folder=folder, recursive=recursive, verbose=verbose, keep_original=keep_original
    )


@app.command(help=scripts.analyze.__doc__)
def analyze(
    path: str = typer.Argument(..., help="Path to file or folder to be analyzed"),
    outdir: Optional[str] = typer.Option(
        None,
        "--outdir",
        "-o",
        help=dedent(
            """
        Output directory for where you want to store the
        output. If not provided a folder with the same
        name as the basename of the input file in the
        currect directory will be used"""
        ),
    ),
    plot: bool = typer.Option(True, help="Plot data"),
    filter_signal: bool = typer.Option(
        True, help="Filter signal using a median filter"
    ),
    alpha: float = typer.Option(
        1.0,
        help=dedent(
            """
            When taking the average over the images include
            only values larger than the bound when all
            values are sorted. If alpha = 1.0, then all
            values will be used when taking the average. If
            alpha = 0.1 then only 10 percent of the pixels
            will be included."""
        ),
    ),
    ead_prom: float = typer.Option(
        0.04,
        help=dedent(
            """
            How prominent a peak should be in order to be
            characterized as an EAD. This value shold be
            between 0 and 1, with a greater value being more
            prominent"""
        ),
    ),
    ead_sigma: float = typer.Option(
        3.0,
        help=dedent(
            """
            Standard deviation in the gaussian smoothing
            kernal applied before estimating EADs."""
        ),
    ),
    std_ex_factor: float = typer.Option(
        1.0,
        help=dedent(
            """
            Exclude values outside this factor times 1
            standard deviation. Default:1.0, meaning exclude
            values outside of 1 std from the mean"""
        ),
    ),
    spike_duration: float = typer.Option(
        0.0,
        help=dedent(
            """
            Remove spikes from signal by setting this value
            greater than 0. This will locate the timing of
            pacing (if available) and delete the signal this
            amount after pacing starts. If 0 or no pacing is
            detected, nothing will be deleted."""
        ),
    ),
    chopping_threshold_factor: float = typer.Option(
        0.3,
        help=dedent(
            """
            Factor of where to synchronize data when
            chopping. Default = 0.3"""
        ),
    ),
    chopping_extend_front: Optional[int] = typer.Option(
        None,
        help=dedent(
            """
            How many milliseconds extra to extend signal at
            the beginning when chopping"""
        ),
    ),
    chopping_extend_end: Optional[int] = typer.Option(
        None,
        help=dedent(
            """
            How many milliseconds extra to extend signal at
            the end when chopping"""
        ),
    ),
    chopping_min_window: Optional[int] = typer.Option(
        None,
        help=dedent(
            """
            Smallest allowed length of beat (in
            milliseconds) to be included in chopping"""
        ),
    ),
    chopping_max_window: Optional[int] = typer.Option(
        None,
        help=dedent(
            """
            Largest allowed length of beat (in
            milliseconds) to be included in chopping"""
        ),
    ),
    ignore_pacing: bool = typer.Option(
        False,
        help=dedent(
            """
            Ignore pacing data, for example if the pacing is
            wrong"""
        ),
    ),
    overwrite: bool = typer.Option(
        True,
        help=dedent(
            """
            If True, overwrite existing data if outdir
            allready exist. If False, then the olddata will
            be copied to a subdirectory with version number
            of the software. If version number is not found
            it will be saved to a folder called "old"."""
        ),
    ),
    reuse_settings: bool = typer.Option(
        False,
        help=dedent(
            """
            If the output folder contains a file called
            settings.json and this flag is turned on, then
            we will use the settings stored in the
            settings.json file. This is handy if you e.g
            want to compare the output of the software
            between different versions, or you to reproduce
            the exact traces from the raw data."""
        ),
    ),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="More verbose"),
):

    scripts.analyze.main(
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
        verbose=verbose,
    )


@app.command(help=scripts.summary.__doc__)
def summary(
    folder: str = typer.Argument(..., help="The folder to be analyzed"),
    filename: str = typer.Option(
        "mps_summary",
        help=dedent(
            """
            Name of the pdf and csv file that is
            the output from the mps_summary script"""
        ),
    ),
    silent: bool = typer.Option(False, help="Turn of printing"),
    ignore_pacing: bool = typer.Option(
        False,
        help=dedent(
            """
            Ignore pacing data, for example if the pacing is
            wrong"""
        ),
    ),
    include_npy: bool = typer.Option(
        False,
        help=dedent(
            """
            If true then try to also open .npy
            files. The default behavious is not to
            include these, because the data that is
            analyzed is also dumped to a .npy
            files, and we do not want to open
            those."""
        ),
    ),
):
    scripts.summary.main(
        folder=folder,
        filename=filename,
        ignore_pacing=ignore_pacing,
        silent=silent,
        include_npy=include_npy,
    )


@app.command(help=scripts.mps2mp4.__doc__)
def mps2mp4(
    path: str = typer.Argument(..., help="Path to the mps file"),
    outfile: Optional[str] = typer.Option(
        None,
        "--outfile",
        "-o",
        help=dedent(
            """
            Output name for where you want to store the output
            movie. If not provided a the same name as the basename
            of the input file will be used"""
        ),
    ),
    synch: bool = typer.Option(
        False, help="Start video at same time as start of pacing"
    ),
):
    scripts.mps2mp4.main(path=path, outfile=outfile, synch=synch)


# def main():
#     """
#     Main execution of the mps package
#     """

#     if len(sys.argv) < 2:
#         print(__doc__)
#         return

#     # Show help message
#     if sys.argv[1] == "-h" or sys.argv[1] == "--help":
#         print(__doc__)

#     elif sys.argv[1] == "-v" or sys.argv[1] == "--version":
#         from mps import __version__

#         print(__version__)

#     elif sys.argv[1] == "-l" or sys.argv[1] == "--license":
#         print(__license__)

#     elif sys.argv[1] == "analyze":
#         bin_utils.analyze_mps.run(sys.argv[2:])

#     elif sys.argv[1] == "summary":
#         bin_utils.mps_summary.run(sys.argv[2:])

#     elif sys.argv[1] == "phase_plot":
#         bin_utils.mps_phase_plot.run(sys.argv[2:])

#     elif sys.argv[1] == "prevalence":
#         bin_utils.mps_prevalence.run(sys.argv[2:])

#     elif sys.argv[1] == "collect":
#         bin_utils.collect_mps.run(sys.argv[2:])

#     elif sys.argv[1] == "mps2mp4":
#         bin_utils.mps2mp4.run(sys.argv[2:])

#     elif sys.argv[1] == "split-pacing":
#         bin_utils.split_pacing(sys.argv[2:])

#     elif sys.argv[1] == "motion":
#         try:
#             import typer
#             from mps_motion_tracking import cli
#         except ImportError:
#             print("Motion tracking software not installed.")
#             print("Please ask Henrik (henriknf@simula.no)")
#             sys.exit()

#         # Run motion tracking
#         sys.argv[1:] = sys.argv[2:]
#         typer.run(cli.main)
#     elif sys.argv[1] == "automate":
#         try:
#             import typer
#             from mps_automation import cli
#         except ImportError:
#             print("Automation scripts are not installed.")
#             print("Please ask Henrik (henriknf@simula.no)")
#             sys.exit()
#         sys.argv[1:] = sys.argv[2:]
#         typer.run(cli.main)

#     else:
#         print("Argument {} not recongnized".format(sys.argv[1]))
#         print(__doc__)


if __name__ == "__main__":
    # main()
    app()
