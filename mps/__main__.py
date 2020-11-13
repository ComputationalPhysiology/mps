#!/usr/bin/env python3
__author__ = "Henrik Finsberg (henriknf@simula.no), 2017--2020"
__maintainer__ = "Henrik Finsberg"
__email__ = "henriknf@simula.no"
__license__ = """
c) 2001-2020 Simula Research Laboratory ALL RIGHTS RESERVED

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
__doc__ = """MPS package for analyzing microphyisological systems


Available arguments
-------------------
All these arguments can be called with '-h' or '--help' to see the
additional options

    analyze
        Analyze mps data (.nd2 or .czi)

    summary
        Create a summary figure and csv file of all files in a folder

    phase_plot
        Make a phase plot with voltage on the x-axis and calcium on the y-axis.

    prevalence
        Estimate the percentage of tissue in the chips,
        and the percentage of beating tissue vs non-beating

    collect
        Gather Voltage and Calcium data into one file that can be
        used as input to the inversion algorithm.

    mps2mp4
        Create movie of data file


Available options
-----------------

    -h, --help      Show this help
    -v, --version   Show version number
    -l, --license   Show license


Contact
-------
Henrik Finsberg (henriknf@simula.no)

"""

import sys

from mps import bin_utils


def main():
    """
    Main execution of the mps package
    """

    if len(sys.argv) < 2:
        print(__doc__)
        return

    # Show help message
    if sys.argv[1] == "-h" or sys.argv[1] == "--help":
        print(__doc__)

    elif sys.argv[1] == "-v" or sys.argv[1] == "--version":
        from mps import __version__

        print(__version__)

    elif sys.argv[1] == "-l" or sys.argv[1] == "--license":
        print(__license__)

    elif sys.argv[1] == "analyze":
        bin_utils.analyze_mps.run(sys.argv[2:])

    elif sys.argv[1] == "summary":
        bin_utils.mps_summary.run(sys.argv[2:])

    elif sys.argv[1] == "phase_plot":
        bin_utils.mps_phase_plot.run(sys.argv[2:])

    elif sys.argv[1] == "prevalence":
        bin_utils.mps_prevalence.run(sys.argv[2:])

    elif sys.argv[1] == "collect":
        bin_utils.collect_mps.run(sys.argv[2:])

    elif sys.argv[1] == "mps2mp4":
        bin_utils.mps2mp4.run(sys.argv[2:])

    else:
        print("Argument {} not recongnized".format(sys.argv[1]))
        print(__doc__)


if __name__ == "__main__":
    main()
