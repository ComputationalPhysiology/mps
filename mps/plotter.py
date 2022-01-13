#!/usr/bin/env python3
__author__ = "Henrik Finsberg (henriknf@simula.no), 2017--2012"
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
import logging

logger = logging.getLogger(__name__)

try:
    import matplotlib as mpl

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        mpl.use("TkAgg")
        import matplotlib.pyplot as plt

    has_mpl = True
except ImportError:
    logger.debug(("Unable to import matplotlib. " "Plotting will not be possible"))
    has_mpl = False
    plt = None

try:
    import seaborn as sns
except ImportError:
    pass
else:
    sns.set_context()


def set_mpl_options(**kwargs):
    if not has_mpl:
        return

    mpl.rcParams["savefig.format"] = kwargs.pop("save_format", "pdf")
    from distutils.spawn import find_executable

    if find_executable("latex"):
        mpl.rcParams["text.usetex"] = kwargs.pop("usetex", True)


set_mpl_options()


def plot_window(self, frames, num_frames, local=None, t=0, fname=None):
    """
    If you can to take average over a smaller region
    you can use this function to plot a frame at the
    given frame

    Arguments
    ---------
    local : array
        Return list of indice to plot as frame
        [startx, endx, starty, endy]
    t : int
        Image number
    fname : str:
        Name of figure to save

    """
    if not has_mpl:
        return

    assert isinstance(t, int)
    assert 0 <= t < num_frames

    i = frames[:, :, t]
    fig, ax = plt.subplots()
    ax.imshow(i, cmap="hot")
    if local:
        ax.plot([local[0], local[0]], [local[2], local[3]], color="b")
        ax.plot([local[1], local[1]], [local[2], local[3]], color="b")
        ax.plot([local[0], local[1]], [local[2], local[2]], color="b")
        ax.plot([local[0], local[1]], [local[3], local[3]], color="b")

    if fname:
        fig.savefig(fname)
    else:
        plt.show()
    plt.close()


def plot_twin_trace(x1, x2, y1, y2, fname):

    if not has_mpl:
        return

    fig, ax1 = plt.subplots()
    ax1.plot(x1, y1, color="r")
    ax1.set_xlabel("Time (ms)")
    ax2 = ax1.twinx()
    ax2.plot(x2, y2, color="b")
    fig.savefig(fname)
    plt.close()


def plot_single_trace(x, y, fname, ylabel="Pixel intensity"):

    if not has_mpl:
        return

    fig, ax = plt.subplots()
    ax.plot(x, y)
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel(ylabel)
    fig.savefig(fname)
    plt.close()


def plot_multiple_traces(xs, ys, fname, titles=None, deep=False, ylabels=None):

    if not has_mpl:
        return

    if titles is None:
        titles = ["" for x in xs]
    if ylabels is None:
        ylabels = ["" for x in xs]

    N = len(xs)
    if N in [2, 4]:
        size = (N // 2, 2)
    else:
        while N % 3 != 0:
            N += 1
        size = (N // 3, 3)
    fig, axs = plt.subplots(*size)

    for i, (x, y, t, ylabel) in enumerate(zip(xs, ys, titles, ylabels)):
        ax = axs.flatten()[i]
        if deep:
            for xi, yi in zip(x, y):
                ax.plot(xi, yi)
        else:
            ax.plot(x, y)

        ax.set_xlabel("Time (ms)")
        ax.set_ylabel(ylabel)
        ax.set_title(" ".join(t.split("_")))

    fig.tight_layout()
    fig.savefig(fname)
    plt.close()


def phase_plots(voltage, calcium, fname):

    fig, ax = plt.subplots()
    ax.plot(voltage, calcium, "o")
    ax.set_xlabel("Voltage")
    ax.set_ylabel("Calcium")
    fig.savefig(fname)
    plt.close()
