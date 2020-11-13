#!/usr/bin/env python3
__author__ = "Henrik Finsberg (henriknf@simula.no), 2017--2019"
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

import functools as ft
import itertools as it

import numpy as np
from scipy.interpolate import UnivariateSpline

from . import utils

logger = utils.get_logger(__name__)


def get_subsignal_average(signals):
    """
    Get average of signals.
    Assume that signals are alinged, but
    they dont have to be of the same length

    Arguments
    ---------
    signals : array
        The data you want to average

    """

    if len(signals) == 0:
        return []

    if len(signals) == 1:
        return signals[0]

    # Check is they have the same lenght
    if all([len(s) == len(signals[0]) for s in signals[1:]]):
        # Then it is easy to take the average
        average = np.mean(signals, 0)

    else:
        # We need to take into account the possibilty
        # the the subsignals have different lenght
        try:
            # Python 2
            from itertools import izip_longest

        except ImportError:
            # Python 3
            from itertools import zip_longest as izip_longest

        def avg(x):
            x = [i for i in x if i]
            if len(x) == 0:
                return 0.0
            return sum(x, 0.0) / len(x)

        average = np.array(tuple(map(avg, izip_longest(*signals))))

    return average


def get_subsignal_average_interpolate(ys, xs=None, N=200, return_list=False):
    """
    Get the avagere of list of signals assuming that
    they align at the same x value

    Arguments
    ---------
    ys : list of arrays
        The signal values
    xs : list of arrays
        The x-values
    N : int
        Lenght of output array (Default: 200)

    Returns
    -------
    Y_avg: array
        The average y values
    X : array
        The new x-values

    """
    if xs is None or len(xs) == 0:
        return get_subsignal_average(ys), range(max([len(i) for i in ys]))

    # Construct new time array
    min_x = np.min([xi[0] for xi in xs])
    max_x = np.max([xi[-1] for xi in xs])
    X = np.linspace(min_x, max_x, N)

    if len(ys) == 0:
        return np.zeros(N), X

    # Check args
    msg = (
        "Expected Xs and Ys has to be of same lenght. " "Got len(xs) = {}, len(ys) = {}"
    ).format(len(xs), len(ys))
    assert len(xs) == len(ys), msg

    for i, (x, y) in enumerate(zip(xs, ys)):
        msg = (
            "Expected X and Y has to be of same lenght. "
            "Got len(x) = {}, len(y) = {} for index {}"
        ).format(len(x), len(y), i)
        assert len(x) == len(y), msg

    Ys = []
    Xs = []

    for i, (x, y) in enumerate(zip(xs, ys)):
        # Take out relevant piece
        idx = next(j + 1 for j, xj in enumerate(X) if xj >= x[-1])

        X_ = X[:idx]
        Xs.append(X_)

        Y = UnivariateSpline(x, y, s=0)(X_)
        Ys.append(Y)

    Y_avg = get_subsignal_average(Ys)

    if return_list:
        return Y_avg, X, Ys, Xs

    return Y_avg, X


def get_average_all(X):
    r"""
    Return the global average of a frame stack, i.e

    .. math::
        y = \begin{pmatrix} y^1, y^2, \cdots , y^T \end{pmatrix}

    with

    .. math::
        y^i = \mathrm{avg}_{\mathrm{all}}(\mathbf{X}^i )
        = \bar{X}^i = \frac{1}{MN} \sum_{m,n} x_{m,n}^i

    Arguments
    ---------
    X : :class:`numpy.ndarray`
        The frame stack (size :math:`M \times N \times T`)

    Returns
    -------
    average : :class:`numpy.ndarray`
        The average vector of length :math:`T`


    """

    N = X.shape[0] * X.shape[1]
    x = ((1.0 / N) * x2 for x1 in X for x2 in x1)
    logger.debug("Compute global average")
    return ft.reduce(np.add, x, np.zeros(X.shape[-1]))


def get_spatial_average(X, alpha=0.5, **kwargs):
    r"""
    Get the spatial average

    Arguments
    ---------
    X : :class:`numpy.ndarray`
        The frame stack (size :math:`M \times N \times T`)
    alpha : float

    Returns
    -------
    average : :class:`numpy.ndarray`
        The average vector of length :math:`T`


    Notes
    -----

    If we assume that the important pixels are those pixels with a high
    pixel value, then we can simply select those pixel with the highest
    pixels values. Formally, for each frame :math:`i`,

    .. math::
        \mathcal{D} = \{ (m,n) : x_{m,n} > f(\mathbf{X}^i) \}

    then,

    .. math::
        y^i = \mathrm{avg}_{\mathrm{spatial}, f}(\mathbf{X}^i )
        = \frac{1}{|\mathcal{D}|} \sum_{(m,n) \in \mathcal{D}} x_{m,n}^i

    where :math:`|\mathcal{D}|` is the number of pixels in
    :math:`\mathcal{D}`.

    """
    N = int((X.shape[0] * X.shape[1]) * (1 - alpha))
    x = (np.mean(np.sort(x.flatten())[N:]).tolist() for x in X.T)
    logger.debug("Compute spatial average")
    return np.array(tuple(it.chain(x)))


def get_temporal_average(X, alpha=0.5, **kwargs):
    r"""
    Get the temporal average

    Arguments
    ---------
    X : :class:`numpy.ndarray`
        The frame stack (size :math:`M \times N \times T`)
    bound : str
        Make average of pixels larger than the bound.
        Choices ('median', 'mean'). Default . 'median'

    Returns
    -------
    average : :class:`numpy.ndarray`
        The average vector of length :math:`T`


    Notes
    -----

    If we assume that pixels that represents noise have pixel values,
    that do not vary much throughout the frame stack compared to signals
    that represents the actualt signals, we can try to include only
    those pixels that have a certain variation compared to all the
    other pixels.
    Let

    .. math::
        \sigma(\mathcal{X}) = \begin{pmatrix}
        \sigma(X)_{1,1} & \sigma(X)_{1,2} & \cdots & \sigma(X)_{1,N} \\
        \sigma(X)_{2,1} & \sigma(X)_{2,2} & \cdots & \sigma(X)_{2,N} \\
        \vdots  & \vdots & \vdots & \vdots \\
        \sigma(X)_{M,1} & \sigma(X)_{M,2} & \cdots & \sigma(X)_{M,N}
        \end{pmatrix}

    be the :math:`M \times N` matrix with each coordinate beeing
    the temporal standard deviation.
    Now, for a given fuctional :math:`f` let

    .. math::
        \mathcal{D} = \{ (m,n) : \sigma(X)_{m,n}
        > f(\sigma(\mathcal{X})) \}

    be the set of pixel coordinates selecting pixels based on the rule
    :math:`f`.
    Then

    .. math::
        y^i = \mathrm{avg}_{\mathrm{temporal}, f}(\mathbf{X}^i )
        = \frac{1}{|\mathcal{D}|} \sum_{(m,n) \in \mathcal{D}} x_{m,n}^i

    where :math:`|\mathcal{D}|` is the number of pixels in $\mathcal{D}$.


    """
    logger.debug("Compute temporal average")
    # Reshape to get a matrix of size num_pixels times num_time_steps
    q = X.T.reshape((X.shape[-1], -1))
    # Find the temporal variation
    temp_std = np.std(X, 2)

    # Get the relevant indices
    X_flat = temp_std.reshape(-1)
    X_sort_arg = np.argsort(X_flat)
    idx = int(np.floor(alpha * len(X_flat)))
    inds = X_sort_arg[idx:]

    # Compute the average
    average_intensity = np.mean(q[:, inds], 1)

    return_indices = kwargs.get("return_indices", False)
    if return_indices:
        return average_intensity, inds

    return average_intensity


def masked_average(X, mask):
    return np.mean(X[mask], axis=0)
