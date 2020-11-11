from pathlib import Path

import pytest

from mps.analysis import apd, integrate_apd, maximum_upstroke_velocity, tau, upstroke
from mps.utils import loadmat

tol = 0.02
normalize = False
use_spline = False


@pytest.fixture
def matlab_data():
    path = Path(__file__).parent.joinpath("data.mat")
    return loadmat(path)


outstr = "{:10} \tPython:{:10.4f}\tMATLAB:{:10.4f}"


def test_ADP30(matlab_data):

    v = matlab_data["v"]
    t_v = matlab_data["t_v"]

    APD30 = apd(v, 0.3, t_v, use_spline=use_spline)
    print(outstr.format("APD30", APD30, matlab_data["APD30"]))
    assert (matlab_data["APD30"] - APD30) / matlab_data["APD30"] < tol


def test_ADP50(matlab_data):

    v = matlab_data["v"]
    t_v = matlab_data["t_v"]

    APD50 = apd(v, 0.5, t_v, use_spline=use_spline)
    print(outstr.format("APD50", APD50, matlab_data["APD50"]))
    assert (matlab_data["APD50"] - APD50) / matlab_data["APD50"] < tol


def test_ADP80(matlab_data):

    v = matlab_data["v"]
    t_v = matlab_data["t_v"]

    APD80 = apd(v, 0.8, t_v, use_spline=use_spline)
    print(outstr.format("APD80", APD80, matlab_data["APD80"]))
    assert (matlab_data["APD80"] - APD80) / matlab_data["APD80"] < tol


def test_CaD30(matlab_data):

    ca = matlab_data["ca"]
    t_ca = matlab_data["t_ca"]

    CaD30 = apd(ca, 0.3, t_ca, use_spline=use_spline)
    print(outstr.format("CaD30", CaD30, matlab_data["CaD30"]))
    assert (matlab_data["CaD30"] - CaD30) / matlab_data["CaD30"] < tol


def test_CaD50(matlab_data):

    ca = matlab_data["ca"]
    t_ca = matlab_data["t_ca"]

    CaD50 = apd(ca, 0.5, t_ca, use_spline=use_spline)
    print(outstr.format("CaD50", CaD50, matlab_data["CaD50"]))
    assert (matlab_data["CaD50"] - CaD50) / matlab_data["CaD50"] < tol


def test_CaD80(matlab_data):

    ca = matlab_data["ca"]
    t_ca = matlab_data["t_ca"]

    CaD80 = apd(ca, 0.8, t_ca, use_spline=use_spline)
    print(outstr.format("CaD80", CaD80, matlab_data["CaD80"]))
    assert (matlab_data["CaD80"] - CaD80) / matlab_data["CaD80"] < tol


def test_dvdt_max(matlab_data):

    v = matlab_data["v"]
    t_v = matlab_data["t_v"]

    dvdt_max = maximum_upstroke_velocity(
        v, t_v, use_spline=use_spline, normalize=normalize
    )
    print(outstr.format("dvdt_max", dvdt_max, matlab_data["dvdt_max"]))
    assert (matlab_data["dvdt_max"] - dvdt_max) / matlab_data["dvdt_max"] < tol


def test_dcdt_max(matlab_data):

    ca = matlab_data["ca"]
    t_ca = matlab_data["t_ca"]

    dcdt_max = maximum_upstroke_velocity(
        ca, t_ca, use_spline=use_spline, normalize=normalize
    )
    print(outstr.format("dcdt_max", dcdt_max, matlab_data["dcdt_max"]))
    assert (matlab_data["dcdt_max"] - dcdt_max) / matlab_data["dcdt_max"] < tol


def test_int_30(matlab_data):

    v = matlab_data["v"]
    t_v = matlab_data["t_v"]

    int_30 = integrate_apd(v, t_v, 0.3, use_spline=use_spline, normalize=normalize)
    print(outstr.format("int_30", int_30, matlab_data["int_30"]))
    assert (matlab_data["int_30"] - int_30) / matlab_data["int_30"] < tol


def test_tau(matlab_data):

    v = matlab_data["v"]
    t_v = matlab_data["t_v"]

    t = tau(t_v, v, 0.75)
    assert (t - 261.494424349009) < 1e-6
    print("tau75 = ", t)


def test_upstroke(matlab_data):

    v = matlab_data["v"]
    t_v = matlab_data["t_v"]

    t = upstroke(t_v, v, 0.8)
    assert (t - 155.54625401610127) < 1e-6
    print("Upstroke80 = ", t)
