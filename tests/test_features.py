import os

from scipy.io import loadmat

from mps.analysis import apd, integrate_apd, maximum_upstroke_velocity, tau, upstroke

current_dir = os.path.abspath(os.path.dirname(__file__))
print("load ", os.path.join(current_dir, "data/data.mat"))
data = loadmat(os.path.join(current_dir, "data/data.mat"))
tol = 0.02
normalize = False
use_spline = False

v = data["v"].flatten()
t_v = data["t_v"].flatten()
ca = data["ca"].flatten()
t_ca = data["t_ca"].flatten()

outstr = "{:10} \tPython:{:10.4f}\tMATLAB:{:10.4f}"


def test_ADP30():

    APD30 = apd(v, 0.3, t_v, use_spline=use_spline)
    print(outstr.format("APD30", APD30, data["APD30"][0][0]))
    assert (data["APD30"][0][0] - APD30) / data["APD30"][0][0] < tol


def test_ADP50():

    APD50 = apd(v, 0.5, t_v, use_spline=use_spline)
    print(outstr.format("APD50", APD50, data["APD50"][0][0]))
    assert (data["APD50"][0][0] - APD50) / data["APD50"][0][0] < tol


def test_ADP80():

    APD80 = apd(v, 0.8, t_v, use_spline=use_spline)
    print(outstr.format("APD80", APD80, data["APD80"][0][0]))
    assert (data["APD80"][0][0] - APD80) / data["APD80"][0][0] < tol


def test_CaD30():

    CaD30 = apd(ca, 0.3, t_ca, use_spline=use_spline)
    print(outstr.format("CaD30", CaD30, data["CaD30"][0][0]))
    assert (data["CaD30"][0][0] - CaD30) / data["CaD30"][0][0] < tol


def test_CaD50():

    CaD50 = apd(ca, 0.5, t_ca, use_spline=use_spline)
    print(outstr.format("CaD50", CaD50, data["CaD50"][0][0]))
    assert (data["CaD50"][0][0] - CaD50) / data["CaD50"][0][0] < tol


def test_CaD80():

    CaD80 = apd(ca, 0.8, t_ca, use_spline=use_spline)
    print(outstr.format("CaD80", CaD80, data["CaD80"][0][0]))
    assert (data["CaD80"][0][0] - CaD80) / data["CaD80"][0][0] < tol


def test_dvdt_max():

    dvdt_max = maximum_upstroke_velocity(
        v, t_v, use_spline=use_spline, normalize=normalize
    )
    print(outstr.format("dvdt_max", dvdt_max, data["dvdt_max"][0][0]))
    assert (data["dvdt_max"][0][0] - dvdt_max) / data["dvdt_max"][0][0] < tol


def test_dcdt_max():

    dcdt_max = maximum_upstroke_velocity(
        ca, t_ca, use_spline=use_spline, normalize=normalize
    )
    print(outstr.format("dcdt_max", dcdt_max, data["dcdt_max"][0][0]))
    assert (data["dcdt_max"][0][0] - dcdt_max) / data["dcdt_max"][0][0] < tol


def test_int_30():

    int_30 = integrate_apd(v, t_v, 0.3, use_spline=use_spline, normalize=normalize)
    print(outstr.format("int_30", int_30, data["int_30"][0][0]))
    assert (data["int_30"][0][0] - int_30) / data["int_30"][0][0] < tol


def test_tau():

    t = tau(t_v, v, 0.75)
    assert (t - 261.494424349009) < 1e-6
    print("tau75 = ", t)


def test_upstroke():

    t = upstroke(t_v, v, 0.8)
    assert (t - 155.54625401610127) < 1e-6
    print("Upstroke80 = ", t)


if __name__ == "__main__":
    # test_ADP30()
    # test_ADP50()
    # test_ADP80()
    # test_CaD30()
    # test_CaD50()
    # test_CaD80()
    # test_dvdt_max()
    # test_dcdt_max()
    # test_int_30()
    # test_tau()
    test_upstroke()
