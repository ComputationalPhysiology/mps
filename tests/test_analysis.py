import ap_features as apf
import numpy as np
import pytest

import mps


@pytest.fixture(params=["with_pacing", "without_pacing"])
def chopped_data(request):
    N = 700
    x = np.linspace(0, 7.0, N)
    time = x * 1000
    alpha = 0.74

    average = np.sin(2 * np.pi * (x + alpha))

    pacing = np.zeros(len(x))
    for r in range(8):
        pacing[1 + 100 * r : 100 * r + 10] = 1
    kwargs = dict(data=average, time=time, extend_front=0)
    if request.param == "with_pacing":
        kwargs["pacing"] = pacing
    else:
        kwargs["extend_front"] = 250

    kwargs["use_pacing_info"] = request.param == "with_pacing"

    yield apf.chopping.chop_data(**kwargs), kwargs


def test_local_averages():

    dt = 10.0
    T = 7000.0  # End time (ms)
    period = 1000
    c = 10

    nx = 150
    length = 500.0
    width = 200.0
    dx = length / nx
    ny = int(width / dx) + 1
    phi = 0.25

    M = int(T / dt + 1)

    x = np.linspace(0, length, nx)
    y = np.linspace(0, width, ny)
    X, Y = np.meshgrid(x, y)

    def Z(t):
        return 0.5 * np.sin(2 * np.pi * ((X / c - t) / period - phi)) + 0.5

    times = np.zeros(M)
    u = np.zeros((M, ny, nx))

    for i, t in enumerate(np.arange(0, T + dt, dt)):
        u[i, :, :] = Z(t)[:, :]
        times[i] = t

    frames = u.T

    t = np.arange(0, T + dt, dt)
    y = 0.5 * np.sin(-2 * np.pi * (phi + t / period)) + 0.5
    avg = mps.analysis.local_averages(frames, times, background_correction=False)[
        0, 0, :
    ]
    assert np.linalg.norm(avg - y) / np.linalg.norm(y) < 0.05


def test_analyze_apds(chopped_data):
    apd_analysis = mps.analysis.analyze_apds(
        chopped_data[0].data,
        chopped_data[0].times,
        plot=False,
    )

    # APD50 should be 500
    assert np.all(np.abs(np.subtract(apd_analysis.apds[50], 500)) < 1.0)


def test_analyze_frequencies(chopped_data):
    freq_analysis = mps.analysis.analyze_frequencies(
        chopped_data[0].data,
        chopped_data[0].times,
        plot=False,
    )
    assert np.all(np.abs(freq_analysis - 1.0) < 0.01)


def test_AnalyzeMPS(mps_data):
    avg = mps.average.get_average_all(mps_data.frames)
    analyzer = mps.analysis.AnalyzeMPS(
        avg,
        mps_data.time_stamps,
        mps_data.pacing,
        outdir="test_outdir",
        plot=True,
    )

    analyzer.analyze_all()
