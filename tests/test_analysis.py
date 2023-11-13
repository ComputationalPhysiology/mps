import ap_features as apf
import numpy as np
import pytest

import mps


@pytest.fixture(params=["with_pacing", "without_pacing"])
def beats(request):
    N = 700
    x = np.linspace(0, 7.0, N)
    time = x * 1000
    alpha = 0.74

    average = np.sin(2 * np.pi * (x + alpha))

    pacing = np.zeros(len(x))
    for r in range(8):
        pacing[1 + 100 * r : 100 * r + 10] = 1

    chopping_options = {"extend_front": 0}
    if request.param == "without_pacing":
        chopping_options["extend_front"] = 250

    chopping_options["ignore_pacing"] = request.param == "without_pacing"
    beat = apf.Beats(
        t=time,
        y=average,
        pacing=pacing,
        chopping_options=chopping_options,
    )

    yield beat.beats


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
        0,
        0,
        :,
    ]
    assert np.linalg.norm(avg - y) / np.linalg.norm(y) < 0.05


def test_analyze_apds(beats):
    apd_analysis = mps.analysis.analyze_apds(
        beats,
        plot=False,
    )

    # APD50 should be 500
    assert np.all(np.abs(np.subtract(apd_analysis.apds[50], 500)) < 1.0)


def test_analyze_frequencies(beats):
    freq_analysis = mps.analysis.analyze_frequencies(
        beats,
        plot=False,
    )
    assert np.all(np.abs(freq_analysis - 1.0) < 0.01)


def test_AnalyzeMPS(mps_data):
    mps.analysis.analyze_mps_func(mps_data, plot=False)


def test_compare_get_average_all_and_local(mps_data):
    analysis_data = mps.analysis.analyze_mps_func(mps_data)
    avg1 = analysis_data["unchopped_data"]["trace"]
    avg2 = mps.analysis.frame2average(
        mps_data.frames,
        mps_data.time_stamps,
        background_correction=True,
    )
    assert np.isclose(avg1, avg2).all()
