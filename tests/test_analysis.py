import numpy as np
import pytest
import scipy.spatial.distance

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

    yield mps.analysis.chop_data(**kwargs), kwargs


def test_filt():

    N = 700
    x = np.linspace(0, 7.0, N)
    average = np.sin(2 * np.pi * (x - 0.25))
    smooth = mps.analysis.filt(average)

    assert np.linalg.norm(average - smooth) / np.linalg.norm(average)


def test_beating_frequency():

    assert (
        abs(mps.analysis.beating_frequency([range(101) for i in range(10)], "s") - 0.01)
        < 1e-12
    )
    assert (
        abs(
            mps.analysis.beating_frequency([range(101) for i in range(10)], "ms") - 10.0
        )
        < 1e-12
    )

    assert (
        abs(
            mps.analysis.beating_frequency(
                [range(10 * i, 101 + 10 * i) for i in range(10)]
            )
            - 10.0
        )
        < 1e-12
    )


def test_background():
    N = 700
    x = np.linspace(0, 7, N)
    y = np.sin(2 * np.pi * 1.2 * x)

    a = 0.1
    b = -1.5
    c = 10.0

    # Test linear and quadratic background
    for a in [0, 0.1]:

        background = a * x ** 2 + b * x + c
        background -= background[0]
        signal = y + background
        estimated_background = mps.analysis.background(x, signal)
        corrected = mps.analysis.correct_background(x, signal)

        z1 = background - background[-1]
        z2 = estimated_background - estimated_background[-1]

        print(abs(np.subtract(z1, z2)))
        assert all(abs(np.subtract(z1, z2)) < 1e-2)
        assert all(abs(np.subtract(corrected, y)))


def test_chop_data(chopped_data):

    chopped_data, kwargs = chopped_data
    # assert len(chopped_data.data) == 7
    print([len(c) for c in chopped_data.data])
    assert chopped_data.parameters.use_pacing_info is kwargs["use_pacing_info"]
    assert chopped_data.parameters.extend_front == kwargs["extend_front"]

    # fig, ax = plt.subplots()
    # for t, c in zip(chopped_data.times, chopped_data.data):
    #     ax.plot(t, c)
    # plt.show()
    N = min([len(d) for d in chopped_data.data])
    data = np.array([d[:N] for d in chopped_data.data])
    q = scipy.spatial.distance.pdist(data, "euclidean") / max(
        [np.linalg.norm(d) for d in data]
    )
    assert all(q < 0.1)

    times = np.array([t[:N] for t in chopped_data.times])
    assert all(scipy.spatial.distance.pdist([t - t[0] for t in times]) < 1e-10)


@pytest.mark.parametrize(
    "start_in, end_in, start_out, end_out",
    (
        ([0, 2, 4], [1, 3, 5], [0, 2, 4], [1, 3, 5]),
        ([0, 2, 4, 6], [1, 3, 5], [0, 2, 4], [1, 3, 5]),
        ([-1, 2, 4, 6], [1, 3, 5], [2, 4], [3, 5]),
        ([4], [1, 3, 5], [4], [5]),
        ([1, 2], [1, 2, 5], [1, 2], [2, 5]),
    ),
)
def test_filter_start_ends_in_chopping(start_in, end_in, start_out, end_out):

    s, e = mps.analysis.filter_start_ends_in_chopping(start_in, end_in)
    print(e)
    print(s)
    assert np.all(s == np.array(start_out))
    assert np.all(e == np.array(end_out))


@pytest.mark.parametrize(
    "start_in, end_in", (([], [1, 3, 5]), ([1, 3], []), ([], []), ([2, 3], [1]))
)
def test_filter_start_ends_in_chopping_raises_EmptyChoppingError(start_in, end_in):

    with pytest.raises(mps.analysis.EmptyChoppingError):
        mps.analysis.filter_start_ends_in_chopping(start_in, end_in)


def test_remove_spikes():
    spike_pt = 2
    N = 10
    pacing = np.zeros(N)
    pacing[spike_pt + 1] = 1
    arr = np.arange(N)

    for spike_dur in [0, 1, 4, 8]:
        arr1 = mps.analysis.remove_spikes(arr, pacing, spike_dur)
        assert all(
            arr1
            == np.concatenate((np.arange(spike_pt), np.arange(spike_pt + spike_dur, N)))
        )


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


def _test_prevalence():
    # from test_load import nd2_name
    path = "/Users/henriknf/Dropbox/Aslak/mps_conduction_velocity/report/notebooks/example_file.nd2"
    mps_data = mps.MPS(path)
    mps.analysis.prevalence(mps_data)


def test_angle2neighhbor():
    eps = 0.001
    tol = 1e-8
    angles = np.arange(-np.pi, np.pi, np.pi / 4)

    for a, a_p, a_m in zip(angles, angles + eps, angles - eps):
        n = mps.med64_utils.angle2neighbor(a)
        n_p = mps.med64_utils.angle2neighbor(a_p)
        n_m = mps.med64_utils.angle2neighbor(a_m)

        assert n == n_p == n_m
        ai = np.arctan2(n[1], n[0])
        msg = f"a = {a} != ai = {ai}"

        if abs(abs(a) - np.pi) < tol:
            # Could be both -pi and pi
            assert abs(a - np.pi) < tol or abs(a + np.pi) < tol, msg
        else:
            assert a == ai, msg
        print(f"{a:.2f}, {n}, {ai:.2f}")


def test_analyze_apds(chopped_data):
    apd_analysis = mps.analysis.analyze_apds(
        chopped_data[0].data, chopped_data[0].times, plot=False
    )

    # APD50 should be 500
    assert np.all(np.abs(np.subtract(apd_analysis.apds[50], 500)) < 1.0)


def test_analyze_frequencies(chopped_data):
    freq_analysis = mps.analysis.analyze_frequencies(
        chopped_data[0].data, chopped_data[0].times, plot=False
    )
    assert np.all(np.abs(freq_analysis - 1.0) < 0.01)


def test_bumpy_signals():
    data = np.load("bumpy_data.npy", allow_pickle=True).item()

    time = data["x"]
    trace = data["y"]
    pacing = data["p"]

    q = mps.analysis.analyze_mps_func_avg(trace, time, pacing)

    assert q["features"]["apd90"] > q["features"]["apd80"]


# def test_poincare_plot(chopped_data):
#     mps.analysis.poincare_plot(chopped_data[0].data, chopped_data[0].times, apds=[0.8])


if __name__ == "__main__":
    # pass
    test_bumpy_signals()
    # test_angle2neighhbor()
    # test_condction_velocity()
    # test_griddata()
    # test_local_averages()

    # test_filter_start_ends_in_chopping()
    # r = lambda: None
    # r.param = "without_pacing"
    # c = chopped_data(r)
    # test_analyze_apds(next(c))
