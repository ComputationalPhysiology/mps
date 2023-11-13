from pathlib import Path

import numpy as np
import pytest

import mps


def ca_transient(t, tstart=0.05):
    tau1 = 0.05
    tau2 = 0.110

    ca_diast = 0.0
    ca_ampl = 1.0

    beta = (tau1 / tau2) ** (-1 / (tau1 / tau2 - 1)) - (tau1 / tau2) ** (-1 / (1 - tau2 / tau1))
    ca = np.zeros_like(t)

    ca[t <= tstart] = ca_diast

    ca[t > tstart] = (ca_ampl - ca_diast) / beta * (
        np.exp(-(t[t > tstart] - tstart) / tau1) - np.exp(-(t[t > tstart] - tstart) / tau2)
    ) + ca_diast
    return ca


@pytest.fixture
def mps_data_path():
    num_beats = 6
    size_beat = 100
    num_pacing_steps = 5
    tstart = 0.01 * num_pacing_steps
    num_frames = num_beats * size_beat
    size_x = 20
    size_y = 20
    frames = np.zeros((size_x, size_y, num_frames))
    time_stamps = np.zeros(num_frames)
    pacing = np.zeros(num_frames)
    pacing[size_beat]

    for beat in range(num_beats):
        t = np.linspace(0, 1, size_beat + 1)
        y = ca_transient(t[:-1], tstart=tstart)

        frames[:, :, size_beat * beat : size_beat * (beat + 1)] = np.tile(
            y,
            (size_x, size_y, 1),
        )
        time_stamps[size_beat * beat : size_beat * (beat + 1)] = t[:-1] * 1000 + beat * 1000
        pacing[beat * size_beat : beat * size_beat + num_pacing_steps] = 1

    info = dict(
        num_frames=num_frames,
        dt=np.mean(np.diff(time_stamps)),
        time_unit="ms",
        size_x=size_x,
        size_y=size_y,
    )

    data = dict(
        frames=frames,
        time_stamps=time_stamps,
        info=info,
        pacing=np.zeros_like(time_stamps),
    )
    path = Path(__file__).parent.joinpath("test_data.npy").absolute()
    np.save(path, data)
    yield path.as_posix()

    path.unlink()


@pytest.fixture
def mps_data(mps_data_path):
    return mps.MPS(mps_data_path)


if __name__ == "__main__":
    pass
    # t = np.linspace(0, 1, 100)
    # import matplotlib.pyplot as plt

    # y = ca_transient(t)
    # plt.plot(t * 1000, y)
    # plt.show()
    # mps_data()
    # mps_data_path()
