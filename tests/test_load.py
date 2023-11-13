from pathlib import Path
from unittest import mock

import numpy as np

import mps

here = Path(__file__).absolute().parent
# czi_name = here.joinpath("data/demo.czi").as_posix()
# nd2_name = here.joinpath("data/voltage.nd2").as_posix()
# mps_data = mps.MPS(nd2_name)


# def test_save_load():
#     fname = "test_save.npy"
#     np.save(fname, mps_data)
#     loaded_data = np.load(fname, allow_pickle=True).item()

#     for k, v in mps_data.__dict__.items():
#         assert k in loaded_data.__dict__.keys()
#         if k == "data":
#             continue
#         if isinstance(v, np.ndarray):
#             assert np.all(v == loaded_data.__dict__.get(k))
#         else:
#             assert v == loaded_data.__dict__.get(k)


def create_dummy_data():
    dt = 10.0
    T = 7000.0  # End time (ms)
    period = 1000
    c = 10

    nx = 150
    length = 500.0
    width = 200.0
    dx = length / nx
    ny = int(width / dx) + 1
    um_per_pixel = dx
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

    pacing_duration = 50
    amp = 1
    pacing = np.zeros_like(times)
    period_rel = int(period / dt)
    for d in range(10, 10 + max(int(pacing_duration / dt), 1)):
        pacing[d::period_rel] = amp

    info = dict(num_frames=M, dt=dt, um_per_pixel=um_per_pixel, size_x=nx, size_y=ny)
    return dict(frames=u.T, time_stamps=times, info=info, pacing=pacing)


def test_from_dict():
    dummy_data = create_dummy_data()
    mps.MPS.from_dict(**dummy_data)


def test_nd2_file():
    dummy_data = create_dummy_data()

    class Frames:
        def __init__(self, frames, times) -> None:
            self._frames = frames
            self._times = times
            self.metadata = {
                "ImageMetadataSeqLV|0": {"SLxPictureMetadata": {"dCalibration": 0.42}},
            }

            self.height, self.width, self.imagecount = frames.shape

        def image(self, t):
            return self._frames[:, :, t]

        def get_time(self, t):
            return self._times[t]

    frames = Frames(dummy_data["frames"], dummy_data["time_stamps"])

    with mock.patch("mps.load.ND2File") as m:
        m.return_value.__enter__.return_value = frames
        data = mps.load.load_nd2("test_file.nd2")

    m.assert_called_with("test_file.nd2")
    assert np.isclose(dummy_data["time_stamps"], data.time_stamps).all()


if __name__ == "__main__":
    pass
    # test_save_load()
    # test_from_dict()
