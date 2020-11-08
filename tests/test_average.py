import os

import numpy as np
import pytest

import mps

current_dir = os.path.abspath(os.path.dirname(__file__))
czi_name = os.path.join(current_dir, "data/demo.czi")
nd2_name = os.path.join(current_dir, "data/voltage.nd2")
average_data = np.load(
    os.path.join(current_dir, "data/average.npy"), allow_pickle=True
).item()


@pytest.fixture
def czi_data():

    data = mps.MPS(czi_name)
    return data


@pytest.fixture
def nd2_data():

    data = mps.MPS(nd2_name)
    return data


@pytest.fixture(params=["all", "spatial", "temporal"])
def average_nd2(request, nd2_data):

    if request.param == "all":
        return (
            mps.average.get_average_all(nd2_data.frames),
            average_data["nd2_average_all"],
        )
    elif request.param == "spatial":
        return (
            mps.average.get_spatial_average(nd2_data.frames),
            average_data["nd2_spatial_average"],
        )
    elif request.param == "temporal":
        return (
            mps.average.get_temporal_average(nd2_data.frames),
            average_data["nd2_temporal_average"],
        )


def test_average_nd2(average_nd2):

    computed, expected = average_nd2
    assert np.linalg.norm(np.subtract(computed, expected)) < 1e-12


@pytest.fixture(params=["all", "spatial", "temporal"])
def average_czi(request, czi_data):

    if request.param == "all":
        return (
            mps.average.get_average_all(czi_data.frames),
            average_data["czi_average_all"],
        )
    elif request.param == "spatial":
        return (
            mps.average.get_spatial_average(czi_data.frames),
            average_data["czi_spatial_average"],
        )
    elif request.param == "temporal":
        return (
            mps.average.get_temporal_average(czi_data.frames),
            average_data["czi_temporal_average"],
        )


def test_average_czi(average_czi):

    computed, expected = average_czi
    assert np.linalg.norm(np.subtract(computed, expected)) < 1e-12


def test_average_ones():

    tol = 1e-12
    data = np.ones((20, 20, 20))

    avg0 = mps.average.get_average_all(data)
    avg1 = mps.average.get_spatial_average(data)
    avg2 = mps.average.get_temporal_average(data)

    for avg in [avg0, avg1, avg2]:
        assert all(avg - 1 < tol)


if __name__ == "__main__":
    test_average_ones()
    # data = czi_data()
    # test_average_czi(data)
