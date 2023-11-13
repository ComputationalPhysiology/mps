import numpy as np
import pytest

import mps


@pytest.mark.parametrize("avg_type", ["all", "spatial", "temporal"])
def test_average(mps_data, avg_type):
    expected = mps_data.frames[0, 0, :]

    if avg_type == "all":
        computed = mps.average.get_average_all(mps_data.frames)
    elif avg_type == "spatial":
        computed = mps.average.get_spatial_average(mps_data.frames, alpha=0.5)
    elif avg_type == "temporal":
        computed = mps.average.get_temporal_average(mps_data.frames, alpha=0.5)

    assert np.linalg.norm(np.subtract(computed, expected)) < 1e-12


def test_average_ones():
    tol = 1e-12
    data = np.ones((20, 20, 20))

    avg0 = mps.average.get_average_all(data)
    avg1 = mps.average.get_spatial_average(data, alpha=0.5)
    avg2 = mps.average.get_temporal_average(data, alpha=0.5)

    for avg in [avg0, avg1, avg2]:
        assert all(avg - 1 < tol)


if __name__ == "__main__":
    test_average_ones()
    # data = czi_data()
    # test_average_czi(data)
