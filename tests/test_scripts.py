import shutil
import subprocess as sp
from pathlib import Path

from distutils.spawn import find_executable

import mps

try:
    import matplotlib  # noqa: F401
except ImportError:
    missing_mpl = True
else:
    missing_mpl = False

try:
    import imageio_ffmpeg  # noqa: F401
except ImportError:
    missing_ffmpeg = True
else:
    missing_ffmpeg = False


import pytest

python = find_executable("python")


def test_mps_analyze(mps_data_path):
    ret = sp.call([python, "-m", "mps", "analyze", mps_data_path])
    assert ret == 0
    shutil.rmtree(Path(mps_data_path).with_suffix(""))


@pytest.mark.skipif(missing_mpl, reason="Requires matplotlib")
def test_mps_phase_plot(mps_data_path):
    out = Path(mps_data_path).parent.joinpath("phase_plot.png")
    ret = sp.call(
        [
            python,
            "-m",
            "mps",
            "phase-plot",
            mps_data_path,
            mps_data_path,
            "-o",
            out.as_posix(),
        ],
    )
    assert ret == 0
    out.unlink()


@pytest.mark.skipif(missing_ffmpeg, reason="Requires imageio-ffmpeg")
def test_mps2mp4(mps_data_path):
    out = Path(mps_data_path).parent.joinpath("movie.mp4")
    ret = sp.call([python, "-m", "mps", "mps2mp4", mps_data_path, "-o", out.as_posix()])
    assert ret == 0
    out.unlink()


@pytest.mark.skipif(missing_mpl, reason="Requires matplotlib")
def test_mps_summary_script(mps_data_path):
    path = Path(mps_data_path)
    new_path = path.parent.joinpath("data").joinpath(path.name)
    another_path = path.parent.joinpath("data").joinpath(f"another_{path.name}")
    new_path.parent.mkdir(exist_ok=True)
    shutil.copy(path, new_path)
    shutil.copy(path, another_path)

    ret = sp.call(
        [
            python,
            "-m",
            "mps",
            "summary",
            "--include-npy",
            new_path.parent.absolute().as_posix(),
        ],
    )

    assert ret == 0

    shutil.rmtree(new_path.parent)


@pytest.mark.skipif(missing_mpl, reason="Requires matplotlib")
def test_mps_summary_function(mps_data_path):
    path = Path(mps_data_path)
    new_path = path.parent.joinpath("data").joinpath(path.name)
    another_path = path.parent.joinpath("data").joinpath(f"another_{path.name}")
    new_path.parent.mkdir(exist_ok=True)
    shutil.copy(path, new_path)
    shutil.copy(path, another_path)
    mps.scripts.summary.main(new_path.parent, include_npy=True)


@pytest.mark.xfail(reason="Not yet implemented properly")
def test_mps_prevalence(mps_data_path):
    ret = sp.call([python, "-m", "mps", "prevalence", mps_data_path])
    assert ret == 0
