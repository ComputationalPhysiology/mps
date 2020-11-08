import os
import subprocess as sp
from distutils.spawn import find_executable

python = find_executable("python")
current_dir = os.path.abspath(os.path.dirname(__file__))
voltage_name = os.path.join(current_dir, "data/voltage.nd2")
calcium_name = os.path.join(current_dir, "data/calcium.nd2")
voltage_dir = os.path.join(current_dir, "data/voltage")
calcium_dir = os.path.join(current_dir, "data/calcium")
brightfield_name = os.path.join(current_dir, "data/brightfield.nd2")


def test_analyze_mps():

    ret = sp.call([python, "-m", "mps", "analyze_mps", voltage_name, "-o", voltage_dir])
    assert ret == 0


def test_collect_mps():

    sp.call([python, "-m", "mps", "analyze_mps", voltage_name, "-o", voltage_dir])
    sp.call([python, "-m", "mps", "analyze_mps", calcium_name, "-o", calcium_dir])

    ret = sp.call([python, "-m", "mps", "collect_mps", voltage_dir, calcium_dir])
    assert ret == 0


if __name__ == "__main__":
    test_analyze_mps()
    test_collect_mps()
