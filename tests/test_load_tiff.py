"""Tests for loading OME-TIFF (e.g. Aurox) and Leica line scan TIFFs.

The files are small synthetic versions mimicking the metadata layout
of the real files.
"""

import csv
import sys

import numpy as np
import pytest
import tifffile

import mps
from mps.utils import file_stem, is_brightfield


def ca_transient(t, period=1000.0, tau1=50.0, tau2=110.0, tstart=50.0):
    s = np.mod(t, period) - tstart
    y = np.where(s > 0, np.exp(-s / tau2) - np.exp(-s / tau1), 0.0)
    return y / y.max()


def read_summary_csv(path):
    """Read the csv written by mps summary as a list of dicts. The
    delimiter is "," on Windows and ";" elsewhere (see utils.to_csv)."""
    delimiter = "," if sys.platform == "win32" else ";"
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def write_ome_tiff(folder, name, frames, delta_t_ms, contrast="Fluorescence"):
    """Write an OME-TIFF that keeps its metadata in a companion file,
    like the Aurox software does."""
    folder.mkdir(parents=True, exist_ok=True)
    companion = f"{name}.companion.ome"
    binary_only = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06">'
        f'<BinaryOnly MetadataFile="{companion}" UUID="urn:uuid:1"/></OME>'
    )
    tiff = folder / f"{name}.ome.tiff"
    tifffile.imwrite(tiff, frames, description=binary_only, metadata=None, photometric="minisblack")

    num_frames, size_y, size_x = frames.shape
    planes = "".join(
        f'<Plane DeltaT="{t}" DeltaTUnit="ms" TheC="0" TheT="{i}" TheZ="0"/>'
        for i, t in enumerate(delta_t_ms)
    )
    (folder / companion).write_text(
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06">'
        '<Image ID="Image:0" Name="pos1">'
        '<Pixels DimensionOrder="XYZCT" ID="Pixels:0" PhysicalSizeX="2.6" '
        'PhysicalSizeY="2.6" SizeC="1" SizeZ="1" '
        f'SizeT="{num_frames}" SizeX="{size_x}" SizeY="{size_y}" '
        'TimeIncrement="0.02" TimeIncrementUnit="ms" Type="uint16">'
        f'<Channel ID="Channel:0:0" Name="Cy5-LED" ContrastMethod="{contrast}"/>'
        f"{planes}</Pixels></Image></OME>",
        encoding="utf-8",
    )
    return tiff


def write_leica_linescan(path, lines, line_time_s=0.0025, finterval_s=None):
    """Write an ImageJ TIFF as exported by Bio-Formats from a Leica .lif line scan."""
    name = "TUS4_DC_BEAT_1"
    other = "TUS4_DC_BEAT_2"
    num_pixels = lines.shape[-1]
    info = "\n".join(
        [
            f"{other} Image #0|ATLConfocalSettingDefinition #0|LineTime = 0.1",
            f"{other} Image #0|ATLConfocalSettingDefinition #0|ScanMode = xt",
            f"{name} Image #0|ATLConfocalSettingDefinition #0|LineTime = {line_time_s}",
            f"{name} Image #0|ATLConfocalSettingDefinition #0|ScanMode = xt",
            f"{name} Image #0|DimensionDescription #2|DimID = 1",
            f"{name} Image #0|DimensionDescription #2|Length = {num_pixels * 0.25e-6}",
            f"{name} Image #0|DimensionDescription #2|NumberOfElements = {num_pixels}",
        ]
    )
    metadata = {"Info": info}
    if lines.ndim == 3:
        metadata["axes"] = "TYX"
        metadata["finterval"] = finterval_s
        metadata["Labels"] = [f"c:1/2 t:{i + 1}/{len(lines)} - {name}" for i in range(len(lines))]
    else:
        metadata["axes"] = "YX"
        metadata["Labels"] = [f"c:1/2 - {name}"]
    tifffile.imwrite(path, lines, imagej=True, metadata=metadata)
    return path


def test_ome_tiff_with_companion(tmp_path):
    rng = np.random.default_rng(1)
    frames = rng.integers(0, 1000, size=(20, 4, 6), dtype=np.uint16)
    delta_t = 1000.0 + np.arange(20) * 20.5
    path = write_ome_tiff(tmp_path, "rec", frames, delta_t)

    data = mps.MPS(path)

    assert data.frames.shape == (6, 4, 20)
    assert np.array_equal(data.frames[:, :, 3], frames[3].T)
    assert np.allclose(data.time_stamps, delta_t - delta_t[0])
    assert data.info["dt"] == pytest.approx(20.5)
    assert data.info["um_per_pixel"] == 2.6
    assert data.metadata["contrast_method"] == "Fluorescence"


def test_ome_tiff_missing_companion(tmp_path):
    frames = np.zeros((3, 4, 6), dtype=np.uint16)
    path = write_ome_tiff(tmp_path, "rec", frames, [0, 20, 40])
    path.with_name("rec.companion.ome").unlink()
    with pytest.raises(FileNotFoundError, match="companion.ome"):
        mps.MPS(path)


def test_leica_linescan_repeats_are_concatenated(tmp_path):
    lines = np.arange(3 * 16 * 5, dtype=np.uint8).reshape(3, 16, 5)
    path = write_leica_linescan(tmp_path / "scan.tif", lines, finterval_s=0.05)

    data = mps.MPS(path, linescan_bin=1)

    assert data.frames.shape == (5, 1, 48)
    assert np.array_equal(data.frames[:, 0, 16 + 2], lines[1, 2])
    expected_t = (np.arange(3)[:, None] * 50.0 + np.arange(16)[None, :] * 2.5).ravel()
    assert np.allclose(data.time_stamps, expected_t)
    assert data.info["um_per_pixel"] == pytest.approx(0.25)
    assert data.info["line_time_ms"] == pytest.approx(2.5)
    assert data.metadata["format"] == "leica-linescan"


def test_leica_linescan_binning(tmp_path):
    lines = np.arange(2 * 16 * 5, dtype=np.uint8).reshape(2, 16, 5)
    path = write_leica_linescan(tmp_path / "scan.tif", lines, finterval_s=0.05)

    data = mps.MPS(path, linescan_bin=4)

    assert data.frames.shape == (5, 1, 8)
    assert np.allclose(data.frames[:, 0, 5], lines[1, 4:8].mean(axis=0))
    assert np.allclose(data.time_stamps[:4], [3.75, 13.75, 23.75, 33.75])
    assert data.time_stamps[4] == pytest.approx(50 + 3.75)
    assert data.info["linescan_bin"] == 4


def test_leica_linescan_single_repeat(tmp_path):
    lines = np.zeros((16, 5), dtype=np.uint8)
    path = write_leica_linescan(tmp_path / "scan.tif", lines)
    data = mps.MPS(path, linescan_bin=1)
    assert data.frames.shape == (5, 1, 16)
    assert np.allclose(data.time_stamps, np.arange(16) * 2.5)


def test_leica_linescan_overlapping_repeats_raises(tmp_path):
    lines = np.zeros((2, 16, 5), dtype=np.uint8)
    path = write_leica_linescan(tmp_path / "scan.tif", lines, finterval_s=0.01)
    with pytest.raises(ValueError, match="Repeat interval"):
        mps.MPS(path)


@pytest.mark.parametrize(
    "name, expected",
    [
        ("E10_bf_10x_15s_50fps.ome.tiff", True),
        ("Control_BF_1.nd2", True),
        ("E10_berst_10x_15s_50fps.ome.tiff", False),
        ("tus4_ogb_berst_WELL2.lif - TUS4_DC_BEAT_1 - C=0.tif", False),
    ],
)
def test_is_brightfield(name, expected):
    assert is_brightfield(name) is expected


def test_file_stem(tmp_path):
    assert file_stem(tmp_path / "rec.ome.tiff") == "rec"
    assert file_stem(tmp_path / "rec.tif") == "rec"


def test_summary_recursive_ome_tiff(tmp_path):
    t = np.arange(0, 6000, 20.0)
    trace = 100 + 50 * ca_transient(t)
    frames = np.tile(trace[:, None, None], (1, 4, 6)).astype(np.uint16)
    write_ome_tiff(tmp_path / "rec1", "E10_berst_rec1", frames, t)
    write_ome_tiff(tmp_path / "rec2", "E10_bf_rec2", frames, t, contrast="Brightfield")

    mps.scripts.summary.main(tmp_path, recursive=True)

    rows = read_summary_csv(tmp_path / "mps_summary.csv")
    assert len(rows) == 1
    assert rows[0]["Filename"] == "rec1/E10_berst_rec1.ome.tiff"
    assert float(rows[0]["Number of beats"]) == 6


def test_analyze_folder_twice_with_unchoppable_file(tmp_path):
    t = np.arange(0, 6000, 20.0)
    beating = np.tile((100 + 50 * ca_transient(t))[:, None, None], (1, 4, 6))
    # Beats longer than max_window -> no beats are found
    slow_beat = ca_transient(t, period=3000, tau1=800, tau2=1200)
    slow = np.tile((100 + 50 * slow_beat)[:, None, None], (1, 4, 6))
    write_ome_tiff(tmp_path / "rec1", "rec1", beating.astype(np.uint16), t)
    write_ome_tiff(tmp_path / "rec2", "rec2", slow.astype(np.uint16), t)

    # Run twice: the second run must not pick up the output (data.npy)
    # from the first run as input
    for _ in range(2):
        mps.scripts.analyze.main(tmp_path, plot=False, max_window=1500)

    assert sorted(p.relative_to(tmp_path).as_posix() for p in tmp_path.rglob("settings.json")) == [
        "rec1/rec1/settings.json",
        "rec2/rec2/settings.json",
    ]
    assert not list(tmp_path.rglob("data/data.npy"))


def test_summary_skips_files_without_beats_and_analyze_output(tmp_path):
    t = np.arange(0, 6000, 20.0)
    beating = np.tile((100 + 50 * ca_transient(t))[:, None, None], (1, 4, 6))
    slow_beat = ca_transient(t, period=3000, tau1=800, tau2=1200)
    slow = np.tile((100 + 50 * slow_beat)[:, None, None], (1, 4, 6))
    write_ome_tiff(tmp_path / "rec1", "rec1", beating.astype(np.uint16), t)
    write_ome_tiff(tmp_path / "rec2", "rec2", slow.astype(np.uint16), t)
    # Output from mps analyze (contains data.npy) must not be used as input
    mps.scripts.analyze.main(tmp_path, plot=False)

    mps.scripts.summary.main(tmp_path, recursive=True, max_window=1500)

    rows = read_summary_csv(tmp_path / "mps_summary.csv")
    assert [r["Filename"] for r in rows] == ["rec1/rec1.ome.tiff"]
    assert (tmp_path / "mps_summary.pdf").is_file()


def test_summary_no_file_with_beats(tmp_path):
    t = np.arange(0, 6000, 20.0)
    slow_beat = ca_transient(t, period=3000, tau1=800, tau2=1200)
    slow = np.tile((100 + 50 * slow_beat)[:, None, None], (1, 4, 6))
    write_ome_tiff(tmp_path, "rec", slow.astype(np.uint16), t)

    mps.scripts.summary.main(tmp_path, max_window=1500)

    assert not (tmp_path / "mps_summary.pdf").exists()
