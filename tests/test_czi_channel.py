import numpy as np
import pytest

from unittest.mock import Mock, patch

from mps.load import _format_czi_frames, load_czi


def test_load_czi_with_mocked_czifile():
    num_frames = 5
    y = 20
    x = 30

    # Fake CZI data with shape: (time, y, x)
    images = np.zeros((num_frames, y, x), dtype=np.uint16)
    for t in range(num_frames):
        images[t, :, :] = t

    # Fake metadata used by load_czi to get um_per_pixel.
    # load_czi multiplies this by 1e6, so 1e-6 m becomes 1.0 um.
    metadata = {"Metadata": {"Scaling": {"Items": {"Distance": {"Value": "0.000001"}}}}}

    class FakeTimeStamps:
        def __init__(self, time_stamps):
            self.time_stamps = time_stamps

    timestamp_data = FakeTimeStamps(np.arange(num_frames, dtype=float))

    timestamp_attachment = Mock()
    timestamp_attachment.data_segment.return_value.data.return_value = timestamp_data

    mock_czi = Mock()
    mock_czi.asarray.return_value = images
    mock_czi.metadata = Mock()
    mock_czi.attachment_directory = [timestamp_attachment]

    mock_context = Mock()
    mock_context.__enter__ = Mock(return_value=mock_czi)
    mock_context.__exit__ = Mock(return_value=None)

    with (
        patch("mps.load.czifile.CziFile", return_value=mock_context),
        patch("mps.load.czifile.elem2dict", return_value=metadata),
        patch("mps.load.czifile.TimeStamps", FakeTimeStamps),
    ):
        data = load_czi("fake_file.czi")

    assert data.frames.shape == (x, y, num_frames)
    assert np.all(data.frames[:, :, 0] == 0)
    assert np.all(data.frames[:, :, 1] == 1)
    assert np.all(data.frames[:, :, 4] == 4)

    assert len(data.time_stamps) == num_frames
    assert data.info["size_x"] == x
    assert data.info["size_y"] == y
    assert data.info["um_per_pixel"] == 1.0
    assert data.info["num_frames"] == num_frames


def test_format_czi_frames_standard_3d_time_first():
    num_frames = 5
    y = 20
    x = 30

    # Shape: (time, y, x)
    images = np.zeros((num_frames, y, x))

    for t in range(num_frames):
        images[t, :, :] = t

    frames = _format_czi_frames(
        images=images,
        num_frames=num_frames,
    )

    assert frames.shape == (x, y, num_frames)
    assert np.all(frames[:, :, 0] == 0)
    assert np.all(frames[:, :, 1] == 1)
    assert np.all(frames[:, :, 4] == 4)


def test_format_czi_frames_selects_channel_when_num_frames_equals_num_channels_time_first():
    num_frames = 5
    num_channels = 5
    y = 20
    x = 30

    # Shape: (time, channel, y, x)
    images = np.zeros((num_frames, num_channels, y, x))

    for t in range(num_frames):
        for c in range(num_channels):
            images[t, c, :, :] = 100 * c + t

    frames = _format_czi_frames(
        images,
        num_frames=num_frames,
        channel=2,
    )

    assert frames.shape == (x, y, num_frames)

    # Selected channel = 2, so each time frame should contain 200 + t.
    assert np.all(frames[:, :, 0] == 200)
    assert np.all(frames[:, :, 1] == 201)
    assert np.all(frames[:, :, 4] == 204)


def test_format_czi_frames_selects_channel_when_num_frames_equals_num_channels_time_last():
    num_frames = 5
    num_channels = 5
    y = 20
    x = 30

    # Shape: (y, channel, x, time)
    images = np.zeros((y, num_channels, x, num_frames))

    for t in range(num_frames):
        for c in range(num_channels):
            images[:, c, :, t] = 100 * c + t

    frames = _format_czi_frames(
        images,
        num_frames=num_frames,
        channel=2,
    )

    assert frames.shape == (x, y, num_frames)

    # Selected channel = 2, so each time frame should contain 200 + t.
    assert np.all(frames[:, :, 0] == 200)
    assert np.all(frames[:, :, 1] == 201)
    assert np.all(frames[:, :, 4] == 204)


def test_format_czi_frames_selects_channel_from_4d_data():
    num_frames = 5

    # Fake CZI-like data with shape:
    # (time, channel, y, x)
    #
    # Use spatial dimensions larger than 10, otherwise the channel-axis
    # heuristic may confuse small spatial dimensions with channel axes.
    images = np.zeros((num_frames, 2, 20, 30), dtype=np.uint16)

    images[:, 0, :, :] = 3
    images[:, 1, :, :] = 7

    frames = _format_czi_frames(
        images=images,
        num_frames=num_frames,
        channel=1,
    )

    assert frames.ndim == 3
    assert frames.shape == (30, 20, num_frames)
    assert np.all(frames == 7)


def test_format_czi_frames_requires_channel_for_multichannel_data():
    num_frames = 5
    images = np.zeros((num_frames, 2, 20, 30), dtype=np.uint16)

    with pytest.raises(ValueError, match="Please select one channel"):
        _format_czi_frames(
            images=images,
            num_frames=num_frames,
            channel=None,
        )


def test_format_czi_frames_rejects_invalid_channel():
    num_frames = 5
    images = np.zeros((num_frames, 2, 20, 30), dtype=np.uint16)

    with pytest.raises(ValueError, match="Requested channel 2"):
        _format_czi_frames(
            images=images,
            num_frames=num_frames,
            channel=2,
        )
