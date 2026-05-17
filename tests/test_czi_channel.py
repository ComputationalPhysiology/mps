import numpy as np
import pytest

from mps.load import _format_czi_frames


def test_format_czi_frames_selects_channel_from_4d_data():
    num_frames = 5

    # Fake CZI-like data with shape:
    # (time, channel, x, y)
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
    assert frames.shape[-1] == num_frames
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
