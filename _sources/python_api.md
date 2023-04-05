# Python API


## Reading a file
The most useful features of this package it to read imaging data which can be done as follow

```python
import mps
# Object containing the frames, time stamps and metadata
data = mps.MPS("file.nd2")
```

The package supports numerous types of formats including `.nd2`, `.czi` and `.tif`.

## Accessing info about the file

You can print some general info about the dataset by accessing the `.info` attribute

```python
>>> print(data.info)
{'num_frames': 591, 'dt': 10.007577876721394, 'time_unit': 'ms', 'um_per_pixel': 2.589535590385804, 'size_x': 415, 'size_y': 154}
```

Most of fields here are self-explanatory. `num_frames` is the number of frames in the image stack, `dt` is the average time difference between two successive frames, `time_unit` is the unit of time (either `'ms'` or `'s'`), `um_per_pixel` is a factor to convert from pixel unit to micro meters. Finally `size_x` and `size_y` is the number of pixels in the x- and y direction.


There is also another attribute called `.metadata` that contains all the metadata that is possible to extract from the file

```python
>>> print(data.metadata)

{'ImageTextInfoLV': {'SLxImageTextInfo': {'TextInfoItem_0': '', 'TextInfoItem_1': '', 'TextInfoItem_2': '', 'TextInfoItem_3': '', 'TextInfoItem_4': '', 'TextInfoItem_5': 'Metadata:\r\nDimensions: T(591) x λ(1)\r\nCamera Name: Flash4.0, SN:001336\r\nNumerical Aperture: 0.3\r\nRefractive Index: 1\r\n Name: Red_VC\r\n Component Count: 1\r\n Modality: Widefield Fluorescence\r\n Camera Settings:   Exposure: 10 ms\r\n  Binning: 4x4\r\n  Scan Mode: Fast\r\n Microscope Settings:
...
```

## Imaging data

`data.frames` gives you access to the frames which are numpy arrays. You can verify that the shape of the frames matches the sizes from the info dictionary
```python
>>> print(data.frames.shape)
(415, 154, 591)
```
You can also check the data type
```python
>>> print(data.frames.dtype)
dtype('uint16')
```

The time stamps are stored as an attribute called `.time_stamps`
```python
>>> print(data.time_stamps.shape)
(591, )
```

Finally, in some cases pacing amplitude is also available in the metadata, and this can be accessed from the attribute `.pacing`

```python
>>> print(data.pacing.shape)
(591, )
```

If pacing information is not available then `data.pacing` will only contain zeros.



## Analyzing imaging data

To convert the data into a trace you can use the library [`ap_features`](https://computationalphysiology.github.io/ap_features)

```python
import ap_features as apf

# Compute the average over all pixels
y = mps.average.get_average_all(data.frames)
# Convert it to an apf.Beats and compute features
# using the ap_features package
trace = apf.Beats(y=y, t=data.time_stamps, pacing=data.pacing)
```

See the documentation of `ap_features` to learn more about how to analyze a trace.
