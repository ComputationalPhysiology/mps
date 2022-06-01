[![CI](https://github.com/ComputationalPhysiology/mps/actions/workflows/main.yml/badge.svg)](https://github.com/ComputationalPhysiology/mps/actions/workflows/main.yml)
[![PyPI version](https://badge.fury.io/py/cardiac-mps.svg)](https://badge.fury.io/py/cardiac-mps)
[![codecov](https://codecov.io/gh/ComputationalPhysiology/mps/branch/master/graph/badge.svg?token=V5DOQ1PUVF)](https://codecov.io/gh/ComputationalPhysiology/mps)
[![github pages](https://github.com/ComputationalPhysiology/mps/actions/workflows/github-pages.yml/badge.svg)](https://github.com/ComputationalPhysiology/mps/actions/workflows/github-pages.yml)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v2.1-blue.svg)](https://www.gnu.org/licenses/lgpl-2.1)


# MPS data

This repository contains the software for reading and analyzing MPS data.
The analysis scripts are heavily based on [ap_features](https://github.com/ComputationalPhysiology/ap_features) which is a package for computing features of action potential traces.

## Installation

To install the package you can simply di
```
python -m pip install cardiac-mps
```

For more options see <https://computationalphysiology.github.io/mps/install.html>

## Usage


### CLI
Once installed you can use the command line script from the terminal.

In the following demo I have the `mps` packaged installed in a [virtual environment](https://realpython.com/python-virtual-environments-a-primer/) which is first activated. Then it shows how to run the analysis script on a single file.

![_](https://github.com/ComputationalPhysiology/mps/blob/master/docs/source/_static/usage.gif)

### Python API
The most useful features of this package it to read imaging data which can be done as follow

```python
import mps
# Object containing the frames, time stamps and metadata
data = mps.MPS("file.nd2")
```
To convert the data into a trace you can do
```python
import ap_features as apf

# Compute the average over all pixels
y = mps.average.get_average_all(data.frames)
# Convert it to an apf.Beats and compute features
# using the ap_features package
trace = apf.Beats(y=y, t=data.time_stamps, pacing=data.pacing)
```



## Documentation

https://computationalphysiology.github.io/mps

## Contributing

See <https://computationalphysiology.github.io/mps/CONTRIBUTING.html>


## Contact

This software is developed by Henrik Finsberg at Simula Research Laboratory.
If you need to get in contact, please send me an email at [henriknf@simula.no](mailto:henriknf@simula.no).

## License

LGPLv2
