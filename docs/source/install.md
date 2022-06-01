# Installation

You can install `cardiac-mps` with pip


```
python -m pip install cardiac-mps
```

## Install from source

If you want the latest version or you want to develop `cardiac-mps` you can install the code on the `master` branch

```
python -m pip install git+https://github.com/ComputationalPhysiology/mps.git@master
```
or clone the repository and install it from there

```
git clone git@github.com:ComputationalPhysiology/mps.git
cd mps
python -m pip install .
```

## Development installation

Developers should use editable install and install the development requirements using the following command
```
python -m pip install -e ".[dev]"
```
It is also recommended to install the `pre-commit` hook that comes with the package
```
pre-commit install
```
Note that linters and formatters will run in the CI system.
