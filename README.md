# MPS data

This repository contains the code for reading and analysing the MPS data.

## Installation

### Install python
In order to run this code you need python version 3.6 or above.
There are several ways to install python. For Windows users I
recommend using the official python distribution as
[https://www.python.org/downloads/](https://www.python.org/downloads/). It
is also possible to use
[anaconda](https://docs.anaconda.com/anaconda/install/) which works on
all platform. On mac it is also possible to use
[Homebrew](https://brew.sh) and on ubuntu you can also use the [Ubuntu
packaging
manager](https://packages.ubuntu.com/search?suite=default&section=all&arch=any&keywords=python3&searchon=names).


### Verification of python installation
When you have installed python you should check that python was
successfully installed and that you are able to run it from the
terminal. To do that you must first open a terminal. If you are on Mac
click on the rocket and type `terminal`. If you are on Windows click on
the start meny and type `cmd` or `PowerShell` (PowerShell is
recommended). If you are on Linux then you should know how to open a
terminal. Next try to type `python` or `python3` in the terminal. If
that works you are all set. If not then you need to update your `PATH`
environment variable to that `python` can be accessed. See
instructions for
[Windows](https://docs.python.org/3/using/windows.html#excursus-setting-environment-variables).



### Create a virtual environment
It is always good practice to install all your dependencies in a
virtual environment. This makes it possible to keep all the packages in the right versions.

#### Using Virtualenv
Another way (which I will be using) is to create a virtual environment
using pythons `virtualenv` package. Install this with pip
```
pip install virtualenv
```
If you don't have pip you download the script
[`get-pip.py`](https://bootstrap.pypa.io/get-pip.py) (save it on you computer) and run

```
python get-pip.py
```
Now you can create a virtual environment using virtualenv by typing
```
python -m virtualenv venv
```
This will create a new folder called `venv` where all your packages will
be installed. (Run this command somewhere that you would like to put
this folder). Now you want to activate the virtual environment.
On Mac / Linux do
```
source venv/bin/activate
```
or on Windows do
```
.\venv\Scripts\activate
```
Now you are ready to install the software.

#### Using Conda
 If you are using `conda` you can use [conda environments](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).
In the conda environment, the commands to perform are quite similar.
In order to create the new environment type:
```
conda create -n yourenvname python=3.7 anaconda
```
Once the environment is created one can activate it in any folder of your personal machine
by prompting:
```
source activate yourenvname
```

### Install the dependencies
Before we install the actual software we will install some packages
that the software depends upon. You can do this by typing
```
pip install -r requirements.txt
```

#### Problems on Windows?
If you are using windows and get problems with compatibility of the
packages then I found out that a good solution is to install the
dependencies using something called `pipwin`. Then you can do
```
pip install pipwin
```
and then do
```
pipwin install -r requirements.txt
```
This will install some prebuilt windows binaries that are hosted at
[https://www.lfd.uci.edu/~gohlke/pythonlibs/](https://www.lfd.uci.edu/~gohlke/pythonlibs/).


### Install the software

Finally you can install the software by typing
```
python setup.py install
```
or (equivalently)
```
pip install .
```


## Usage

While it is possible to use the `mps` software through the python
interface, there is also a command-line interface. You can run all the
scripts from the terminal (GUIs will come later). There are three
main scripts


* analyze_mps

	- type `analyze_mps --help` for more info


* analyze_mps_all

	- type `analyze_mps_all --help` for more info

* mps_summary

	- type `mps_summary --help` for more info


## Demos

Demos can be found in the folder ``demo``.
`
## Automated testing

There are tests located in the `test` folder that are run every time
there is a change in the repository.

## Documentation

For more documentation you can build the documentation in html by
doing
```
cd doc
make html
```

Now you can open the file `doc/build/html/index.html` in your favorite web browser to see the documentation.

## Maintainer

Henrik Finsberg (henriknf@simula.no)
