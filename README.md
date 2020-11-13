# MPS data

This repository contains the software for reading and analyzing MPS data.

## Instructions for regular users

### Installation

If you only want to use this software (and not develop it) then you can download the binaries

- [Windows binaries](https://github.com/finsberg/mps/suites/1493516025/artifacts/25997813)
- [MacOSX binaries](https://github.com/finsberg/mps/suites/1493516025/artifacts/25997811)
- [Linux binaries](https://github.com/finsberg/mps/suites/1493516025/artifacts/25997812)

The binary is called `mps` (with a `.exe` extension if you are running on Windows)
Note that the binaries ships with python, so to don't even need to install python to run the binaries.
Note, however that the size of the binaries are large.
In order to run the binaries you just need to download the zip file, and extract the executable to a location that your operating system can find it.
For example, if you are on windows you can copy this file to `C:\Windows\System32` directory, and if on linux you could copy it to `usr/bin`.
This will, however require admin access. If you don't have admin access you can put it in a new folder and add that folder to your `PATH` environment variable.
On unix you also need to make the file executable:
```
chmod +x mps
```
You can also run the file from its current location.

### Usage

If you execute the binary without any arguments, i.e

on Unix
```
./mps
```
or on Windows
```
.\mps.exe
```

you should get the following output

```
MPS package for analyzing microphyisological systems

Available arguments
-------------------
All these arguments can be called with '-h' or '--help' to see the
additional options

    analyze
        Analyze mps data (.nd2 or .czi)

    summary
        Create a summary figure and csv file of all files in a folder

    phase_plot
        Make a phase plot with voltage on the x-axis and calcium on the y-axis.

    prevalence
        Estimate the percentage of tissue in the chips,
        and the percentage of beating tissue vs non-beating

    collect
        Gather Voltage and Calcium data into one file that can be
        used as input to the inversion algorithm.

    mps2mp4
        Create movie of data file


Available options
-----------------

    -h, --help      Show this help
    -v, --version   Show version number


Contact
-------
Henrik Finsberg (henriknf@simula.no)
```

This shows the available commands that you can use with the script.
For example if yoo have an mps data file called `file.nd2` then running

(Unix)
```
./mps analyze file.nd2
```

(Windows)
```
.\mps.exe analyze file.nd2
```

will analyze that file. You can also type

```
.\mps.exe analyze --help
```
to see all the available options.


### Known issues

1. The `prevalence` script is not yet working

2. To use the `mps2mp4` script you also need to install `ffmpeg` separatly


## Instructions for developers

TBW

### Installations instructions

TBW


## Documentation

If you run

```
make docs
```

you will generate documentation that can be viewed in the browser.
