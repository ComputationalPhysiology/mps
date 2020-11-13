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

## If you find a bug?

If the scripts behave in an unexpected manner, or you encounter any bugs, please submit this as an issue.
Click on the issue tab (on the top of this page). Write a descriptive title and paste in the output from your console.


## Contact

This software is developed by Henrik Finsberg at Simula Research Laboratory.
If you need to get in contact, please send me an email at [henriknf@simula.no](mailto:henriknf@simula.no).

## License

c) 2001-2020 Simula Research Laboratory ALL RIGHTS RESERVED

END-USER LICENSE AGREEMENT
PLEASE READ THIS DOCUMENT CAREFULLY. By installing or using this
software you agree with the terms and conditions of this license
agreement. If you do not accept the terms of this license agreement
you may not install or use this software.

Permission to use, copy, modify and distribute any part of this
software for non-profit educational and research purposes, without
fee, and without a written agreement is hereby granted, provided
that the above copyright notice, and this license agreement in its
entirety appear in all copies. Those desiring to use this software
for commercial purposes should contact Simula Research Laboratory AS:
post@simula.no

IN NO EVENT SHALL SIMULA RESEARCH LABORATORY BE LIABLE TO ANY PARTY
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE
"MPS" EVEN IF SIMULA RESEARCH LABORATORY HAS BEEN ADVISED
OF THE POSSIBILITY OF SUCH DAMAGE. THE SOFTWARE PROVIDED HEREIN IS
ON AN "AS IS" BASIS, AND SIMULA RESEARCH LABORATORY HAS NO OBLIGATION
TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
SIMULA RESEARCH LABORATORY MAKES NO REPRESENTATIONS AND EXTENDS NO
WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESSED, INCLUDING, BUT
NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS
