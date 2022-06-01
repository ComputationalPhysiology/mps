# Binary executables

There are three files in this directory. One for Windows, one for Mac and one for Linux (Ubuntu).
These are currently located in a `.zip` file, so the first thing you need to do is to unzip the binary that you want to use.

The binary is called `mps` (with a `.exe` extension if you are running on Windows)
Note that the binaries ships with python, so to don't even need to install python to run the binaries.
Note, however that the size of the binaries are large.

## Usage

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
./mps analyze file.czi
```

(Windows)
```
.\mps.exe analyze file.czi
```

will analyze that file. You can also type

```
.\mps.exe analyze --help
```
to see all the available options.

### Adding the binary to your path

It is inconvenient to run the binary directly from the extracted folder. Therefore it is recommended that you move this file to some location that is in your [`PATH` environment variable](https://docs.oracle.com/en/database/oracle/r-enterprise/1.5.1/oread/creating-and-modifying-environment-variables-on-windows.html).
For example, if you are on windows you can copy this file to `C:\Windows\System32` directory, and if on linux you could copy it to `usr/bin`.
This will, however require admin access. If you don't have admin access you can put it in a new folder and add that folder to your `PATH` environment variable.
On linux and mac you also need to make the file executable:
```
chmod +x mps
```
