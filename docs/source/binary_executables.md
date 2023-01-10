# Pre-built binaries

You can also get pre-built binaries for Windows, MaxOSX and Linux in case you don't want to install anything.

- [Windows binaries](https://github.com/ComputationalPhysiology/mps/suites/7730035079/artifacts/324042066)
- [MacOSX binaries](https://github.com/ComputationalPhysiology/mps/suites/7730035079/artifacts/324042063)
- [Linux (ubuntu) binaries](https://github.com/ComputationalPhysiology/mps/suites/7730035079/artifacts/324042065)


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
./mps --help
Usage: mps [OPTIONS] COMMAND [ARGS]...

Options:
  --version                       Show version
  --license                       Show license
  --install-completion [bash|zsh|fish|powershell|pwsh]
                                  Install completion for the specified shell.
  --show-completion [bash|zsh|fish|powershell|pwsh]
                                  Show completion for the specified shell, to
                                  copy it or customize the installation.
  --help                          Show this message and exit.

Commands:
  analyze       Analyze flourecense data
  motion        Estimate motion in stack of images
  mps2mp4       Create movie of data file
  phase-plot    Create movie of data file
  split-pacing  Run script on a folder with files and this will copy the...
  summary       Create a summary pdf of all files in the a directory.
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



## Known issues

1. The `prevalence` script is not yet working

2. To use the `mps2mp4` script you also need to install `ffmpeg` separately
