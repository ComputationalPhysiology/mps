# Command line interface - `mps-summary`

If you have a folder containing several imaging files, then you can run `mps-summary` to get one figure with all the traces and one spreadsheet with some summary statistics.


Say that we have a folder called `folder` containing several imaging files. Then you can run the command
```
mps-summary folder
```
or
```
python -m mps summary folder
```
and it will produce two new files called `mps_summary.csv` and `mps_summary.pdf` inside the folder.

You can also list all the possible option using the command
```
mps-summary --help
```
