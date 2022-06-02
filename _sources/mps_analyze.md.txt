# Command line interface - `mps-analyze`

If you have an imaging file you might want to run a number of different analysis operations on it. The `mps-analyze` script takes as input an imaging file and output a folder with figures and spreadsheets.

Say that we have a file `file.nd2`. Then we can analyze it using the command
```
mps-analyze file.nd2
```

The script will display some info in the terminal

```
2022-06-01 15:31:34,135 - mps.scripts.analyze - INFO - Run analysis script
2022-06-01 15:31:34,138 - mps.load - INFO - Load nd2 file /Users/henriknf/file.nd2
2022-06-01 15:31:34,309 - mps.load - INFO - Loaded 591 frames in 0.17 seconds
2022-06-01 15:31:37,174 - mps.scripts.analyze - INFO - Finished analyzing MPS data. Data stored folder 'file'.
Total elapsed time: 3.04 seconds
```
and a new folder called `file` contains the following files

```
file
├── EAD_analysis.png
├── README.txt
├── apd_analysis.png
├── average.png
├── average_pacing.png
├── background.png
├── chopped_data.csv
├── chopped_data.png
├── chopped_data_aligned.png
├── data.npy
├── data.txt
├── metadata.json
├── original.png
├── original_pacing.png
├── settings.json
├── sliced_filtered.png
├── trace.png
└── unchopped_data.csv
```

The file called `README.txt` contains a description of all the files in the folder.
For easy reference we also provide this description below

- `apd_analysis.png`
    Figure showing plots of action potential durations
    and corrected actions potential durations (using
    the friderica correction formula). Top panel shows
    bars of the APD for the different beats where the *
    indicated if the APD for that beat is significantly different.
    The middle panel show the APD for each each plotted as a
    line as well as a linear fit of the APD80 and cAPD80 line.
    The intention behind this panel is to see it there is any
    correlation between the APD and the beat number
    The lower panel shows where the cut is made to compute
    the different APDs
- `average_pacing.png`
    These are the average trace with pacing
- `average.png`
    These are the average of the traces in chopped_data
- `background.png`
    This plots the original trace and the background that we
    subtract in the corrected trace.
- `chopped_data_aligned.png`
    All beats plotted on with the time starting at zero
- `chopped_data.csv`
    A csv file with the chopped data
- `chopped_data.png`
    Left panel: All the beats that we were able to extract from the corrected trace
    Right panel: The intersection of all beats where APD80 and APD30 are within 1
    standard deviation of the mean.
- `data.npy`
    A file containing all the results. You can load the results
    in python as follows
    ```python
    >> import numpy as np
    >> data = np.load('data.npy', allow_pickle=True).item()
    ```
    data is now a python dictionary.
- `data.txt`
    This contains a short summary of analysis.
- `EAD_analysis.png`
    A plot of the beats containing EADs
- `frequency_analysis.png`
    A plot of the frequencies of the different beats computed
    as the time between peaks.
- `metadata.json`
    Metadata stored inside the raw data
- `original_pacing.png`
    This is the the raw trace obtained after averaging the frames
    in the stack without any background correction or filtering
    together with the pacing amplitude.
- `original.png`
    This is the the raw trace obtained after averaging the frames
    in the stack without any background correction or filtering
- `sliced_filtered.png`
    Original trace where we have performed slicing and filtering
- `trace.png`
    The is the corrected version of the original trace where we
    have performed a background correction and filtering
- `unchopped_data.csv`
    A csv file with all the unchopped data (original and corrected)
