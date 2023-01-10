# Command line interface - `mps-motion`

## Integration with motion tracking

If you also have installed the [motion tracking script](https://github.com/ComputationalPhysiology/mps_motion) should should be able to also run

```
python -m mps motion file.nd2
```

or simply

```
mps-motion file.nd2
```
in order to run the motion analysis. As before you can do
```
mps-motion --help
```
to see all the available options. The script will produce a folder called `motion` containing the following files
```
motion
├── features.csv
├── results.csv
├── settings.json
├── u_norm.png
└── v_norm.png
```
Here
- `settings.json` are the settings used for the motion analysis
- `u_norm.png` is a plot of the displacement norm averaged over the entire chip as a function of time.
- `v_norm.png` is a plot of the velocity norm averaged over the entire chip as a function of time.
- `results.csv` is a spreadsheet containing the traces
- `features.csv` is a spreadsheet containing some features of each individual trace.
