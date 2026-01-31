# ESDP-FinalProject

Earth System Data Processing at the University of Cologne. Final project for Winter Semester 2025/2026. Collaborators: Johanna Kasischke, Nagibe Maroun González.



## Load IMERG data

### Folder structure
```
├── ESDP-FinalProject
│   └── GPM
│       └── RS
│   └── configfile.py
│   └── load_IMERG-data.ipynb
```

I followed the following tutorial for downloading the IMERG data:
[IMERG tutorial](https://gpm-api.readthedocs.io/en/latest/tutorials/tutorial_02_IMERG.html)

1. create an NASA PPS or NASA Earthdata Account
2. put username and passwort in the code in the `configfile.py`file (! I filled in my login information here. But we need to ask the prof how we work with it when we submit the homework, because I don't want to share my login information with him. You can still use it this way, because it doesn't make sense, to change it everytime we start the script.)
More information can be found here [gpm-api](https://gpm-api.readthedocs.io/en/latest/03_quickstart.html). 
3. to download and process the data and make a first plot, you need to run the `load_IMERG-data.ipynb` notebook. therefore you need the following modules:
    - `datetime`
    - `cartopy.crs` 
    - `matplotlib.pyplot`
    - `xarray as xr`
    - `ximage` &rarr; install via `conda install ximage`
    - `gpm-api` &rarr; install via `pip install gpm-api` (there might be dependency issues, but for me it worked until now)
    -  `pyresample` &rarr; via `conda install -c conda-forge pyresample` (there might be dependency issues, but for me it worked until now)







