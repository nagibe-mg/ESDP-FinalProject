# Earth System Data Processing: Precipitation Comparison Chain
## Final Project - Earth System Data Processing - WiSE 25/26
**Authors: Johanna Kasishcke and Nagibe Maroun González**

Collection of notebooks and information regarding the final project for Earth System Data Processing. 
Main responsible:
- IMERG data processing: Johanna
- ERA5 data processing: Nagibe


### Overview
The project consists of a robust automated processing chain to compare precipitation data from two different sources: ERA5 dataset and IMERG data. This is relevant because depending on how the data was acquired, it can be different. More specifically, ERA5 deals with modelled data (form a physical model), while IMERG is satellite data. 

The system performs the following key tasks:

Automated Download: Fetches daily batch data with "catch-up" logic to handle gaps.
Spatial Harmonization: Regrids distinct source grids (Lat/Lon) onto a common HEALPix grid (NSIDE 128) to enable direct pixel-to-pixel comparison.
Data Optimization: Processes raw, heavy satellite data into lightweight, region-specific daily aggregates.
Analysis & Visualization: Generates daily comparison maps, calculates inter-annual variability (e.g., 2014 vs. 2015), and produces animated GIFs of weather patterns over Mexico.

### Repository Structure
```
ESDP-FinalProject/
└── ERA5_Data/
    ├── ERA5_Precipitation.ipynb': The main control notebook and visualization interface for data from 2014.
    ├── ERA5_Precipitation2015.ipynb': The main control notebook and visualization interface for data from 2015.
    ├── RealDailyRoutineERA5.py': Contains the core logic for downloading, regridding, saving and some plotting routines.
    ├── Imerg_ERA5_Comparison.ipynb': Comparison for interannual data. (misleading name)
└── GPM/
    ├── data/              # Downloaded datasets
    ├── figures/           # Generated plots and animations
    ├── notebooks/         # Jupyter notebooks for analysis
    │   ├── analyse_stats.ipynb
    │   └── compare_era5_IMERG.ipynb
    ├── scripts/           # Main processing scripts
    │   ├── configfile.py       # Your NASA credentials go here
    │   ├── ControlFlowIMERG.py # Core data processing
    │   ├── runIMERG.py         # Run this to start processing
    │   ├── plot.py             # Create visualizations
    │   └── statistics.py       # Generate statistics
    └── ZARR/              # Processed data storage
        └── IMERG_PRECIP_nside128.zarr
```        

### Requirements
- xarray 
- zarr 
- healpy (not working on Windows)
- cdsapi (with a valid user account)
- matplotlib
- imageio
- netCDF4
- datetime
- json
- imegeio.v2
- gpm-api  &rarr; install via `pip install gpm-api` (there might be dependency issues, but for me it worked until now)
- Possibly: pyresample via `conda install -c conda-forge pyresample` (there might be dependency issues, but for me it worked until now)


### Dataset ERA5
ERA5-Land is a reanalysis dataset providing a consistent view of the evolution of land variables over several decades at an enhanced resolution compared to ERA5. ERA5-Land has been produced by replaying the land component of the ECMWF ERA5 climate reanalysis. Reanalysis combines model data with observations from across the world into a globally complete and consistent dataset using the laws of physics. Reanalysis produces data that goes several decades back in time, providing an accurate description of the climate of the past.
The chosen variable is total precipitation, daily sum. Data is given in m/day and the format is NETCDF and later processed to Zarr.

### Dataset IMERG
IMERG (Integrated Multi-satellitE Retrievals for GPM) provides satellite-based precipitation estimates across most of Earth's surface. It's especially useful for oceans and remote areas without ground-based weather stations.
Learn more: [IMERG Overview](https://gpm.nasa.gov/data/imerg).
We use half-hourly data from GES DISC's "IMERG Final Run" and convert it to daily totals:
```ds = (ds['precipitation'] * 0.5).sum(dim='time')```

### Implementation Details for ERA5 Data

1. Control Flow:
The workflow is managed by a control function that scans a user-defined date range. Before processing a day, the script checks for a marker file. If the marker is missing, it triggers the pipeline: Download → Regrid → Store. This ensures that if the process is interrupted it will resume exactly where it left off.
2. Dataset & Variable Selection: 
The routine is fully flexible in terms of the data parameters. All user settings (dates, variables, levels, and hours) are passed as arguments from the notebook.
3. Spatial Regridding: 
The processing chain converts the raw ERA5 Latitude-Longitude grid into a HEALPix (Hierarchical Equal Area isoLatitude Pixelization) grid (remapped to NSIDE 128) using linear interpolation.
4. Storage Strategy: 
Data is archived in a Zarr store using a group-based hierarchy and chunked by the time dimension (1 day per chunk) to optimize for time-series analysis.
5. Visualization: 
It loads data for the selected spatial domain and plots using a scatter plot function. To ensure scientific comparability, all plots use a fixed color scale (0 to 60 mm/day) and consistent colormaps (blue scale). With these frames, a gif animation was created. 

### Implementation Details for IMERG data

**Helpful resources:**
- tutorial for downloading IMERG data: [IMERG tutorial](https://gpm-api.readthedocs.io/en/latest/tutorials/tutorial_02_IMERG.html)
- set-up of GPM API: [gpm-api](https://gpm-api.readthedocs.io/en/latest/03_quickstart.html)

#### Run the code:
1. To work with the given code, you need to install the following modules preferably in a virtual environment or via conda:
    - `datetime`, `timedelta`
    - `os`
    - `cartopy.crs` 
    - `matplotlib`
    - `xarray`
    - `healpy`
2. The actual code is written in the `ControlFlowIMERG.py` script. This script contains the downloading and processing routine. To execute this code, the file `runIMERG.py`needs to be run. In this file, the timeperiod of interest can be defined.
3. To plot the data and create an animation, you need to run the file `plot.py` in the scripts folder.
    - To plot single day, you need to change day index (day_idx) for different day
    - keep in mind, that python starts indexing with 0
    - images and other graphics are stored in the `figures` folder
4. To create the statistics of the data, you need to run the `statistics.py` script. 
5. To compare both datasets, the ipython-notebook `compare_era5_IMERG.py` needs to be executed.

### Technical Notes:

- Large files are managed using Git LFS
- Data is stored in ZARR format for efficient access
- All outputs are saved to the figures/ folder

### Challenges
* Main issue was the IMERG API. It was difficult to constrain variables such as spatial domain, selected parameters, data frequency, etc. This caused very large download volume. We had to rethink our strategy. 

* Zarr store created by one author could not be opened by the other author by unknown reasons. 

### References
Copernicus Climate Change Service (C3S)(2019): ERA5-Land hourly data from 1950 to present. Copernicus Climate Change Service (C3S) Climate Data Store (CDS). DOI: 10.24381/cds.e2161bac

Huffman, G.J., E.F. Stocker, D.T. Bolvin, E.J. Nelkin, Jackson Tan (2023), GPM IMERG Final Precipitation L3 1 day 0.1 degree x 0.1 degree V07, Edited by Andrey Savtchenko, Greenbelt, MD, Goddard Earth Sciences Data and Information Services Center (GES DISC), Accessed: [04.02.2026], [10.5067/GPM/IMERGDF/DAY/07](https://doi.org/10.5067/GPM/IMERGDF/DAY/07)
