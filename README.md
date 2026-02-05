# Earth System Data Processing: Precipitation Comparison Chain
## Final Project - Earth System Data Processing - WiSE 25/26
**Authors: Johanna Kasishcke and Nagibe Maroun González**

Collection of notebooks and information regarding the final project for Earth System Data Processing.

### Overview
The project consists of a robust automated processing chain to compare precipitation data from two different sources: ERA5 dataset and IMERG data. This is relevant because depending on how the data was acquired, it can be different. More specifically, ERA5 deals with modelled data (form a physical model), while IMERG is satellite data. 

The system performs the following key tasks:

Automated Download: Fetches daily batch data with "catch-up" logic to handle gaps.
Spatial Harmonization: Regrids distinct source grids (Lat/Lon) onto a common HEALPix grid (NSIDE 128) to enable direct pixel-to-pixel comparison.
Data Optimization: Processes raw, heavy satellite data into lightweight, region-specific daily aggregates.
Analysis & Visualization: Generates daily comparison maps, calculates inter-annual variability (e.g., 2014 vs. 2015), and produces animated GIFs of weather patterns over Mexico.

### Repository Structure
- 'ERA5_Data/ERA5_Precipitation.ipynb': The main control notebook and visualization interface for data from 2014.
- 'ERA5_Data/ERA5_Precipitation2015.ipynb': The main control notebook and visualization interface for data from 2015.
- 'ERA5_Data/RealDailyRoutineERA5.py': Contains the core logic for downloading, regridding, saving and some plotting routines.
- 'ERA5_Data/Imerg_ERA5_Comparison.ipynb': Comparison for interannual data. (misleading name)

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


### Dataset ERA5
ERA5-Land is a reanalysis dataset providing a consistent view of the evolution of land variables over several decades at an enhanced resolution compared to ERA5. ERA5-Land has been produced by replaying the land component of the ECMWF ERA5 climate reanalysis. Reanalysis combines model data with observations from across the world into a globally complete and consistent dataset using the laws of physics. Reanalysis produces data that goes several decades back in time, providing an accurate description of the climate of the past.
The chosen variable is total precipitation, daily sum. Data is given in m/day and the format is NETCDF and later processed to Zarr.

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



### Challenges
* Main issue was the IMERG API. It was difficult to constrain variables such as spatial domain, selected parameters, data frequency, etc. This caused very large download volume. We had to rethink our strategy. 

* Zarr store created by one author could not be opened by the other author by unknown reasons. 

### References
Copernicus Climate Change Service (C3S)(2019): ERA5-Land hourly data from 1950 to present. Copernicus Climate Change Service (C3S) Climate Data Store (CDS). DOI: 10.24381/cds.e2161bac
GPM IMERG: Huffman, G. et al. (2019). GPM IMERG Final Precipitation L3 1 day 0.1 degree x 0.1 degree V06, Greenbelt, MD, Goddard Earth Sciences Data and Information Services Center (GES DISC)
