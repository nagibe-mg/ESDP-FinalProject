## Project Proposal: Comparative Analysis of ERA5-Land and GPM-IMERG Precipitation (Mexico, May 2014)

##### Members: Johanna Kasischke, Nagibe Maroun-Gonzalez



**Task Description and Objectives**

The goal of this project is to develop a scalable data processing pipeline to compare two leading global precipitation datasets: ERA5-Land (reanalysis) and GPM-IMERG (satellite-derived). We will focus on Mexico during May 2014 as a case study to evaluate how these products capture regional precipitation patterns.



**Core Objectives:**

* Automated Acquisition: Implement a script to download ERA5-Land and IMERG data for a specific bounding box and time range.
* Data Harmonization: Address the differing resolutions—ERA5-Land (0.1°×0.1°) and IMERG (0.1°×0.1° but on a different grid) through regridding and coordinate matching.
* Statistical Evaluation
* Scalability: Ensure the tool can be easily adjusted for different regions or timeframes by modifying a single configuration file.
* Parallelization: processing of a big dataset can be expensive and slow. Parallelization of this can make it more efficient. 



**Tasks** 

* data Engineering \& Pipeline- Set up the cdsapi for ERA5-Land and NASA Earthdata scripts for IMERG
* handle file format conversions (GRIB/HDF5 to NetCDF)
* implement spatial subsetting
* analysis \& visualization- Develop the regridding logic
* compute statistical metrics
* create the final comparative maps and time-series plots
* contribute to the final README documentation and presentation



**Implementation Plan \& Methodology**

Download: Use Python to fetch ERA5-Land precipitation and IMERG "Final Run".

Processing: Align temporal resolutions by aggregating both to daily totals. Potentially use parallelization to make this more efficient. 

Regrid IMERG data to match the ERA5-Land grid using interpolation to ensure comparison.

Analysis: Generate spatial difference maps (ERA5 minus IMERG) to identify where the models diverge most.

Documentation: Use Git for version control, maintaining a clear history of implementation choices.

