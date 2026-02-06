# Radial Array Tools

This project uses Python scripts and ArcGIS geoprocessing tools to perform visibility analyses using radial arrays. A radial array is a GIS vector array for which all vectors have the same origin, and are separated by an equal angular interval. These scripts use radial arrays to assess visibility, treating rays in each array as sight-lines for an observer situated at the array's origin point.

These scripts incorporate elements of scripts that I designed during my Masters and PhD research. I have continued to develop them while working for the Maritime Encounters and The Practical Mariner projects. See 'Bibliography' below for relevant links.

This repository is currently being updated to prepare for initial release.

## VIS_NETCDF

This script produces a multidimensional NetCDF file containing angular visibility indices representing the visibility of land as viewed from the sea. For more details, see Smith & Hulin 2026, in the bibliography below.

### Quick Start Guide

To run the visibility model, run *vis_netcdf.py* in the *source* folder of the repository. This will prompt the user to use a default path (which may be altered in the .py file), or
enter a path manually. The path must be to a valid .json file containing model parameters (see *Script Inputs and Outputs* for the parameters list).

The script is currently configured to produce a multidimensional NetCDF file (see bibliography, below). If the input 'out_ncdf' path already exists, the script will ask you to confirm before overwriting it.

### Required Script Inputs

 - *out_ncdf*: The output file produced by the script. This parameter should be a string containing a valid path to a NetCDF output file. This path should use double-backslashes and have the extension '.nc'.


 - *pr_gdb*: The geodatabase used for scratch datasets. This parameter should be a string containing a valid path to an ArcGIS geodatabase. This path should use double-backslashes and have the extension '.gdb'.


 - *in_dem*: The elevation model used to determine whether sight-lines are obstructed by terrain. The vertical units (i.e. metres above sea level) used by this raster must conform to those used elsewhere in the input parameters. This parameter should be a string containing a valid path to georeferenced raster image. This path should use double-backslashes. Currently the script has only been tested with .tif files.


 - *pt_mask_fc*: Features which serve as a geographic constraint on the output (i.e. the 'outer bounds' of the model area). This parameter should be a string containing a valid path to ESRI shapefile. This path should use double-backslashes and have the extension '.shp'.


 - *land_fc*: Features which allow the script to eliminate points on land. This parameter should be a string containing a valid path to ESRI shapefile. This path should use double-backslashes and have the extension '.shp'. Currently the script requires a valid file.


 - *xy_spacing*: Controls the cartesian resolution of the output (its 'XY' resolution). It should be a positive numeric value (currently the script has only been tested with integers).


 - *dist_range*: Controls the 'maximum visible distances' checked for visibility by the script. It should be a list of numeric values (currently the script has only been tested with integers).


 - *densify_dist*: Controls the density at which individual vectors in the radial array will sample the input DEM. Lower values will produce better results, at the cost of increased computation. It should be a list of numeric values (currently the script has only been tested with integers).


 - *obs_z_range*: Controls observer heights used by the script. It should be a list of numeric values (currently the script has only been tested with integers).


### Optional Script Inputs

 - *obs_z_offset*: Adds an offset (in vertical units) to the observer. This is meant to represent the height of an observer above the land/sea surface. It should be a numeric value. By default, it is set to zero.


 - *sample_raster*: WIP


 - *landmark_fc*: WIP


 - *override_dist_list*: WIP


 - *pt_mask_json*: WIP


 - *write_log*: WIP


### Demo

WIP

### Bibliography

WIP
