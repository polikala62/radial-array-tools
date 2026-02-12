# Radial Array Tools



*Karl Smith*

*Readme last updated February 2026*





This project uses Python scripts and ArcGIS geoprocessing tools to perform visibility analyses using radial arrays. A radial array is a GIS vector array for which all vectors have the same origin, and are separated by an equal angular interval. These scripts use radial arrays to assess visibility, treating rays in each array as sight-lines for an observer situated at the array's origin point. They incorporate elements of scripts that I designed during my Masters and PhD research. I have continued to develop them while working for the Maritime Encounters and The Practical Mariner projects. See 'Bibliography' below for relevant links. See also the ArcGIS story map "[Making Landfall in the Ancient World](https://storymaps.arcgis.com/stories/e150f03929c049a6b9cce3ba73d98ce3)", which describes outputs from an earlier version of the script.




There is currently one fully-functional, (mostly) commented script in this 'toolbox', called 'vis\_netcdf'. It uses functions in 'radial-array-tools\\source\\functions'. Not all functions in 'radial-array-tools\\source\\functions' are used by 'vis\_netcdf' - many have been used in other projects. If you come across one that looks useful, let me know and I will consider making a dedicated tool for it!





## 1. The *vis\_netcdf* Script

This script produces a multidimensional NetCDF file containing angular visibility indices representing the visibility of land as viewed from the sea. For more details, see Smith \& Hulin (forthcoming), in the bibliography below.



#### 1.1 *vis\_netcdf* Quick Start Guide

To run the visibility model, run *vis\_netcdf.py* in the *source* folder of the repository. This will prompt the user to use a default path (which may be altered in the .py file), or
enter a path manually. The path must be to a valid .json file containing model parameters (see *Script Inputs and Outputs* for the parameters list).

The script is currently configured to produce a multidimensional NetCDF file (see bibliography, below). If the input 'out\_ncdf' path already exists, the script will ask you to confirm before overwriting it.



#### 1.2 *vis\_netcdf* Required Script Inputs



* *out\_ncdf*: The output file produced by the script. This parameter should be a string containing a valid path to a NetCDF output file. This path should use double-backslashes and have the extension '.nc'.



* *pr\_gdb*: The geodatabase used for scratch datasets. This parameter should be a string containing a valid path to an ArcGIS geodatabase. This path should use double-backslashes and have the extension '.gdb'.



* *in\_dem*: The elevation model used to determine whether sight-lines are obstructed by terrain. The vertical units (i.e. metres above sea level) used by this raster must conform to those used elsewhere in the input parameters. This parameter should be a string containing a valid path to georeferenced raster image. This path should use double-backslashes. Currently the script has only been tested with .tif files.



* *pt\_mask\_fc*: Features which serve as a geographic constraint on the output (i.e. the 'outer bounds' of the model area). This parameter should be a string containing a valid path to ESRI shapefile. This path should use double-backslashes and have the extension '.shp'.



* *land\_fc*: Features which allow the script to eliminate points on land. This parameter should be a string containing a valid path to ESRI shapefile. This path should use double-backslashes and have the extension '.shp'. Currently the script requires a valid file.



* *xy\_spacing*: Controls the cartesian resolution of the output (its 'XY' resolution). It should be a positive numeric value (currently the script has only been tested with integers).



* *dist\_range*: Controls the 'maximum visible distances' checked for visibility by the script. It should be a list of numeric values (currently the script has only been tested with integers).



* *densify\_dist*: Controls the density at which individual vectors in the radial array will sample the input DEM. Lower values will produce better results, at the cost of increased computation. It should be a list of numeric values (currently the script has only been tested with integers).



* *obs\_z\_range*: Controls observer heights used by the script. It should be a list of numeric values (currently the script has only been tested with integers).



#### 1.3 *vis\_netcdf* Optional Script Inputs



* *array\_angular\_increment*: Controls the angular resolution of radial arrays produced by the script. Measured in degrees of arc. The default value is 1.



* *obs\_z\_offset*: Adds an offset (in vertical units) to the observer. This is meant to represent the height of an observer above the land/sea surface. It should be a numeric value. By default, it is set to zero.


* *output\_null\_value*: Assigns a numeric value for null results (i.e. summary values could not be calculated for an array at that observer point). By default, the value is zero.


* *sample\_raster*: A raster containing values which will be 'sampled' by visible points determined by the script. This should be a valid path to a raster file. This parameter is an empty string by default - meaning no raster values will be calculated.



* *landmark\_fc*: Features containing 'landmarks'. If a valid path is supplied, then the script will calculate how many landmarks can be seen from each observer point / elevation combination. This will be added as a separate variable in the output NetCDF.



* *override\_dist\_list*: A list containing maximum visible distances to be used by the script. This list overrides the *dist\_range* parameter. Use this if you want to specify distance ranges that are not evenly spaced (i.e. \[500, 2500, 5000] instead of \[500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]).



* *pt\_mask\_json*: A valid path with a .json file extension. If such a path is present, the script will write the results of grid / point mask intersection as a dictionary to the .json file, and will load that file when repeating the analysis. This was implemented to save processing time on large grids. Note that the script will NOT recompute the .json file if other input parameters are changed - changing the size or resolution of the grid after creating a pt\_mask\_json will probably crash the script.



* *write\_log*: A boolean controlling whether text logs are written by the script. Lots are written as .json files, in the same directory as the output netcdf. By default this parameter is *False*.



#### 1.4 *vis\_netcdf* Demo



Currently two demos have been prepared in the directory '\\radial-array-tools\\demo\\vis\_netcdf'. Instructions and further description can be found in a separate README in this directory.





## Bibliography



Lonergan, C., Hedley, N., 2016. *Unpacking isovists: a framework for 3D spatial visibility analysis*. Cartography and Geographic Information Science 43, 87–102.



Pouncett, J., 2013. *Expanding Horizons: Visibility, Monuments and Topography*, in: *Across Space and Time. Papers from the 41st Annual Conference of Computer Applications and Quantitative Methods in Archaeology (CAA), Perth, 25-28 March 2013*. Presented at the CAA 2013, Perth.



Rew, R., Davis, G., 1990. *NetCDF: an interface for scientific data access*. IEEE Computer Graphics and Applications 10, 76–82.



Smith, K., Hulin, L. *Mediterranean Maritime Visibility: Old Limits and New Approaches*. Antiquity. Forthcoming; accepted for publication.



Wheatley, D., Gillings, M., 2000. *Vision, Perception and GIS,* in: Lock, G. (Ed.), *Beyond the Map: Archaeology and Spatial Technologies*. IOS Press, Amsterdam, pp. 1–27.



