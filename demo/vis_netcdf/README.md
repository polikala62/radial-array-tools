\# VIS\_NETCDF\_DEMO README

Last updated February 2026



\## Overview



This folder contains two demos set in the San Juan Islands, an archipelago of 200-odd islands in the northwestern United States. If you're wondering why, it's because I grew up there and had the datasets handy. It's useful to have a testing area where I already know what visibility conditions should like - at least for a few areas!



Demo 1 demonstrates the operation of the 'vis\_netcdf' script over the entire archipelago at a low resolution (1km). It also samples a set of landmarks representing topographic prominences (see below).



Demo 2 demonstrates the operation of the 'vis\_netcdf' script over a small area at a higher resolution (30m). It focusses on Waldron Island, in the NW corner of the San Juans.



Demo 3 is still a work in progress - it will demonstrate the 'raster sampling' aspect of the script. At the moment raster sampling is very processing-intensive and could not be reasonably incorporated into a demo. I am working to make raster sampling quicker, and when it has been improved I hope to add this demo.



\## Data Sources



The DEM used in this demo is derived from the USGS GTOPO30 satellite DEM - it has been reprojected into UTM Zone 10 North and bilinearly interpolated to a resolution of 30m. The 'vis\_netcdf\_demo1\_landmarks.shp' shapefile contains 'peaks' identified using the ArcGIS Pro tool 'Geomorphon Landforms'. This dataset is a derivative of the reprojected DEM - features have been filtered for size and names have been added.



\## Directions



Follow the following steps to complete the demo:



1. Open 'vis\_netcdf\_demo1\_parameters.json' or 'vis\_netcdf\_demo2\_parameters.json', and check that filepaths in the .json files are correct. All paths should be to files in subfolders within this directory. Remember to use '\\' rather than '' for folder separators (otherwise Python might mistake parts of the path for string literals).
2. Run the visibility netcdf script at 'source\\vis\_netcdf\_run.py'. You will be prompted to input a parameter file path. Use the path to either 'vis\_netcdf\_demo1\_parameters.json' or 'vis\_netcdf\_demo2\_parameters.json' to start the script. You will be prompted to overwrite the output, enter 'y' and press ENTER to do so.
3. Wait for the script to finish. At the time of writing, demo 1 took about 45 minutes to complete, and demo 2 took about 4.5 hours.
4. Find the output in the 'output' folder in this directory. The output for demo 1 should be 'vis\_netcdf\_demo1\_output.nc', and the output for demo 2 should be 'vis\_netcdf\_demo2\_output.nc'. NetCDF is a multidimensional data storage format. At the time of writing it could be opened as a raster in QGIS, and as a multidimensional raster in ESRI ArcPro.



See the main README file in this repository for descriptions of the output.

