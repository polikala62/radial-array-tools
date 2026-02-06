VIS_NETCDF_DEMO README

Last updated February 04, 2026


This folder contains two demos, one demonstrating the operation of the script over a large area (45 km), and the other showing it operating at higher resolutions (c. 6 km). Demo 1 is set in the San Juan Islands, an archipelago of 200-odd islands in the northwestern United States. Demo 2 is set around Waldron Island, which is part of the San Juans. The DEM used in this demo is derived from the USGS GTOPO30 satellite DEM - it has been reprojected into UTM Zone 10 North and bilinearly interpolated to a resolution of 30m.

Follow the following steps to complete the demo:

1. Open 'vis_netcdf_demo1_parameters.json' or 'vis_netcdf_demo2_parameters.json', and check that filepaths in the .json files are correct. All paths should be to files in subfolders within this directory. Remember to use '\\' rather than '\' for folder separators (otherwise Python might mistake parts of the path for string literals).

2. Run the visibility netcdf script at 'source\vis_netcdf_run.py'. You will be prompted to input a parameter file path. Use the path to either 'vis_netcdf_demo1_parameters.json' or 'vis_netcdf_demo2_parameters.json' to start the script. You will be prompted to overwrite the output, enter 'y' and press ENTER to do so.

3. Wait for the script to finish. At the time of writing, demo 1 took about 45 minutes to complete, and demo 2 took about 4.5 hours.

4. Find the output in the 'output' folder in this directory. The output for demo 1 should be 'vis_netcdf_demo1_output.nc', and the output for demo 2 should be 'vis_netcdf_demo2_output.nc'. NetCDF is a multidimensional data storage format. At the time of writing it could be opened as a raster in QGIS, and as a multidimensional raster in ESRI ArcPro.

See the main README file in this repository for descriptions of the output.